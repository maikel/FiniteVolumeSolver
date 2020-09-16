// Copyright (c) 2019 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "fub/AMReX.hpp"
#include "fub/AMReX_CutCell.hpp"
#include "fub/Solver.hpp"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Union.H>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

#include <cmath>
#include <iostream>

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 3;

static constexpr int scratch_gcw = 8;
static constexpr int flux_gcw = scratch_gcw - 2;

static constexpr double r_tube = 0.015;
static constexpr double r_inner = 0.5 * 0.130;
static constexpr double r_outer = 0.5 * 0.389;
static constexpr double r_tube_center = 2.0 * r_inner;
static constexpr double alpha = 2. * M_PI / 6.;

auto Center(double x, double phi) -> ::amrex::RealArray {
  using std::cos;
  using std::sin;
  return {x, r_tube_center * cos(phi), r_tube_center * sin(phi)};
}

auto MakeTubeSolver(fub::Burke2012& mechanism,
                    const fub::ProgramOptions& options, int k,
                    const std::shared_ptr<fub::CounterRegistry>& counters) {
  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << fmt::format(
      "==================== Tube #{} =========================", k);

  using namespace fub::amrex;
  std::vector<pybind11::dict> dicts{};
  dicts = fub::GetOptionOr(options, "Tubes", dicts);
  if (dicts.size() < size_t(k)) {
    throw std::runtime_error("You need to specify options for each tube.");
  }
  fub::ProgramOptions tube_options = fub::ToMap(dicts[k]);

  CartesianGridGeometry grid_geometry(
      fub::GetOptions(tube_options, "GridGeometry"));
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  PatchHierarchyOptions hierarchy_options(
      fub::GetOptions(tube_options, "PatchHierarchy"));
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  using Complete = fub::IdealGasMix<Tube_Rank>::Complete;
  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 1e-3),
                            std::make_pair(&Complete::pressure, 1e-2),
                            std::make_pair(&Complete::temperature, 1e-1)};

  ::amrex::Box refine_box{{grid_geometry.cell_dimensions[0] - 5, 0, 0},
                          {grid_geometry.cell_dimensions[0] - 1, 0, 0}};
  ConstantBox constant_box{refine_box};

  Complete state(equation);
  {
    using namespace std::literals;
    const fub::ProgramOptions initial_options =
        fub::GetOptions(tube_options, "InitialCondition");
    const std::string moles =
        fub::GetOptionOr(initial_options, "moles", "N2:79,O2:21"s);
    const double temperature =
        fub::GetOptionOr(initial_options, "temperature", 300.0);
    const double pressure =
        fub::GetOptionOr(initial_options, "pressure", 101325.0);
    const double velocity = fub::GetOptionOr(initial_options, "velocity", 0.0);
    equation.GetReactor().SetMoleFractions(moles);
    equation.GetReactor().SetTemperature(temperature);
    equation.GetReactor().SetPressure(pressure);
    equation.CompleteFromReactor(state, fub::Array<double, 1, 1>(velocity));
    equation.CompleteFromCons(state, state);
  }
  ConstantData initial_data{equation, state};

  PressureValveOptions valve_options =
      fub::GetOptions(tube_options, "PressureValveBoundary");
  BOOST_LOG(log) << "PressureValveBoundary:";
  valve_options.Print(log);

  PressureValveBoundary valve{equation, valve_options};
  IsentropicPressureBoundary right_boundary(equation, 101325.0,
                                            fub::Direction::X, 1);
  BoundarySet boundaries{{valve, right_boundary}};

  // If a checkpoint path is specified we will fill the patch hierarchy with
  // data from the checkpoint file, otherwise we will initialize the data by
  // the initial data function.
  std::shared_ptr<GriddingAlgorithm> gridding = [&] {
    std::string checkpoint{};
    if (options.count("checkpoint")) {
      checkpoint = options.at("checkpoint").cast<std::string>();
    }
    if (!checkpoint.empty()) {
      MPI_Comm comm = MPI_COMM_WORLD;
      std::string input = fub::ReadAndBroadcastFile(
          fmt::format("{}/Valve_{}", checkpoint, k), comm);
      auto ifs = std::istringstream(input);
      boost::archive::text_iarchive ia(ifs);
      ia >> *(valve.GetSharedState());
    }
    if (checkpoint.empty()) {
      BOOST_LOG(log) << "Initialize grid by initial condition...";
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(equation, grid_geometry, hierarchy_options),
          initial_data, TagAllOf(gradient, constant_box), boundaries);
      gridding->GetPatchHierarchy().SetCounterRegistry(counters);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint = fmt::format("{}/Tube_{}", checkpoint, k);
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
      fub::amrex::DataDescription desc =
          fub::amrex::MakeDataDescription(equation);
      PatchHierarchy h = ReadCheckpointFile(checkpoint, desc, grid_geometry,
                                            hierarchy_options);
      h.SetCounterRegistry(counters);
      std::shared_ptr<GriddingAlgorithm> gridding =
          std::make_shared<GriddingAlgorithm>(std::move(h), initial_data,
                                              TagAllOf(gradient, constant_box),
                                              boundaries);
      return gridding;
    }
  }();

  fub::ideal_gas::MusclHancockPrimMethod<Tube_Rank> flux_method{equation};
  HyperbolicMethod method{FluxMethodAdapter(flux_method),
                          EulerForwardTimeIntegrator(),
                          Reconstruction(equation)};

  IntegratorContext context(gridding, method, scratch_gcw, flux_gcw);

  BOOST_LOG(log) << fmt::format(
      "==================== End Tube #{} =========================", k);

  return std::pair{std::move(context), valve};
}

auto MakePlenumSolver(fub::Burke2012& mechanism,
                      const fub::ProgramOptions& options) {
  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << "==================== Plenum =========================";
  const fub::ProgramOptions plenum_options = fub::GetOptions(options, "Plenum");

  fub::amrex::CartesianGridGeometry grid_geometry(
      fub::GetOptions(plenum_options, "GridGeometry"));
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  auto MakePolygon = [](auto... points) {
    fub::Polygon::Vector xs{};
    fub::Polygon::Vector ys{};

    (xs.push_back(std::get<0>(points)), ...);
    (ys.push_back(std::get<1>(points)), ...);

    return fub::Polygon(std::move(xs), std::move(ys));
  };

  double r_inlet_start = r_tube;
  double r_inlet_end = 0.0225;
  {
    const fub::ProgramOptions inlet_options =
        fub::GetOptions(plenum_options, "InletGeometry");
    r_inlet_start =
        fub::GetOptionOr(inlet_options, "r_start", r_inlet_start);
    r_inlet_end = fub::GetOptionOr(inlet_options, "r_end", r_inlet_end);
  }

  auto DivergentInlet = [&](double height,
                            const std::array<double, 3>& center) {
    const double xlo = center[0] - height;
    const double xhi = center[0];
    const double xdiv = xhi - 0.075;
    const double r = r_inlet_start;
    const double r2 = r_inlet_end;
    auto polygon =
        MakePolygon(std::pair{xlo, r}, std::pair{xdiv, r}, std::pair{xhi, r2},
                    std::pair{xhi, -r2}, std::pair{xdiv, -r},
                    std::pair{xlo, -r}, std::pair{xlo, r});
    auto tube_in_zero = fub::Invert(fub::RotateAxis(polygon));
    amrex::RealArray real_center{center[0], center[1], center[2]};
    auto tube_in_center = amrex::EB2::TranslationIF(
        fub::amrex::Geometry(tube_in_zero), real_center);
    return tube_in_center;
  };

  auto Cylinder = [&](double radius, double height,
                      const std::array<double, 3>& center) {
    return amrex::EB2::CylinderIF(radius, height, 0, center, true);
  };

  auto embedded_boundary = amrex::EB2::makeIntersection(
      amrex::EB2::makeUnion(
          amrex::EB2::makeIntersection(
              Cylinder(r_outer, 0.5, {0.25, 0.0, 0.0}),
              DivergentInlet(0.2, Center(0.0, 0.0 * alpha)),
              DivergentInlet(0.2, Center(0.0, 1.0 * alpha)),
              DivergentInlet(0.2, Center(0.0, 2.0 * alpha)),
              DivergentInlet(0.2, Center(0.0, 3.0 * alpha)),
              DivergentInlet(0.2, Center(0.0, 4.0 * alpha)),
              DivergentInlet(0.2, Center(0.0, 5.0 * alpha))),
          amrex::EB2::CylinderIF(r_inner, 1.0, 0, {0.25, 0.0, 0.0}, false)),
      amrex::EB2::PlaneIF({0.5, 0.0, 0.0}, {1.0, 0.0, 0.0}, false));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};
  fub::Complete<fub::IdealGasMix<Plenum_Rank>> state(equation);
  {
    using namespace std::literals;
    const fub::ProgramOptions initial_options =
        fub::GetOptions(plenum_options, "InitialCondition");
    const std::string moles =
        fub::GetOptionOr(initial_options, "moles", "N2:79,O2:21"s);
    const double temperature =
        fub::GetOptionOr(initial_options, "temperature", 300.0);
    const double pressure =
        fub::GetOptionOr(initial_options, "pressure", 101325.0);
    std::array<double, 3> velocity{};
    velocity = fub::GetOptionOr(initial_options, "velocity", velocity);
    equation.GetReactor().SetMoleFractions(moles);
    equation.GetReactor().SetTemperature(temperature);
    equation.GetReactor().SetPressure(pressure);
    equation.CompleteFromReactor(state,
                                 {velocity[0], velocity[1], velocity[2]});
    equation.CompleteFromCons(state, state);
  }
  using namespace fub::amrex::cutcell;
  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.0),
                              state, state);

  PatchHierarchyOptions hierarchy_options(
      fub::GetOptions(plenum_options, "PatchHierarchy"));
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  ::amrex::Geometry coarse_geometry =
      fub::amrex::GetCoarseGeometry(grid_geometry);

  BOOST_LOG(log) << "Compute EB level set data...";
  std::shared_ptr<fub::CounterRegistry> registry =
      std::make_shared<fub::CounterRegistry>();
  {
    fub::Timer timer = registry->get_timer("MakeIndexSpaces");
    hierarchy_options.index_spaces = MakeIndexSpaces(
        shop, coarse_geometry, hierarchy_options.max_number_of_levels);
  }

  using State = fub::Complete<fub::IdealGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.01}};

  ::amrex::RealBox inlet{{-0.1, -0.5, -0.5}, {0.05, +0.5, +0.5}};
  const ::amrex::Box refine_box =
      fub::amrex::BoxWhichContains(inlet, coarse_geometry);
  ConstantBox constant_box{refine_box};

  // ::amrex::RealBox outlet{{0.5, -0.5, -0.5}, {0.54, +0.5, +0.5}};
  // const ::amrex::Box outlet_box =
  //     fub::amrex::BoxWhichContains(outlet, coarse_geometry);

  IsentropicPressureBoundaryOptions boundary_options =
      fub::GetOptions(plenum_options, "IsentropicPressureBoundary");
  BOOST_LOG(log) << "IsentropicPressureBoundary:";
  boundary_options.Print(log);

  BoundarySet boundary_condition{
      {TransmissiveBoundary{fub::Direction::X, 0},
       IsentropicPressureBoundary{equation, boundary_options},
       TransmissiveBoundary{fub::Direction::Y, 0},
       TransmissiveBoundary{fub::Direction::Y, 1},
       TransmissiveBoundary{fub::Direction::Z, 0},
       TransmissiveBoundary{fub::Direction::Z, 1}}};

  // If a checkpoint path is specified we will fill the patch hierarchy with
  // data from the checkpoint file, otherwise we will initialize the data by
  // the initial data function.
  std::shared_ptr gridding = [&] {
    std::string checkpoint{};
    if (options.count("checkpoint")) {
      checkpoint = options.at("checkpoint").cast<std::string>();
    }
    if (checkpoint.empty()) {
      BOOST_LOG(log) << "Initialize grid by initial condition...";
      PatchHierarchy h(equation, grid_geometry, hierarchy_options);
      h.SetCounterRegistry(registry);
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
          boundary_condition);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint += "/Plenum";
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), grid_geometry,
          hierarchy_options);
      h.SetCounterRegistry(registry);
      return std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
          boundary_condition);
    }
  }();

  // Make Solver

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::ideal_gas::MusclHancockPrimMethod<Plenum_Rank> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method{flux_method, hll_method};

  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  IntegratorContext context(gridding, method, scratch_gcw, flux_gcw);

  BOOST_LOG(log) << "==================== End Plenum =========================";
  return context;
}

void WriteCheckpoint(
    const std::string& path,
    const fub::amrex::MultiBlockGriddingAlgorithm& grid,
    fub::span<const std::shared_ptr<fub::amrex::PressureValve>> valves,
    int rank, const fub::amrex::MultiBlockIgniteDetonation& ignition) {
  auto tubes = grid.GetTubes();
  std::ptrdiff_t cycles = grid.GetCycles();
  int k = 0;
  for (auto& tube : tubes) {
    std::string name = fmt::format("{}/{:09}/Tube_{}", path, cycles, k);
    fub::amrex::WriteCheckpointFile(name, tube->GetPatchHierarchy());

    if (rank == 0) {
      std::string valve = fmt::format("{}/{:09}/Valve_{}", path, cycles, k);
      std::ofstream valve_checkpoint(valve);
      boost::archive::text_oarchive oa(valve_checkpoint);
      oa << *valves[k];
    }

    k = k + 1;
  }
  std::string name = fmt::format("{}/{:09}/Plenum", path, cycles);
  fub::amrex::cutcell::WriteCheckpointFile(
      name, grid.GetPlena()[0]->GetPatchHierarchy());
  if (rank == 0) {
    name = fmt::format("{}/{:09}/Ignition", path, cycles);
    std::ofstream ignition_checkpoint(name);
    boost::archive::text_oarchive oa(ignition_checkpoint);
    oa << ignition.GetNextIgnitionTimePoints();
  }
}

struct CheckpointOutput : fub::OutputAtFrequencyOrInterval<
                              fub::amrex::MultiBlockGriddingAlgorithm> {
  CheckpointOutput(
      const fub::ProgramOptions& options,
      std::vector<std::shared_ptr<fub::amrex::PressureValve>> valves,
      const fub::amrex::MultiBlockIgniteDetonation* ignition)
      : OutputAtFrequencyOrInterval(options), valves_{std::move(valves)},
        ignition_{ignition} {
    directory_ = fub::GetOptionOr(options, "directory", directory_);
    fub::SeverityLogger log = fub::GetInfoLogger();
    BOOST_LOG(log) << "CheckpointOutput:";
    OutputAtFrequencyOrInterval::Print(log);
    BOOST_LOG(log) << fmt::format(" - directory = '{}'", directory_);
  }

  void
  operator()(const fub::amrex::MultiBlockGriddingAlgorithm& grid) override {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    boost::log::sources::severity_logger<boost::log::trivial::severity_level>
        log(boost::log::keywords::severity = boost::log::trivial::info);
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", grid.GetTimePoint().count());
    BOOST_LOG(log) << fmt::format("Write checkpoint to '{}'.", directory_);
    WriteCheckpoint(directory_, grid, valves_, rank, *ignition_);
  }

  std::vector<std::shared_ptr<fub::amrex::PressureValve>> valves_{};
  const fub::amrex::MultiBlockIgniteDetonation* ignition_{};
  std::string directory_{"./Checkpoint/"};
};

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard scope_guard{};

  fub::Burke2012 mechanism{};

  fub::SeverityLogger log = fub::GetInfoLogger();

  std::vector<fub::amrex::cutcell::IntegratorContext> plenum{};
  std::vector<fub::amrex::IntegratorContext> tubes{};
  std::vector<fub::amrex::BlockConnection> connectivity{};
  std::vector<std::shared_ptr<fub::amrex::PressureValve>> valves{};

  BOOST_LOG(log) << "Make Plenum Solver...";
  plenum.push_back(MakePlenumSolver(mechanism, options));
  auto counter_database = plenum[0].GetCounterRegistry();

  auto MakeConnection = [&](int k) {
    BOOST_LOG(log) << "Make Tube Solver " << k << "...";
    auto&& [tube, valve] =
        MakeTubeSolver(mechanism, options, k, counter_database);
    tubes.push_back(std::move(tube));
    valves.push_back(valve.GetSharedState());
    fub::amrex::BlockConnection connection;
    connection.direction = fub::Direction::X;
    connection.side = 0;
    connection.ghost_cell_width = scratch_gcw;
    connection.plenum.id = 0;
    connection.tube.id = k;
    connection.tube.mirror_box = tubes[k].GetGeometry(0).Domain();
    const double xlo = plenum[0].GetGeometry(0).ProbLo(0);
    connection.plenum.mirror_box = fub::amrex::BoxWhichContains(
        fub::amrex::DomainAroundCenter(Center(xlo, k * alpha),
                                       {0.03, r_tube, r_tube}),
        plenum[0].GetGeometry(0));
    BOOST_LOG(log) << "Mirror Box = " << connection.plenum.mirror_box;
    return connection;
  };

  connectivity.push_back(MakeConnection(0));
  connectivity.push_back(MakeConnection(1));
  connectivity.push_back(MakeConnection(2));
  connectivity.push_back(MakeConnection(3));
  connectivity.push_back(MakeConnection(4));
  connectivity.push_back(MakeConnection(5));

  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};
  fub::IdealGasMix<Plenum_Rank> plenum_equation{mechanism};

  fub::amrex::MultiBlockIntegratorContext context(
      fub::FlameMasterReactor(mechanism), std::move(tubes), std::move(plenum),
      std::move(connectivity));

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<Plenum_Rank>, std::move(context), fub::StrangSplitting{});
      // fub::int_c<Plenum_Rank>, std::move(context), fub::GodunovSplitting{});

  const std::size_t n_tubes = system_solver.GetContext().Tubes().size();
  const int max_number_of_levels = system_solver.GetContext()
                                       .Tubes()[0]
                                       .GetPatchHierarchy()
                                       .GetMaxNumberOfLevels();

  fub::amrex::IgniteDetonationOptions ignite_options(
      fub::GetOptions(options, "IgniteDetonation"));
  BOOST_LOG(log) << "IgniteDetonation:";
  ignite_options.Print(log);
  fub::amrex::MultiBlockIgniteDetonation ignition{
      tube_equation, n_tubes, max_number_of_levels, ignite_options};

  std::string checkpoint{};
  if (options.count("checkpoint")) {
    checkpoint = options.at("checkpoint").cast<std::string>();
  }
  if (!checkpoint.empty()) {
    MPI_Comm comm = context.GetMpiCommunicator();
    std::string input =
        fub::ReadAndBroadcastFile(checkpoint + "/Ignition", comm);
    std::istringstream ifs(input);
    boost::archive::text_iarchive ia(ifs);
    std::vector<fub::Duration> last_ignitions;
    ia >> last_ignitions;
    ignition.SetNextIgnitionTimePoints(last_ignitions);
  }

  fub::SplitSystemSourceLevelIntegrator ign_solver(
      std::move(system_solver), std::move(ignition), fub::GodunovSplitting{});

  fub::amrex::MultiBlockKineticSouceTerm source_term(tube_equation);

  fub::SplitSystemSourceLevelIntegrator level_integrator(
      std::move(ign_solver), std::move(source_term),
      fub::StrangSplittingLumped{});

  // fub::NoSubcycleSolver solver(std::move(level_integrator));
  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  fub::OutputFactory<fub::amrex::MultiBlockGriddingAlgorithm> factory{};
  factory.RegisterOutput<fub::amrex::MultiWriteHdf5>("HDF5");
  factory.RegisterOutput<fub::amrex::MultiBlockPlotfileOutput>("Plotfile");
  factory.RegisterOutput<fub::amrex::LogProbesOutput>("LogProbes");
  factory.RegisterOutput<CheckpointOutput>(
      "Checkpoint", valves,
      &solver.GetLevelIntegrator().GetSystem().GetSource());
  using CounterOutput =
      fub::CounterOutput<fub::amrex::MultiBlockGriddingAlgorithm,
                         std::chrono::milliseconds>;
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  fub::MultipleOutputs<fub::amrex::MultiBlockGriddingAlgorithm> outputs(
      std::move(factory), fub::GetOptions(options, "Output"));

  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);

  outputs(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
}

int main(int argc, char** argv) {
  // fub::EnableFloatingPointExceptions();
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
  pybind11::scoped_interpreter interpreter{};
  std::optional<fub::ProgramOptions> opts = fub::ParseCommandLine(argc, argv);
  if (opts) {
    MyMain(*opts);
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}
