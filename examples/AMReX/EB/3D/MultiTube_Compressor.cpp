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
#include <AMReX_EB2_IF_Union.H>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>

#include <cmath>
#include <iostream>

#include <hdf5.h>

#include <xmmintrin.h>

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 3;

static constexpr double r_tube = 0.015;
static constexpr double r_inner = 0.5 * 0.130;
static constexpr double r_outer = 0.5 * 0.385;
static constexpr double r_tube_center = 0.5 * r_inner + 0.5 * r_outer;
static constexpr double alpha = 2. * M_PI / 5.;

auto Center(double x, double phi) -> ::amrex::RealArray {
  using std::cos;
  using std::sin;
  return {x, r_tube_center * sin(phi), r_tube_center * cos(phi)};
}

struct NoInit {
  static void InitializeData(fub::amrex::PatchLevel&,
                             const fub::amrex::GriddingAlgorithm& /*grid*/,
                             int /*level*/, fub::Duration /*time*/) noexcept {}
};

auto MakeTubeSolver(fub::Burke2012& mechanism,
                    const fub::ProgramOptions& options, int k) {
  using namespace fub::amrex;
  std::vector<pybind11::dict> dicts{};
  dicts = fub::GetOptionOr(options, "Tubes", dicts);
  if (dicts.size() < size_t(k)) {
    throw std::runtime_error("You need to specify options for each tube.");
  }
  fub::ProgramOptions tube_options = fub::ToMap(dicts[k]);
  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};

  CartesianGridGeometry grid_geometry(
      fub::GetOptions(tube_options, "GridGeometry"));
  PatchHierarchyOptions hierarchy_options(
      fub::GetOptions(tube_options, "PatchHierarchy"));
  DataDescription desc = MakeDataDescription(equation);

  using Complete = fub::IdealGasMix<Tube_Rank>::Complete;
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
    equation.GetReactor().SetMoleFractions(moles);
    equation.GetReactor().SetTemperature(temperature);
    equation.GetReactor().SetPressure(pressure);
    equation.CompleteFromReactor(state);
  }
  ConstantData initial_data{equation, state};

  PressureValveOptions valve_options =
      fub::GetOptions(tube_options, "PressureValveBoundary");
  PressureValveBoundary valve{equation, valve_options};
  IsentropicPressureBoundary right{equation, 101325.0, fub::Direction::X, 1};
  BoundarySet boundaries{{valve, right}};

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
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(desc, grid_geometry, hierarchy_options), initial_data,
          TagAllOf(gradient, constant_box), boundaries);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint = fmt::format("{}/Tube_{}", checkpoint, k);
      PatchHierarchy h = ReadCheckpointFile(checkpoint, desc, grid_geometry,
                                            hierarchy_options);
      std::shared_ptr<GriddingAlgorithm> gridding =
          std::make_shared<GriddingAlgorithm>(std::move(h), NoInit{},
                                              TagAllOf(gradient, constant_box),
                                              boundaries);
      return gridding;
    }
  }();

  fub::ideal_gas::MusclHancockPrimMethod<Tube_Rank> flux_method{equation};
  HyperbolicMethod method{FluxMethodAdapter(flux_method), EulerForwardTimeIntegrator(),
                          Reconstruction(equation)};

  const int scratch_gcw = 4;
  const int flux_gcw = 2;

  IntegratorContext context(gridding, method, scratch_gcw, flux_gcw);

  return std::pair{std::move(context), valve};
}

auto MakePlenumSolver(fub::Burke2012& mechanism,
                      const fub::ProgramOptions& options) {
  using namespace fub::amrex::cutcell;
  auto Cylinder = [&](double radius, double height,
                      const std::array<double, 3>& center) {
    return amrex::EB2::CylinderIF(radius, height, 0, center, true);
  };

  auto embedded_boundary = amrex::EB2::makeUnion(
      amrex::EB2::makeIntersection(
          Cylinder(r_outer, 0.5, {0.25, 0.0, 0.0}),
          Cylinder(r_tube, 0.3, Center(0.6, 0.0 * alpha)),
          Cylinder(r_tube, 0.3, Center(0.6, 1.0 * alpha)),
          Cylinder(r_tube, 0.3, Center(0.6, 2.0 * alpha)),
          Cylinder(r_tube, 0.3, Center(0.6, 3.0 * alpha)),
          Cylinder(r_tube, 0.3, Center(0.6, 4.0 * alpha))),
      amrex::EB2::CylinderIF(r_inner, 1.0, 0, {0.25, 0.0, 0.0}, false));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  const fub::ProgramOptions plenum_options = fub::GetOptions(options, "Plenum");

  fub::amrex::CartesianGridGeometry grid_geometry(
      fub::GetOptions(plenum_options, "GridGeometry"));

  amrex::Geometry coarse_geom = fub::amrex::GetCoarseGeometry(grid_geometry);

  PatchHierarchyOptions hierarchy_options =
      fub::GetOptions(plenum_options, "PatchHierarchy");
  hierarchy_options.index_spaces =
      MakeIndexSpaces(shop, grid_geometry, hierarchy_options);

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};
  using Complete = fub::Complete<fub::IdealGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&Complete::pressure, 0.05},
                             std::pair{&Complete::density, 0.01}};

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
    equation.GetReactor().SetMoleFractions(moles);
    equation.GetReactor().SetTemperature(temperature);
    equation.GetReactor().SetPressure(pressure);
    equation.CompleteFromReactor(state);
  }

  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04),
                              state, state);

  ::amrex::RealBox inlet{{-0.1, -0.5, -0.5}, {0.05, +0.5, +0.5}};
  ::amrex::RealBox outlet{{0.5, -0.5, -0.5}, {0.54, +0.5, +0.5}};
  //const ::amrex::Box inlet_box =
  //    fub::amrex::BoxWhichContains(inlet, coarse_geom);
  const ::amrex::Box outlet_box =
      fub::amrex::BoxWhichContains(outlet, coarse_geom);
  ConstantBox constant_box{outlet_box};

  MassflowBoundaryOptions massflow_options =
      fub::GetOptions(plenum_options, "MassflowBoundary");
  BoundarySet boundary_condition{{MassflowBoundary{equation, massflow_options},
                                  TransmissiveBoundary{fub::Direction::X, 1}}};

  // If a checkpoint path is specified we will fill the patch hierarchy with
  // data from the checkpoint file, otherwise we will initialize the data by
  // the initial data function.
  std::shared_ptr gridding = [&] {
    std::string checkpoint{};
    if (options.count("checkpoint")) {
      checkpoint = options.at("checkpoint").cast<std::string>();
    }
    if (checkpoint.empty()) {
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(equation, grid_geometry, hierarchy_options),
          initial_data,
          TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
          boundary_condition);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint += "/Plenum";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), grid_geometry,
          hierarchy_options);
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
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  const int scratch_gcw = 4;
  const int flux_gcw = 2;

  return IntegratorContext(gridding, method, scratch_gcw, flux_gcw);
}

struct CheckpointOptions {
  CheckpointOptions() = default;

  CheckpointOptions(const fub::ProgramOptions& options) {
    checkpoint = fub::GetOptionOr(options, "checkpoint", checkpoint);
  }

  template <typename Logger> void Print(Logger& log) const {
    if (!checkpoint.empty()) {
      BOOST_LOG(log) << "Restart simulation from checkpoint '" << checkpoint
                     << "'!";
    } else {
      BOOST_LOG(log) << "No Checkpoint given.";
    }
  }

  std::string checkpoint{};
};

void WriteCheckpoint(
    const std::string& path,
    const fub::amrex::MultiBlockGriddingAlgorithm& grid,
    fub::span<const std::shared_ptr<fub::amrex::PressureValve>> valves,
    int rank, const fub::amrex::MultiBlockIgniteDetonation& ignition) {
  auto tubes = grid.GetTubes();
  int k = 0;
  for (auto& tube : tubes) {
    std::string name = fmt::format("{}/Tube_{}", path, k);
    fub::amrex::WriteCheckpointFile(name, tube->GetPatchHierarchy());

    if (rank == 0) {
      std::string valve = fmt::format("{}/Valve_{}", path, k);
      std::ofstream valve_checkpoint(valve);
      boost::archive::text_oarchive oa(valve_checkpoint);
      oa << *valves[k];
    }

    k = k + 1;
  }
  std::string name = fmt::format("{}/Plenum", path);
  fub::amrex::cutcell::WriteCheckpointFile(
      name, grid.GetPlena()[0]->GetPatchHierarchy());
  if (rank == 0) {
    name = fmt::format("{}/Ignition", path);
    std::ofstream ignition_checkpoint(name);
    boost::archive::text_oarchive oa(ignition_checkpoint);
    oa << ignition.GetNextIgnitionTimePoints();
  }
}

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard scope_guard{};
  fub::Burke2012 mechanism{};

  std::vector<fub::amrex::cutcell::IntegratorContext> plenum{};
  std::vector<fub::amrex::IntegratorContext> tubes{};
  std::vector<fub::amrex::BlockConnection> connectivity{};
  std::vector<std::shared_ptr<fub::amrex::PressureValve>> valves{};

  plenum.push_back(MakePlenumSolver(mechanism, options));

  auto MakeConnection = [&](int k) {
    auto&& [tube, valve] = MakeTubeSolver(mechanism, options, k);
    fub::amrex::BlockConnection connection;
    connection.direction = fub::Direction::X;
    connection.side = 1;
    connection.ghost_cell_width = 4;
    connection.plenum.id = 0;
    connection.tube.id = k;
    connection.tube.mirror_box = tubes[k]
                                     .GetGriddingAlgorithm()
                                     ->GetPatchHierarchy()
                                     .GetGeometry(0)
                                     .Domain();
    connection.plenum.mirror_box = fub::amrex::BoxWhichContains(
        fub::amrex::DomainAroundCenter(Center(+0.53, k * alpha),
                                       {0.03, r_tube, r_tube}),
        plenum[0].GetGeometry(0));

    tubes.push_back(std::move(tube));
    valves.push_back(valve.GetSharedState());
    connectivity.push_back(connection);
  };

  MakeConnection(0);
  MakeConnection(1);
  MakeConnection(2);
  MakeConnection(3);
  MakeConnection(4);

  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};

  fub::amrex::MultiBlockIntegratorContext context(
      fub::FlameMasterReactor(mechanism), std::move(tubes), std::move(plenum),
      std::move(connectivity), valves);

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<Plenum_Rank>, std::move(context), fub::GodunovSplitting{});

  const std::size_t n_tubes = system_solver.GetContext().Tubes().size();
  const int max_number_of_levels = system_solver.GetContext()
                                       .Tubes()[0]
                                       .GetPatchHierarchy()
                                       .GetMaxNumberOfLevels();

  fub::amrex::MultiBlockIgniteDetonation ignition{
      tube_equation, n_tubes, max_number_of_levels,
      fub::amrex::IgniteDetonationOptions(options, "IgniteDetonation")};

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
      std::move(ign_solver), std::move(source_term), fub::StrangSplitting{});

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  fub::OutputFactory<fub::amrex::MultiBlockGriddingAlgorithm> factory{};
  factory.RegisterOutput<fub::amrex::MultiWriteHdf5>("HDF5");
  factory.RegisterOutput<fub::amrex::MultiBlockPlotfileOutput>("Plotfile");
  factory.RegisterOutput<fub::amrex::LogProbesOutput>("LogProbes");
  using CounterOutput =
      fub::CounterOutput<fub::amrex::MultiBlockGriddingAlgorithm,
                         std::chrono::milliseconds>;
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  fub::MultipleOutputs<fub::amrex::MultiBlockGriddingAlgorithm> outputs(
      std::move(factory), fub::GetOptions(options, "Output"));

  outputs(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, fub::GetOptions(options, "RunOptions"),
                     wall_time_reference, outputs);
}

int main(int argc, char** argv) {
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
