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
#include <AMReX_EB_LSCore.H>

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

auto DomainAroundCenter(const ::amrex::RealArray& x, double rx)
    -> ::amrex::RealBox {
  return ::amrex::RealBox{{x[0] - rx, x[1] - r_tube, x[2] - r_tube},
                          {x[0] + rx, x[1] + r_tube, x[2] + r_tube}};
}

struct TubeSolverOptions {
  int n_cells{200};
  int max_refinement_level{1};
  std::array<double, 2> x_domain{0.53, 2.0};
  double phi{0.0};
};

struct NoInit {
  static void InitializeData(const ::amrex::MultiFab&,
                             const ::amrex::Geometry&) noexcept {}
};

auto MakeTubeSolver(fub::Burke2012& mechanism, const TubeSolverOptions& opts,
                    const std::map<std::string, pybind11::object>& options,
                    int k) {
  const std::array<int, AMREX_SPACEDIM> n_cells{opts.n_cells, 1, 1};
  const double x_lo = opts.x_domain[0];
  const double x_up = opts.x_domain[1];
  const double x_len = x_up - x_lo;
  const double r_len = 0.5 * x_len;
  const double x_mid = 0.5 * x_lo + 0.5 * x_up;
  amrex::RealBox xbox = DomainAroundCenter(Center(x_mid, opts.phi), r_len);
  const std::array<int, AMREX_SPACEDIM> periodicity{0, 0, 0};

  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};

  using namespace fub::amrex;

  CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = xbox;

  DataDescription desc = MakeDataDescription(equation);

  PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = opts.max_refinement_level;
  hier_opts.refine_ratio = amrex::IntVect{2, 1, 1};

  amrex::Geometry geom(
      amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1}}, &xbox,
      -1, periodicity.data());
  geom.refine(hier_opts.refine_ratio);

  using Complete = fub::IdealGasMix<1>::Complete;
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 1e-3),
                            std::make_pair(&Complete::pressure, 1e-2),
                            std::make_pair(&Complete::temperature, 1e-1)};

  ::amrex::Box refine_box{{opts.n_cells - 5, 0, 0}, {opts.n_cells - 1, 0, 0}};
  ConstantBox constant_box{refine_box};

  equation.GetReactor().SetMoleFractions("N2:79,O2:21");
  equation.GetReactor().SetTemperature(300.0);
  equation.GetReactor().SetPressure(101325.0);
  fub::Complete<fub::IdealGasMix<Tube_Rank>> state(equation);
  equation.CompleteFromReactor(state);
  ConstantData initial_data{equation, state};

  PressureValveOptions valve_opts(options, fmt::format("valve{}", k));
  PressureValveBoundary valve{equation, valve_opts};
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
          PatchHierarchy(desc, geometry, hier_opts), initial_data,
          TagAllOf(gradient, constant_box), boundaries);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint = fmt::format("{}/Tube_{}", checkpoint, k);
      PatchHierarchy h =
          ReadCheckpointFile(checkpoint, desc, geometry, hier_opts);
      std::shared_ptr<GriddingAlgorithm> gridding =
          std::make_shared<GriddingAlgorithm>(std::move(h), NoInit{},
                                              TagAllOf(gradient, constant_box),
                                              boundaries);
      return gridding;
    }
  }();

  fub::ideal_gas::MusclHancockPrimMethod<Tube_Rank> flux_method{equation};
  // fub::EinfeldtSignalVelocities<fub::IdealGasMix<Tube_Rank>> signals{};
  // fub::HllMethod hll_method{equation, signals};
  HyperbolicMethod method{FluxMethod(fub::execution::openmp, flux_method),
                          ForwardIntegrator(fub::execution::openmp),
                          Reconstruction(fub::execution::openmp, equation)};

  return std::make_pair(fub::amrex::IntegratorContext(gridding, method), valve);
}

::amrex::Box BoxWhichContains(const ::amrex::RealBox& xbox,
                              const ::amrex::Geometry& geom) {
  ::amrex::Box domain = geom.Domain();
  ::amrex::IntVect lo = domain.smallEnd();
  ::amrex::IntVect up = domain.bigEnd();
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    for (int i = domain.smallEnd(d); i < domain.bigEnd(d); ++i) {
      const double x = geom.CellCenter(i, d);
      if (x < xbox.lo(d)) {
        lo[d] = std::max(lo[d], i);
      }
      if (x > xbox.hi(d)) {
        up[d] = std::min(up[d], i);
      }
    }
  }
  return ::amrex::Box{lo, up};
}

auto MakePlenumSolver(fub::Burke2012& mechanism, int num_cells, int n_level,
                      const std::map<std::string, pybind11::object>& options) {
  const std::array<int, Plenum_Rank> n_cells{num_cells, num_cells, num_cells};
  const std::array<double, Plenum_Rank> xlower{0.03, -0.5 * 0.50, -0.5 * 0.50};
  const std::array<double, Plenum_Rank> xupper{+0.53, +0.5 * 0.50, +0.5 * 0.50};
  const std::array<int, Plenum_Rank> periodicity{0, 0, 0};

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(amrex::Box{{}, ::amrex::IntVect(num_cells - 1)},
                              &xbox, -1, periodicity.data());

  auto embedded_boundary = amrex::EB2::makeUnion(
      amrex::EB2::makeIntersection(
          amrex::EB2::CylinderIF(r_outer, 0.5, 0, {0.25, 0.0, 0.0}, true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(0.6, 0.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(0.6, 1.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(0.6, 2.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(0.6, 3.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(0.6, 4.0 * alpha),
                                 true)),
      amrex::EB2::CylinderIF(r_inner, 1.0, 0, {0.25, 0.0, 0.0}, false));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};

  // Make Gridding Algorithm

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  equation.GetReactor().SetMoleFractions("N2:79,O2:21");
  equation.GetReactor().SetTemperature(300.0);
  equation.GetReactor().SetPressure(101325.0);
  fub::Complete<fub::IdealGasMix<Plenum_Rank>> right(equation);
  equation.CompleteFromReactor(right);

  using namespace fub::amrex::cutcell;

  fub::amrex::cutcell::RiemannProblem initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04), right, right);

  PatchHierarchyOptions hier_opts{};
  hier_opts.max_number_of_levels = n_level;
  hier_opts.index_spaces = MakeIndexSpaces(shop, coarse_geom, n_level);

  using State = fub::Complete<fub::IdealGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.01}};

  ::amrex::RealBox inlet{{-0.1, -0.5, -0.5}, {0.05, +0.5, +0.5}};
  ::amrex::RealBox outlet{{0.5, -0.5, -0.5}, {0.54, +0.5, +0.5}};
  const ::amrex::Box inlet_box = BoxWhichContains(inlet, coarse_geom);
  const ::amrex::Box outlet_box = BoxWhichContains(outlet, coarse_geom);
  ConstantBox constant_box{outlet_box};

  const double required_massflow = 0.5; // [kg / s]
  const double surface_area = M_PI * (r_outer * r_outer - r_inner * r_inner);
  BoundarySet boundary_condition{
      {MassflowBoundary{"Massflow", equation, inlet_box, required_massflow,
                        surface_area, fub::Direction::X, 0},
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
          PatchHierarchy(equation, geometry, hier_opts), initial_data,
          TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
          boundary_condition);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint += "/Plenum";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), geometry,
          hier_opts);
      return std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
          boundary_condition);
    }
  }();

  // Make Solver

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  // fub::ideal_gas::MusclHancockPrimMethod<Plenum_Rank> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(hll_method, hll_method);

  HyperbolicMethod method{
      FluxMethod{fub::execution::openmp_simd, cutcell_method},
      fub::amrex::cutcell::TimeIntegrator{},
      Reconstruction{fub::execution::openmp_simd, equation}};

  return fub::amrex::cutcell::IntegratorContext(gridding, method);
}

struct ProgramOptions {
  ProgramOptions() = default;

  explicit ProgramOptions(const std::map<std::string, pybind11::object>& vm) {
    plenum_n_cells = fub::GetOptionOr(vm, "plenum_n_cells", plenum_n_cells);
    max_refinement_level =
        fub::GetOptionOr(vm, "max_number_of_levels", max_refinement_level);
    checkpoint = fub::GetOptionOr(vm, "checkpoint", checkpoint);
    constexpr double tube_len_over_plenum_len = 1.47 / 0.56;
    tube_n_cells = static_cast<int>(tube_len_over_plenum_len *
                                    static_cast<double>(plenum_n_cells));
    tube_n_cells = tube_n_cells - tube_n_cells % 8;
  }

  template <typename Logger> void Print(Logger& log) const {
    BOOST_LOG(log) << "Grid Options:";
    BOOST_LOG(log) << "  - plenum_n_cells = " << plenum_n_cells;
    BOOST_LOG(log) << "  - tube_n_cells = " << tube_n_cells;
    BOOST_LOG(log) << "  - max_refinement_level = " << max_refinement_level;

    if (!checkpoint.empty()) {
      BOOST_LOG(log) << "Restart simulation from checkpoint '" << checkpoint
                     << "'!";
    }
  }

  int plenum_n_cells{128};
  int tube_n_cells{};
  int max_refinement_level{1};
  std::string checkpoint{};
};

std::optional<std::map<std::string, pybind11::object>>
ParseCommandLine(int argc, char** argv) {
  namespace po = boost::program_options;
  po::options_description desc{};
  std::string config_path{};
  desc.add_options()("config", po::value<std::string>(&config_path),
                     "Path to the config file which can be parsed.");
  po::variables_map vm;
  std::map<std::string, pybind11::object> options{};
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("config")) {
      config_path = vm["config"].as<std::string>();
      options = fub::ParsePythonScript(config_path, MPI_COMM_WORLD);
    }
    po::notify(vm);
  } catch (std::exception& e) {
    amrex::Print()
        << "[Error] An Error occured while reading program options:\n";
    amrex::Print() << e.what();
    return {};
  }

  if (vm.count("help")) {
    amrex::Print() << desc << "\n";
    return {};
  }

  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);

  fub::RunOptions(options).Print(log);
  ProgramOptions(options).Print(log);
  for (int i = 0; i < 5; ++i) {
    fub::amrex::PressureValveOptions(options, fmt::format("valve{}", i))
        .Print(log);
  }
  fub::amrex::IgniteDetonationOptions(options, "ignite").Print(log);
  return options;
}

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
    oa << ignition.GetLastIgnitionTimePoints();
  }
}

void MyMain(const std::map<std::string, pybind11::object>& vm) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::Burke2012 mechanism{};

  ProgramOptions po(vm);

  const int n_level = po.max_refinement_level;
  auto plenum = MakePlenumSolver(mechanism, po.plenum_n_cells, n_level, vm);

  int tube_n_cells = po.tube_n_cells;
  std::vector<fub::amrex::IntegratorContext> tubes;

  std::vector<fub::amrex::BlockConnection> connectivity{};
  std::vector<std::shared_ptr<fub::amrex::PressureValve>> valves{};

  auto MakeConnection = [&](int k) {
    TubeSolverOptions opts{};
    opts.max_refinement_level = po.max_refinement_level;
    opts.n_cells = tube_n_cells;
    opts.phi = k * alpha;
    auto&& [tube, valve] = MakeTubeSolver(mechanism, opts, vm, k);
    tubes.push_back(std::move(tube));
    valves.push_back(valve.GetSharedState());
    fub::amrex::BlockConnection connection;
    connection.direction = fub::Direction::X;
    connection.side = 1;
    connection.plenum.id = 0;
    connection.tube.id = k;
    connection.tube.mirror_box = tubes[k]
                                     .GetGriddingAlgorithm()
                                     ->GetPatchHierarchy()
                                     .GetGeometry(0)
                                     .Domain();
    connection.plenum.mirror_box =
        BoxWhichContains(DomainAroundCenter(Center(+0.53, k * alpha), 0.03),
                         plenum.GetGeometry(0));
    return connection;
  };

  connectivity.push_back(MakeConnection(0));
  connectivity.push_back(MakeConnection(1));
  connectivity.push_back(MakeConnection(2));
  connectivity.push_back(MakeConnection(3));
  connectivity.push_back(MakeConnection(4));

  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};
  fub::IdealGasMix<Plenum_Rank> equation{mechanism};

  fub::amrex::MultiBlockIntegratorContext context(
      fub::FlameMasterReactor(mechanism), std::move(tubes), {std::move(plenum)},
      std::move(connectivity), valves);

  fub::DimensionalSplitLevelIntegrator system_solver(fub::int_c<Plenum_Rank>,
                                                     context);

  fub::amrex::MultiBlockIgniteDetonation ignition{
      tube_equation, context.GetGriddingAlgorithm(),
      fub::amrex::IgniteDetonationOptions(vm, "ignite")};

  std::string checkpoint{};
  if (vm.count("checkpoint")) {
    checkpoint = vm.at("checkpoint").cast<std::string>();
  }
  if (!checkpoint.empty()) {
    MPI_Comm comm = context.GetMpiCommunicator();
    std::string input =
        fub::ReadAndBroadcastFile(checkpoint + "/Ignition", comm);
    std::istringstream ifs(input);
    boost::archive::text_iarchive ia(ifs);
    std::vector<fub::Duration> last_ignitions;
    ia >> last_ignitions;
    ignition.SetLastIgnitionTimePoints(last_ignitions);
  }

  fub::SplitSystemSourceLevelIntegrator ign_solver(system_solver, ignition,
                                                   fub::GodunovSplitting{});

  fub::amrex::MultiBlockKineticSouceTerm source_term{
      fub::IdealGasMix<Tube_Rank>{mechanism}, context.GetGriddingAlgorithm()};

  fub::SplitSystemSourceLevelIntegrator level_integrator{ign_solver,
                                                         source_term};

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  fub::OutputFactory<fub::amrex::MultiBlockGriddingAlgorithm> factory{};
  factory.RegisterOutput<fub::amrex::MultiWriteHdf5>("HDF5");
  fub::MultipleOutputs<fub::amrex::MultiBlockGriddingAlgorithm> outputs(
      std::move(factory),
      fub::ToMap(fub::GetOptionOr(vm, "output", pybind11::dict{})));

  outputs(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, fub::RunOptions(vm), wall_time_reference, outputs);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
  pybind11::scoped_interpreter interpreter{};
  {
    fub::amrex::ScopeGuard _{};
    auto vm = ParseCommandLine(argc, argv);
    if (vm) {
      MyMain(*vm);
    }
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}
