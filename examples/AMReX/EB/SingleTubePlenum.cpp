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
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/algorithm/string/replace.hpp>
#include <boost/core/null_deleter.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/manipulators/add_value.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/file.hpp>

#include <pybind11/stl.h>

#include <cmath>
#include <iostream>

#include <xmmintrin.h>

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = AMREX_SPACEDIM;

static constexpr double r_tube = 0.015;

template <typename T>
using ProbesView =
    fub::basic_mdspan<T, fub::extents<AMREX_SPACEDIM, fub::dynamic_extent>>;

struct ProgramOptions {
  ProgramOptions() = default;

  explicit ProgramOptions(const std::map<std::string, pybind11::object>& vm) {
    std::map<std::string, pybind11::object> grid =
        fub::ToMap(fub::GetOptionOr(vm, "grid", pybind11::dict()));
    std::vector xs =
        fub::GetOptionOr(grid, "x_range", std::vector{-0.03, 0.57});
    std::vector ys =
        fub::GetOptionOr(grid, "y_range", std::vector{-0.15, 0.15});
    std::vector zs =
        fub::GetOptionOr(grid, "z_range", std::vector{-0.15, 0.15});
    auto check_size = [](auto& xs, int n, const char* name) {
      if (int(xs.size()) != n) {
        throw std::runtime_error(
            fmt::format("Option '{}' need exactly {} values.", name, n));
      }
    };
    check_size(xs, 2, "grid.x_range");
    check_size(ys, 2, "grid.y_range");
    check_size(zs, 2, "grid.z_range");
    plenum_xbox = amrex::RealBox(xs[0], ys[0], zs[0], xs[1], ys[1], zs[1]);
    std::vector n_cells =
        fub::GetOptionOr(grid, "n_cells", std::vector{64, 64, 64});
    check_size(n_cells, 3, "grid.n_cells");
    plenum_n_cells = amrex::IntVect{n_cells[0], n_cells[1], n_cells[2]};
    n_levels = fub::GetOptionOr(grid, "max_number_of_levels", n_levels);
    checkpoint = fub::GetOptionOr(vm, "checkpoint", checkpoint);

    tube_xbox = amrex::RealBox(-1.5, -r_tube, -r_tube, xs[0], +r_tube, r_tube);
    const double tube_length = tube_xbox.hi(0) - tube_xbox.lo(0);
    const double plenum_domain_length = plenum_xbox.hi(0) - plenum_xbox.lo(0);
    const double t_over_p = tube_length / plenum_domain_length;
    tube_n_cells = int(double(plenum_n_cells[0]) * t_over_p);
    tube_n_cells = tube_n_cells - tube_n_cells % 8;

    auto plenum = fub::ToMap(fub::GetOptionOr(vm, "plenum", pybind11::dict()));
    plenum_temperature =
        fub::GetOptionOr(plenum, "temperature", plenum_temperature);
    plenum_radius = fub::GetOptionOr(plenum, "radius", plenum_radius);
    plenum_outlet_radius =
        fub::GetOptionOr(plenum, "outlet_radius", plenum_outlet_radius);
    plenum_outlet_exit_radius =
        fub::GetOptionOr(plenum, "outlet_exit_radius", plenum_outlet_exit_radius);
    plenum_outlet_exit_length =
        fub::GetOptionOr(plenum, "outlet_exit_length", plenum_outlet_exit_length);
    plenum_outlet_length =
        fub::GetOptionOr(plenum, "outlet_length", plenum_outlet_length);
    plenum_length = fub::GetOptionOr(plenum, "length", plenum_length);

    auto output = fub::ToMap(fub::GetOptionOr(vm, "output", pybind11::dict()));
    output_directory = fub::GetOptionOr(output, "directory", output_directory);
  }

  template <typename Logger> void Print(Logger& log) const {
    BOOST_LOG(log) << "Grid Options:";
    BOOST_LOG(log) << fmt::format(
        "  - plenum_xbox = {{{{{}, {}, {}}}, {{{}, {}, {}}}}} [m]",
        plenum_xbox.lo(0), plenum_xbox.lo(1), plenum_xbox.lo(2),
        plenum_xbox.hi(0), plenum_xbox.hi(1), plenum_xbox.hi(2));
    std::array<int, 3> n_cells{plenum_n_cells[0], plenum_n_cells[1],
                               plenum_n_cells[2]};
    BOOST_LOG(log) << fmt::format("  - plenum_n_cells = {{{}}} [-]",
                                  fmt::join(n_cells, ", "));
    BOOST_LOG(log) << fmt::format(
        "  - tube_xbox = {{{{{}, {}, {}}}, {{{}, {}, {}}}}} [m]",
        tube_xbox.lo(0), tube_xbox.lo(1), tube_xbox.lo(2), tube_xbox.hi(0),
        tube_xbox.hi(1), tube_xbox.hi(2));
    BOOST_LOG(log) << "  - tube_n_cells = " << tube_n_cells << " [-]";
    BOOST_LOG(log) << "  - n_levels = " << n_levels << " [-]";

    BOOST_LOG(log) << "Problem Options:";
    BOOST_LOG(log) << "  - plenum.radius = " << plenum_radius << " [m]";
    BOOST_LOG(log) << "  - plenum.temperature = " << plenum_temperature
                   << " [K]";
    BOOST_LOG(log) << "  - plenum.outlet_radius= " << plenum_outlet_radius
                   << " [m]";
    BOOST_LOG(log) << "  - plenum.outlet_length = " << plenum_outlet_length << " [m]";
    BOOST_LOG(log) << "  - plenum.outlet_exit_radius = " << plenum_outlet_exit_radius << " [m]";
    BOOST_LOG(log) << "  - plenum.outlet_exit_length = " << plenum_outlet_exit_length << " [m]";
    BOOST_LOG(log) << "  - plenum.length = " << plenum_length << " [m]";

    if (!checkpoint.empty()) {
      BOOST_LOG(log) << "Restart simulation from checkpoint '" << checkpoint
                     << "'!";
    }
  }

  ::amrex::RealBox plenum_xbox{};
  ::amrex::IntVect plenum_n_cells{64, 64, 64};
  ::amrex::RealBox tube_xbox{};
  int tube_n_cells{};
  int n_levels{1};
  std::string checkpoint{};
  double plenum_radius{0.5};
  double plenum_temperature{300};
  double plenum_outlet_radius{r_tube};
  double plenum_outlet_exit_radius{plenum_outlet_radius};
  double plenum_outlet_exit_length{plenum_outlet_radius};
  double plenum_length{0.25};
  double plenum_outlet_length{0.04+0.04+0.03};
  std::string output_directory{"SingleTubePlenum"};
};

struct TubeSolverOptions {
  int n_cells{200};
  int max_refinement_level{1};
  std::array<double, 2> x_domain{-1.5, -0.03};
  double phi{0.0};
};

auto DomainAroundCenter(const ::amrex::RealArray& x, double rx)
    -> ::amrex::RealBox {
  return ::amrex::RealBox{{x[0] - rx, x[1] - r_tube, x[2] - r_tube},
                          {x[0] + rx, x[1] + r_tube, x[2] + r_tube}};
}

auto MakeTubeSolver(const ProgramOptions& po, fub::Burke2012& mechanism,
                    const std::map<std::string, pybind11::object>& vm) {
  const std::array<int, AMREX_SPACEDIM> n_cells{po.tube_n_cells, 1, 1};
  amrex::RealBox xbox = po.tube_xbox;
  const std::array<int, AMREX_SPACEDIM> periodicity{0, 0, 0};

  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};

  using namespace fub::amrex;

  CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = xbox;

  DataDescription desc = MakeDataDescription(equation);

  PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = po.n_levels;
  hier_opts.refine_ratio = amrex::IntVect{2, 1, 1};

  amrex::Geometry geom(amrex::Box{{}, {n_cells[0] - 1, 0, 0}}, &xbox, -1,
                       periodicity.data());
  geom.refine(hier_opts.refine_ratio);

  using Complete = fub::IdealGasMix<1>::Complete;
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 1e-3),
                            std::make_pair(&Complete::pressure, 1e-2),
                            std::make_pair(&Complete::temperature, 1e-1)};

  ::amrex::Box refine_box{{n_cells[0] - 5, 0, 0}, {n_cells[0] - 1, 0, 0}};
  ConstantBox constant_box{refine_box};

  equation.GetReactor().SetMoleFractions("N2:79,O2:21");
  equation.GetReactor().SetTemperature(po.plenum_temperature);
  equation.GetReactor().SetPressure(101325.0);
  fub::Complete<fub::IdealGasMix<Tube_Rank>> state(equation);
  equation.CompleteFromReactor(state);
  ConstantData initial_data{equation, state};

  PressureValveOptions valve_opts(vm, "valve");
  PressureValveBoundary valve{equation, valve_opts};
  BoundarySet boundaries{{valve}};

  // If a checkpoint path is specified we will fill the patch hierarchy with
  // data from the checkpoint file, otherwise we will initialize the data by
  // the initial data function.
  std::shared_ptr<GriddingAlgorithm> gridding = [&] {
    std::string checkpoint = fub::GetOptionOr(vm, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(desc, geometry, hier_opts), initial_data,
          TagAllOf(gradient, constant_box), boundaries);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint = fmt::format("{}/Tube", checkpoint);
      PatchHierarchy h =
          ReadCheckpointFile(checkpoint, desc, geometry, hier_opts);
      std::shared_ptr<GriddingAlgorithm> gridding =
          std::make_shared<GriddingAlgorithm>(std::move(h), initial_data,
                                              TagAllOf(gradient, constant_box),
                                              boundaries);
      return gridding;
    }
  }();

  fub::ideal_gas::MusclHancockPrimMethod<Tube_Rank> flux_method{equation};
  HyperbolicMethod method{FluxMethod(fub::execution::openmp, flux_method),
                          ForwardIntegrator(fub::execution::openmp),
                          Reconstruction(fub::execution::openmp, equation)};

  return std::pair{fub::amrex::IntegratorContext(gridding, method), valve};
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

auto MakePlenumSolver(const ProgramOptions& po, fub::Burke2012& mechanism,
                      const std::map<std::string, pybind11::object>& vm) {
  const std::array<int, Plenum_Rank> n_cells{
      po.plenum_n_cells[0], po.plenum_n_cells[1], po.plenum_n_cells[2]};
  const int n_level = po.n_levels;

  const std::array<int, Plenum_Rank> periodicity{0, 0, 0};

  amrex::RealBox xbox = po.plenum_xbox;
  amrex::Geometry coarse_geom(
      amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1}}, &xbox,
      -1, periodicity.data());

  auto embedded_boundary = amrex::EB2::makeIntersection(
      amrex::EB2::CylinderIF(po.plenum_radius, po.plenum_length, 0,
                             {0.5 * po.plenum_length, 0.0, 0.0}, true),
      amrex::EB2::CylinderIF(r_tube, 0.2, 0, {1e-6, 0.0, 0.0}, true),
      amrex::EB2::CylinderIF(po.plenum_outlet_radius, 0.2, 0,
                             {po.plenum_length, 0.0, 0.0}, true),
      fub::amrex::Geometry(fub::ConeIF({po.plenum_length, 0.0, 0.0},
                                       po.plenum_radius, 0.04, true)),
      amrex::EB2::translate(amrex::EB2::rotate(fub::amrex::Geometry(fub::ConeIF({0.0, 0.0, 0.0},
                                       po.plenum_outlet_exit_radius, po.plenum_outlet_exit_length, true)), M_PI, 2), {po.plenum_length+po.plenum_outlet_length, 0.0, 0.0}));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};

  // Make Gridding Algorithm

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = po.plenum_xbox;
  geometry.periodicity = periodicity;

  equation.GetReactor().SetMoleFractions("N2:79,O2:21");
  equation.GetReactor().SetTemperature(po.plenum_temperature);
  equation.GetReactor().SetPressure(101325.0);
  fub::Complete<fub::IdealGasMix<Plenum_Rank>> right(equation);
  equation.CompleteFromReactor(right);

  using namespace fub::amrex::cutcell;

  fub::amrex::cutcell::RiemannProblem initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.0), right, right);

  PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces = MakeIndexSpaces(shop, coarse_geom, n_level);

  using State = fub::Complete<fub::IdealGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.01},
                             std::pair{&State::density, 0.05}};

  ::amrex::RealBox inlet{{xbox.lo(0), -r_tube, -r_tube},
                         {0.01, +r_tube, +r_tube}};
  ::amrex::Box refine_box = BoxWhichContains(inlet, coarse_geom);
  ConstantBox constant_in_box{refine_box};

  ::amrex::RealBox outlet{{0.5, xbox.lo(1), xbox.lo(2)},
                          {xbox.hi(0), xbox.hi(1), xbox.hi(2)}};
  refine_box = BoxWhichContains(outlet, coarse_geom);
  ConstantBox constant_out_box{refine_box};

  outlet = amrex::RealBox{{xbox.hi(0) - 0.01, xbox.lo(1), xbox.lo(2)},
                          {xbox.hi(0), xbox.hi(1), xbox.hi(2)}};
  const ::amrex::Box outlet_box = BoxWhichContains(outlet, coarse_geom);

  //  const double p0 = 101325.0;
  BoundarySet boundary_condition{
      {TransmissiveBoundary{fub::Direction::X, 0},
       IsentropicPressureBoundary{"RightPlenumBoundary", equation, outlet_box,
                                  101325.0, fub::Direction::X, 1}}};

  std::shared_ptr gridding = [&] {
    std::string checkpoint = fub::GetOptionOr(vm, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(equation, geometry, options), initial_data,
          TagAllOf(TagCutCells(), gradients, constant_in_box, constant_out_box,
                   TagBuffer(2)),
          boundary_condition);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint += "/Plenum";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), geometry,
          options);
      return std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), gradients, constant_in_box, constant_out_box,
                   TagBuffer(2)),
          boundary_condition);
    }
  }();

  // Make Solver

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  //  fub::ideal_gas::MusclHancockPrimMethod<Plenum_Rank> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(hll_method, hll_method);

  HyperbolicMethod method{
      FluxMethod{fub::execution::openmp_simd, cutcell_method},
      fub::amrex::cutcell::TimeIntegrator{},
      Reconstruction{fub::execution::openmp_simd, equation}};

  return fub::amrex::cutcell::IntegratorContext(gridding, method);
}

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
  fub::amrex::PressureValveOptions(options).Print(log);
  fub::amrex::IgniteDetonationOptions(options, "ignite").Print(log);
  return options;
}

void MyMain(const std::map<std::string, pybind11::object>& vm);

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

void WriteCheckpoint(const std::string& path,
                     const fub::amrex::MultiBlockGriddingAlgorithm& grid,
                     std::shared_ptr<fub::amrex::PressureValve> valve, int rank,
                     const fub::amrex::MultiBlockIgniteDetonation& ignition) {
  auto tubes = grid.GetTubes();
  std::string name = fmt::format("{}/Tube", path);
  fub::amrex::WriteCheckpointFile(name, tubes[0]->GetPatchHierarchy());
  if (rank == 0) {
    name = fmt::format("{}/Valve", path);
    std::ofstream valve_checkpoint(name);
    boost::archive::text_oarchive oa(valve_checkpoint);
    oa << *valve;
  }
  name = fmt::format("{}/Plenum", path);
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

  const ProgramOptions po(vm);

  auto plenum = MakePlenumSolver(po, mechanism, vm);
  auto [tube, valve] = MakeTubeSolver(po, mechanism, vm);
  auto valve_state = valve.GetSharedState();

  fub::amrex::BlockConnection connection;
  connection.direction = fub::Direction::X;
  connection.side = 0;
  connection.plenum.id = 0;
  connection.plenum.mirror_box = plenum.GetGriddingAlgorithm()
                                     ->GetPatchHierarchy()
                                     .GetGeometry(0)
                                     .Domain();
  connection.tube.id = 0;
  connection.tube.mirror_box =
      tube.GetGriddingAlgorithm()->GetPatchHierarchy().GetGeometry(0).Domain();

  fub::IdealGasMix<Plenum_Rank> plenum_equation{mechanism};
  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};

  fub::amrex::MultiBlockIntegratorContext context(
      fub::FlameMasterReactor(mechanism), {std::move(tube)},
      {std::move(plenum)}, {connection});

  fub::DimensionalSplitLevelIntegrator system_solver(fub::int_c<Plenum_Rank>,
                                                     std::move(context));

  fub::amrex::MultiBlockIgniteDetonation ignition{
      tube_equation, context.GetGriddingAlgorithm(),
      fub::amrex::IgniteDetonationOptions(vm, "ignite")};

  std::string checkpoint = fub::GetOptionOr(vm, "checkpoint", std::string{});
  if (!checkpoint.empty()) {
    MPI_Comm comm = context.GetMpiCommunicator();
    std::string input =
        fub::ReadAndBroadcastFile(checkpoint + "/Ignition", comm);
    std::istringstream ifs(input);
    {
      boost::archive::text_iarchive ia(ifs);
      std::vector<fub::Duration> last_ignitions;
      ia >> last_ignitions;
      ignition.SetLastIgnitionTimePoints(last_ignitions);
    }
    input = fub::ReadAndBroadcastFile(checkpoint + "/Valve", comm);
    ifs = std::istringstream(input);
    boost::archive::text_iarchive ia(ifs);
    ia >> *valve.GetSharedState();
  }

  fub::DimensionalSplitSystemSourceSolver ign_solver(system_solver, ignition,
                                                     fub::GodunovSplitting{});

  fub::amrex::MultiBlockKineticSouceTerm source_term{
      fub::IdealGasMix<Tube_Rank>{mechanism}, context.GetGriddingAlgorithm()};

  fub::DimensionalSplitSystemSourceSolver solver{ign_solver, source_term};

  std::string base_name = po.output_directory;

  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  using namespace fub::amrex;
  fub::OutputFactory<MultiBlockGriddingAlgorithm> factory{};
  factory.RegisterFactory<MultiWriteHdf5>("HDF5");
  factory.RegisterFactory<MultiBlockPlotfileOutput>("Plotfiles");
  factory.RegisterFactory<LogProbesOutput>("Probes");
  factory.RegisterFactory<fub::InvokeFunction<MultiBlockGriddingAlgorithm>>(
      "Checkpoint", [&](const MultiBlockGriddingAlgorithm& grid) {
        std::string name =
            fmt::format("{}/Checkpoint/{:05}", base_name, grid.GetCycles());
        amrex::Print() << "Write Checkpoint to '" << name << "'!\n";
        WriteCheckpoint(name, grid, valve_state, rank, ignition);
      });
  fub::MultipleOutputs<MultiBlockGriddingAlgorithm> outputs(
      std::move(factory),
      fub::ToMap(fub::GetOptionOr(vm, "output", pybind11::dict{})));

  outputs(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, fub::RunOptions(vm), wall_time_reference, outputs);
}
