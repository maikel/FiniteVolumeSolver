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

  explicit ProgramOptions(const boost::program_options::variables_map& vm) {
    auto GetOptionOr = [&](const char* opt, auto default_value) {
      if (vm.count(opt)) {
        return vm[opt].as<std::decay_t<decltype(default_value)>>();
      }
      return default_value;
    };
    std::vector xs = GetOptionOr("grid.x_range", std::vector{-0.03, 0.57});
    std::vector ys = GetOptionOr("grid.y_range", std::vector{-0.15, 0.15});
    std::vector zs = GetOptionOr("grid.z_range", std::vector{-0.15, 0.15});
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
    std::vector n_cells = GetOptionOr("grid.n_cells", std::vector{64, 64, 64});
    check_size(n_cells, 3, "grid.n_cells");
    plenum_n_cells = amrex::IntVect{n_cells[0], n_cells[1], n_cells[2]};
    n_levels = GetOptionOr("grid.max_number_of_levels", n_levels);
    checkpoint = GetOptionOr("grid.checkpoint", checkpoint);

    tube_xbox = amrex::RealBox(-1.5, -r_tube, -r_tube, xs[0], +r_tube, +r_tube);
    const double tube_length = tube_xbox.hi(0) - tube_xbox.lo(0);
    const double plenum_domain_length = plenum_xbox.hi(0) - plenum_xbox.lo(0);
    const double t_over_p = tube_length / plenum_domain_length;
    tube_n_cells = int(double(plenum_n_cells[0]) * t_over_p);
    tube_n_cells = tube_n_cells - tube_n_cells % 8;

    plenum_temperature = GetOptionOr("plenum.temperature", plenum_temperature);
    plenum_outlet_radius = GetOptionOr("plenum.outlet_radius", plenum_outlet_radius);
    plenum_length = GetOptionOr("plenum.length", plenum_length);
    plenum_jump = GetOptionOr("plenum.jump", plenum_jump);

    output_directory = GetOptionOr("output.directory", output_directory);
    checkpoint_interval = fub::Duration(GetOptionOr("output.checkpoint_interval", checkpoint_interval.count()));
    x_probes = GetOptionOr("output.x_probes", x_probes);
  }

  static boost::program_options::options_description GetCommandLineOptions() {
    namespace po = boost::program_options;
    po::options_description desc{"Grid Options"};
    // clang-format off
    desc.add_options()
    ("grid.n_cells", po::value<std::vector<int>>()->multitoken(), "Base number of cells in the plenum for the coarsest level")
    ("grid.x_range", po::value<std::vector<double>>()->multitoken(), "Coordinate range for the x coordinate direction.")
    ("grid.y_range", po::value<std::vector<double>>()->multitoken(), "Coordinate range for the y coordinate direction.")
    ("grid.z_range", po::value<std::vector<double>>()->multitoken(), "Coordinate range for the z coordinate direction.")
    ("grid.max_number_of_levels", po::value<int>(), "Maximal number of refinement levels across tube and plenum domains.")
    ("grid.checkpoint", po::value<std::string>(), "The path to the checkpoint files. Used to restart a simulation.");

    po::options_description prob_desc{"Problem Options"};
    prob_desc.add_options()
    ("plenum.temperature", po::value<double>(), "Temperature value for the plenum")
    ("plenum.outlet_radius", po::value<double>(), "Radius of plenum outlet")
    ("plenum.length", po::value<double>(), "Length of plenum body")
    ("plenum.jump", po::value<double>(), "Diameter Jump from plenum inlet to plenum body")
    ("output.directory", po::value<double>(), "Output directory")
    ("output.checkpoint_interval", po::value<double>(), "Time interval to make checkpoint ")
    ("output.x_probes", po::value<std::vector<double>>()->multitoken(), "Coordinates in x direction for locations where the state is measured in each time step.");
    // clang-format on
    desc.add(prob_desc);
    return desc;
  }

  template <typename Logger> void Print(Logger& log) const {
    BOOST_LOG(log) << "Grid Options:";
    BOOST_LOG(log) << fmt::format(
        "  - plenum_xbox = {{{{{}, {}, {}}}, {{{}, {}, {}}}}}",
        plenum_xbox.lo(0), plenum_xbox.lo(1), plenum_xbox.lo(2),
        plenum_xbox.hi(0), plenum_xbox.hi(1), plenum_xbox.hi(2));
    std::array<int, 3> n_cells{plenum_n_cells[0], plenum_n_cells[1],
                               plenum_n_cells[2]};
    BOOST_LOG(log) << fmt::format("  - plenum_n_cells = {{{}}}",
                                  fmt::join(n_cells, ", "));
    BOOST_LOG(log) << fmt::format(
        "  - tube_xbox = {{{{{}, {}, {}}}, {{{}, {}, {}}}}}", tube_xbox.lo(0),
        tube_xbox.lo(1), tube_xbox.lo(2), tube_xbox.hi(0), tube_xbox.hi(1),
        tube_xbox.hi(2));
    BOOST_LOG(log) << "  - tube_n_cells = " << tube_n_cells;
    BOOST_LOG(log) << "  - n_levels = " << n_levels;

    BOOST_LOG(log) << "Problem Options:";
    BOOST_LOG(log) << "  - plenum.temperature = " << plenum_temperature;
    BOOST_LOG(log) << "  - plenum.outlet_radius= " << plenum_outlet_radius;
    BOOST_LOG(log) << "  - plenum.length= " << plenum_length;
    BOOST_LOG(log) << "  - plenum.jump= " << plenum_jump;
    BOOST_LOG(log) << "  - output.directory = '" << output_directory << "'";
    BOOST_LOG(log) << "  - output.checkpoint_interval = " << checkpoint_interval.count() << " [s]";
    BOOST_LOG(log) << fmt::format("  - output.x_probes = {{{}}}",
                                  fmt::join(x_probes, ", "));

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
  double plenum_temperature{300};
  double plenum_outlet_radius{r_tube};
  double plenum_jump{2 * r_tube};
  double plenum_length{0.25};
  fub::Duration checkpoint_interval{0.0001};
  std::string output_directory{"SingleTube"};
  std::vector<double> x_probes{};
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
                    const boost::program_options::variables_map& vm) {
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
    std::string checkpoint{};
    if (vm.count("grid.checkpoint")) {
      checkpoint = vm["grid.checkpoint"].as<std::string>();
    }
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
                      const boost::program_options::variables_map& vm) {
  const std::array<int, Plenum_Rank> n_cells{
      po.plenum_n_cells[0], po.plenum_n_cells[1], po.plenum_n_cells[2]};
  const int n_level = po.n_levels;

  const std::array<int, Plenum_Rank> periodicity{0, 0, 0};

  amrex::RealBox xbox = po.plenum_xbox;
  amrex::Geometry coarse_geom(
      amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1}}, &xbox,
      -1, periodicity.data());

  
auto MakePlenum2D = [](auto... points) {
  fub::Polygon::Vector xs{};
  fub::Polygon::Vector ys{};

  (xs.push_back(std::get<0>(points)), ...);
  (ys.push_back(std::get<1>(points)), ...);

  return fub::Polygon(std::move(xs), std::move(ys));
};

const double r_inner = r_tube;
const double d_tube = 2 * r_tube;

auto plenum2D = MakePlenum2D(
        std::pair{+0.00, r_inner},
        std::pair{+po.plenum_length, r_inner},
        std::pair{+po.plenum_length + 0.03, r_inner + po.plenum_jump + r_tube - po.plenum_outlet_radius},
        std::pair{+1.00, r_inner + po.plenum_jump + r_tube - po.plenum_outlet_radius},
        std::pair{+1.00, r_inner + po.plenum_jump + r_tube + po.plenum_outlet_radius},
        std::pair{+po.plenum_length + 0.03, r_inner + po.plenum_jump + r_tube + po.plenum_outlet_radius},
        std::pair{+po.plenum_length, r_inner + 2*po.plenum_jump + d_tube},
        std::pair{+0.00, r_inner + 2*po.plenum_jump + d_tube},
        std::pair{+0.00, r_inner});

const double r_tube_center = 0.5 * (r_inner + r_inner + 2*po.plenum_jump + d_tube);

auto Center = [r_tube_center](double x, double phi) -> ::amrex::RealArray {
  using std::cos;
  using std::sin;
  return {x, r_tube_center * sin(phi), r_tube_center * cos(phi)};
};

  auto embedded_boundary = amrex::EB2::makeUnion(
      amrex::EB2::makeIntersection(
          fub::amrex::Geometry(fub::Invert(fub::RotateAxis(plenum2D))),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 0.0), true)));
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
  const ::amrex::Box refine_box = BoxWhichContains(inlet, coarse_geom);
  ConstantBox constant_box{refine_box};

  //  const double p0 = 101325.0;
  ::amrex::RealBox outlet{{xbox.hi(0) - 0.01, xbox.lo(1), xbox.lo(2)}, {xbox.hi(0), xbox.hi(1), xbox.hi(2)}};
  const ::amrex::Box outlet_box = BoxWhichContains(outlet, coarse_geom);
  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  IsentropicPressureBoundary{"RightPlenumBoundary", equation, outlet_box,
                                  101325.0, fub::Direction::X, 1}}};

  std::shared_ptr gridding = [&] {
    std::string checkpoint{};
    if (vm.count("grid.checkpoint")) {
      checkpoint = vm["grid.checkpoint"].as<std::string>();
    }
    if (checkpoint.empty()) {
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(equation, geometry, options), initial_data,
          TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
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
          TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
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

void MyMain(const boost::program_options::variables_map& vm);

std::optional<boost::program_options::variables_map>
ParseCommandLine(int argc, char** argv) {
  namespace po = boost::program_options;
  po::options_description desc = fub::RunOptions::GetCommandLineOptions();
  std::string config_path{};
  desc.add_options()("config", po::value<std::string>(&config_path),
                     "Path to the config file which can be parsed.");
  desc.add(ProgramOptions::GetCommandLineOptions());
  desc.add(fub::amrex::PressureValveOptions::GetCommandLineOptions());
  desc.add(
      fub::amrex::IgniteDetonationOptions::GetCommandLineOptions("ignite"));
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("config")) {
      config_path = vm["config"].as<std::string>();
      po::store(po::parse_config_file(config_path.c_str(), desc), vm);
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

  fub::RunOptions(vm).Print(log);
  ProgramOptions(vm).Print(log);
  fub::amrex::PressureValveOptions(vm).Print(log);
  fub::amrex::IgniteDetonationOptions(vm, "ignite").Print(log);
  return vm;
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
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

void WriteCheckpoint(
    const std::string& path,
    const fub::amrex::MultiBlockGriddingAlgorithm& grid,
    std::shared_ptr<fub::amrex::PressureValve> valve,
    int rank, const fub::amrex::MultiBlockIgniteDetonation& ignition) {
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


template <typename Logger>
void LogTubeProbes(Logger& log, ProbesView<const double> probes,
                   const fub::amrex::PatchHierarchy& hierarchy, MPI_Comm comm) {
  std::vector<double> buffer = GatherStates(hierarchy, probes, comm);
  fub::mdspan<const double, 2> states(buffer.data(), probes.extent(1),
                                      buffer.size() / probes.extent(1));
  fub::Burke2012 burke2012{};
  std::array<double, 11> molar_mass{};
  burke2012.getMolarMass(molar_mass);
  auto sH2 = fub::Burke2012::sH2;
  auto sO2 = fub::Burke2012::sO2;
  auto sH2O = fub::Burke2012::sH2O;
  for (int i = 0; i < probes.extent(1); ++i) {
    const double rho = states(i, 0);
    const double u = states(i, 1) / rho;
    const double h2 = states(i, 3 + sH2) / molar_mass[size_t(sH2)];
    const double o2 = states(i, 3 + sO2) / molar_mass[size_t(sO2)];
    const double h2o = states(i, 3 + sH2O) / molar_mass[size_t(sH2O)];
    const double T = states(i, 16);
    const double p = states(i, 14);
    const double a = states(i, 15);
    const double x = probes(0, i);
    using boost::log::add_value;
    using boost::log::trivial::trace;
    BOOST_LOG_SEV(log, trace)
        << add_value("X", x) << add_value("Density", rho)
        << add_value("VelocityX", u) << add_value("Temperature", T)
        << add_value("Pressure", p) << add_value("SpeedOfSound", a)
        << add_value("H2", h2) << add_value("O2", o2) << add_value("H2O", h2o);
  }
}

template <typename Logger>
void LogPlenumProbes(Logger& log, ProbesView<const double> probes,
                     const fub::amrex::cutcell::PatchHierarchy& hierarchy,
                     MPI_Comm comm) {
  std::vector<double> buffer = GatherStates(hierarchy, probes, comm);
  fub::mdspan<const double, 2> states(buffer.data(), probes.extent(1),
                                      buffer.size() / probes.extent(1));
  fub::Burke2012 burke2012{};
  std::array<double, 11> molar_mass{};
  burke2012.getMolarMass(molar_mass);
  auto sH2 = fub::Burke2012::sH2;
  auto sO2 = fub::Burke2012::sO2;
  auto sH2O = fub::Burke2012::sH2O;
  for (int i = 0; i < probes.extent(1); ++i) {
    const double rho = states(i, 0);
    const double u = states(i, 1) / rho;
    const double v = states(i, 2) / rho;
    const double w = states(i, 3) / rho;
    const double h2 = states(i, 5 + sH2) / molar_mass[size_t(sH2)];
    const double o2 = states(i, 5 + sO2) / molar_mass[size_t(sO2)];
    const double h2o = states(i, 5 + sH2O) / molar_mass[size_t(sH2O)];
    const double a = states(i, 17);
    const double T = states(i, 18);
    const double p = states(i, 16);
    const double x = probes(0, i);
    using boost::log::add_value;
    using boost::log::trivial::trace;
    BOOST_LOG_SEV(log, trace)
        << add_value("X", x) << add_value("Density", rho)
        << add_value("VelocityX", u) << add_value("VelocityY", v)
        << add_value("VelocityZ", w) << add_value("Temperature", T)
        << add_value("Pressure", p) << add_value("SpeedOfSound", a)
        << add_value("H2", h2) << add_value("O2", o2) << add_value("H2O", h2o);
  }
}

void InitializeProbeOutput(const ProgramOptions& po) {
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    auto init_probe_output = [](std::string file_name, std::string header,
                                auto filter, auto formatter) {
      using text_sink = boost::log::sinks::synchronous_sink<
          boost::log::sinks::text_ostream_backend>;
      boost::filesystem::path path{file_name};
      boost::filesystem::path dir = path.parent_path();
      if (!boost::filesystem::exists(dir)) {
        boost::filesystem::create_directory(dir);
      }
      boost::shared_ptr<text_sink> file = boost::make_shared<text_sink>();
      boost::shared_ptr stream = boost::make_shared<std::ofstream>(file_name);
      if (!*stream) {
        throw std::runtime_error(
            fmt::format("Could not create '{}'.", file_name));
      }
      *stream << header;
      file->locked_backend()->add_stream(stream);
      file->set_formatter(formatter);
      file->set_filter(filter);
      boost::log::core::get()->add_sink(file);
    };
    init_probe_output(
        fmt::format("{}/ProbesTube.dat", po.output_directory),
        fmt::format("{:24}{:24}{:24}{:24}{:24}{:24}{:24}{:24}{:24}{:24}\n",
                    "Time", "X", "Density", "VelocityX", "Temperature",
                    "Pressure", "SpeedOfSound", "H2", "O2", "H2O"),
        [](const boost::log::attribute_value_set& attrs) -> bool {
          if (attrs.count("Channel")) {
            auto channel = attrs["Channel"].extract<std::string>();
            return channel && channel.get() == "ProbesTube";
          }
          return false;
        },
        [](const boost::log::record_view& rec,
           boost::log::formatting_ostream& stream) {
          namespace log = boost::log;
          double t = log::extract<double>("Time", rec).get();
          double x = log::extract<double>("X", rec).get();
          double rho = log::extract<double>("Density", rec).get();
          double u = log::extract<double>("VelocityX", rec).get();
          double T = log::extract<double>("Temperature", rec).get();
          double p = log::extract<double>("Pressure", rec).get();
          double a = log::extract<double>("SpeedOfSound", rec).get();
          double h2 = log::extract<double>("H2", rec).get();
          double o2 = log::extract<double>("O2", rec).get();
          double h2o = log::extract<double>("H2O", rec).get();
          stream << fmt::format(
              "{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}"
              "{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}",
              t, x, rho, u, T, p, a, h2, o2, h2o);
        });
    init_probe_output(
        fmt::format("{}/ProbesPlenum.dat", po.output_directory),
                      fmt::format("{:24}{:24}{:24}{:24}{:24}{:24}{:24}{:24}{:24}{:24}{:24}{:24}\n",
                    "Time", "X", "Density", "VelocityX", "VelocityY",
                    "VelocityZ", "Temperature", "Pressure", "SpeedOfSound",
                    "H2", "O2", "H2O"),
        [](const boost::log::attribute_value_set& attrs) -> bool {
          if (attrs.count("Channel")) {
            auto channel = attrs["Channel"].extract<std::string>();
            return channel && channel.get() == "ProbesPlenum";
          }
          return false;
        },
        [](const boost::log::record_view& rec,
           boost::log::formatting_ostream& stream) {
          namespace log = boost::log;
          double t = log::extract<double>("Time", rec).get();
          double x = log::extract<double>("X", rec).get();
          double rho = log::extract<double>("Density", rec).get();
          double u = log::extract<double>("VelocityX", rec).get();
          double v = log::extract<double>("VelocityY", rec).get();
          double w = log::extract<double>("VelocityZ", rec).get();
          double T = log::extract<double>("Temperature", rec).get();
          double p = log::extract<double>("Pressure", rec).get();
          double a = log::extract<double>("SpeedOfSound", rec).get();
          double h2 = log::extract<double>("H2", rec).get();
          double o2 = log::extract<double>("O2", rec).get();
          double h2o = log::extract<double>("H2O", rec).get();
          stream << fmt::format(
              "{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}"
              "{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}"
              "{:< 24.15g}{:< 24.15g}",
              t, x, rho, u, v, w, T, p, a, h2, o2, h2o);
        });
  }
}

void MyMain(const boost::program_options::variables_map& vm) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::Burke2012 mechanism{};

  const ProgramOptions po(vm);
  InitializeProbeOutput(po);

  auto plenum = MakePlenumSolver(po, mechanism, vm);
  auto [tube, valve] = MakeTubeSolver(po, mechanism, vm);
  auto valve_state = valve.GetSharedState();

  ::amrex::RealBox inlet{{-0.1, -r_tube, -r_tube}, {0.05, +r_tube, +r_tube}};

  fub::amrex::BlockConnection connection;
  connection.direction = fub::Direction::X;
  connection.side = 0;
  connection.plenum.id = 0;
  connection.plenum.mirror_box = 
      plenum.GetGriddingAlgorithm()->GetPatchHierarchy().GetGeometry(0).Domain();  
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

  std::string checkpoint{};
  if (vm.count("grid.checkpoint")) {
    checkpoint = vm["grid.checkpoint"].as<std::string>();
  }
  if (!checkpoint.empty()) {
    MPI_Comm comm = context.GetMpiCommunicator();
    std::string input = ReadAndBroadcastFile(checkpoint + "/Ignition", comm);
    std::istringstream ifs(input);
    {
      boost::archive::text_iarchive ia(ifs);
      std::vector<fub::Duration> last_ignitions;
      ia >> last_ignitions;
      ignition.SetLastIgnitionTimePoints(last_ignitions);
    }
    input =
        ReadAndBroadcastFile(checkpoint + "/Valve", comm);
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

  std::vector<double> xprobes = po.x_probes;
  std::sort(xprobes.begin(), xprobes.end());
  auto lb =
      std::lower_bound(xprobes.begin(), xprobes.end(), po.plenum_xbox.lo(0));
  const int n_tube_probes =
      static_cast<int>(std::distance(xprobes.begin(), lb));
  const int n_plenum_probes = static_cast<int>(xprobes.size()) - n_tube_probes;
  FUB_ASSERT(n_tube_probes + n_plenum_probes == int(xprobes.size()));

  std::vector<double> tube_probes_buffer(n_tube_probes * 3);
  ProbesView<double> tube_probes(tube_probes_buffer.data(), n_tube_probes);
  std::for_each(xprobes.begin(), lb, [&, i = 0](double xpos) mutable {
    tube_probes(0, i) = xpos;
    i = i + 1;
  });

  const double r_inner = r_tube;
  const double d_tube = 2 * r_tube;
  const double r_tube_center = 0.5 * (r_inner + r_inner + 2*po.plenum_jump + d_tube);

  std::vector<double> probes_buffer(n_plenum_probes * 3);
  ProbesView<double> plenum_probes(probes_buffer.data(), n_plenum_probes);
  std::for_each(lb, xprobes.end(), [&, i = 0](double xpos) mutable {
    plenum_probes(0, i) = xpos;
    plenum_probes(1, i) = r_tube_center;
    i = i + 1;
  });

  std::vector<double> slice_xs = {-3e-3, 3e-3, 0.1, 0.2, 0.245, 0.25, 0.31};
  std::vector<::amrex::Box> output_boxes{};
  output_boxes.reserve(slice_xs.size() + 1);

  std::transform(
      slice_xs.begin(), slice_xs.end(), std::back_inserter(output_boxes),
      [&](double x0) {
        const auto& plenum =
            context.GetGriddingAlgorithm()->GetPlena()[0]->GetPatchHierarchy();
        const int finest_level = plenum.GetNumberOfLevels() - 1;
        const ::amrex::Geometry& geom = plenum.GetGeometry(finest_level);
        const ::amrex::RealBox& probDomain = geom.ProbDomain();
        const double xlo[] = {x0, probDomain.lo(1), probDomain.lo(2)};
        const double* xhi = probDomain.hi();
        const ::amrex::RealBox slice_x(xlo, xhi);
        ::amrex::Box slice_box = BoxWhichContains(slice_x, geom);
        slice_box.setBig(0, slice_box.smallEnd(0));
        return slice_box;
      });

  output_boxes.push_back([&](double y0) {
    const auto& plenum =
        context.GetGriddingAlgorithm()->GetPlena()[0]->GetPatchHierarchy();
    const int finest_level = plenum.GetNumberOfLevels() - 1;
    const ::amrex::Geometry& geom = plenum.GetGeometry(finest_level);
    const ::amrex::RealBox& probDomain = geom.ProbDomain();
    const double xlo[] = {probDomain.lo(0), y0, probDomain.lo(2)};
    const double* xhi = probDomain.hi();
    const ::amrex::RealBox slice_x(xlo, xhi);
    ::amrex::Box slice_box = BoxWhichContains(slice_x, geom);
    slice_box.setBig(1, slice_box.smallEnd(1));
    return slice_box;
  }(0.0));

  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);
  MPI_Comm comm = context.GetMpiCommunicator();
  int rank = -1;
  MPI_Comm_rank(comm, &rank);
  auto output =
      [&](std::shared_ptr<fub::amrex::MultiBlockGriddingAlgorithm> gridding,
          auto cycle, auto timepoint, int output_num) {
        BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", timepoint.count());
        if (output_num >= 0) {
          BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "ProbesTube");
          LogTubeProbes(log, tube_probes,
                        gridding->GetTubes()[0]->GetPatchHierarchy(), comm);
        }
        if (output_num >= 0) {
          BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "ProbesPlenum");
          LogPlenumProbes(log, plenum_probes,
                          gridding->GetPlena()[0]->GetPatchHierarchy(), comm);
        }
        if (output_num == 0) {
          BOOST_LOG(log) << fmt::format("Write Plotfiles to '{}/Plotfiles'.",
                                        base_name);
          std::string name =
              fmt::format("{}/Plotfiles/Plenum/plt{:05}", base_name, cycle);
          fub::amrex::cutcell::WritePlotFile(
              name, gridding->GetPlena()[0]->GetPatchHierarchy(),
              plenum_equation);
          name = fmt::format("{}/Plotfiles/Tube/plt{:05}", base_name, cycle);
          fub::amrex::WritePlotFile(
              name, gridding->GetTubes()[0]->GetPatchHierarchy(),
              tube_equation);
          name = fmt::format("{}/Checkpoint/{:05}", base_name, cycle);
          WriteCheckpoint(name, *gridding, valve_state, rank, ignition);
        }
        //
        // Ouput MATLAB files on each output interval
        //
        if (output_num < 2) {
          auto tubes = gridding->GetTubes();
          int k = 0;
          for (auto& tube : tubes) {
            std::string name =
                fmt::format("{}/Matlab/Tube_{}/plt{:05}.dat", base_name, k, cycle);
            fub::amrex::WriteTubeData(name, tube->GetPatchHierarchy(),
                                      tube_equation, timepoint, cycle, comm);
            k = k + 1;
          }
          k = 0;
          std::for_each(output_boxes.begin(), output_boxes.end() - 1,
                        [&](const ::amrex::Box& out_box) {
                          std::string name =
                              fmt::format("{}/Matlab/Plenum_x{}/plt{:05}.dat",
                                          base_name, k, cycle);
                          auto& plenum = gridding->GetPlena()[0];
                          fub::amrex::cutcell::Write2Dfrom3D(
                              name, plenum->GetPatchHierarchy(), out_box, plenum_equation,
                              timepoint, cycle, comm);
                          k = k + 1;
                        });
          const ::amrex::Box out_box = output_boxes.back();
          std::string name =
              fmt::format("{}/Matlab/Plenum_y0/plt{:05}.dat", base_name, cycle);
          auto& plenum = gridding->GetPlena()[0];
          fub::amrex::cutcell::Write2Dfrom3D(name, plenum->GetPatchHierarchy(),
                                            out_box, plenum_equation, timepoint, cycle,
                                            comm);
        }
      };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint(), 0);
  fub::RunOptions run_options(vm);
  run_options.output_interval = std::vector<fub::Duration>{po.checkpoint_interval, run_options.output_interval[0]};
  run_options.output_frequency = std::vector{0, 0, 1};
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}
