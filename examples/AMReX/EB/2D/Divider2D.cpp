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

#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>

#include <iostream>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>

#include <boost/container/pmr/vector.hpp>
#include <boost/filesystem.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

struct DividerOptions {
  std::vector<std::string> wall_filenames{};
  double mach_number{1.77};
  std::array<double, 2> x_range{0.0, 0.042};
  std::array<double, 2> y_range{-0.016, +0.026};
  std::array<int, 2> n_cells{200, 200};
  int n_level{1};

  DividerOptions() = default;

  DividerOptions(const std::map<std::string, pybind11::object>& map) {
    wall_filenames = fub::GetOptionOr(map, "wall_filenames", wall_filenames);
    x_range = fub::GetOptionOr(map, "x_range", x_range);
    y_range = fub::GetOptionOr(map, "y_range", y_range);
    n_cells = fub::GetOptionOr(map, "n_cells", n_cells);
    mach_number = fub::GetOptionOr(map, "Mach_number", mach_number);
    n_level = fub::GetOptionOr(map, "n_level", n_level);
  }

  template <typename Logger> void Print(Logger& log) {
    BOOST_LOG(log) << "Divider Options:"
                   << "\n  - mach_number = " << mach_number << " [-]"
                   << "\n  - x_range = {" << x_range[0] << ", " << x_range[1]
                   << "} [m]"
                   << "\n  - y_range = {" << y_range[0] << ", " << y_range[1]
                   << "} [m]"
                   << "\n  - n_cells = {" << n_cells[0] << ", " << n_cells[1]
                   << "} [-]"
                   << "\n  - n_level = " << n_level << " [-]"
                   << fmt::format("\n  - wall_filenames = {{{}}}",
                                  fmt::join(wall_filenames, ", "));
  }
};

std::optional<std::map<std::string, pybind11::object>>
ParseCommandLine(int argc, char** argv) {
  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);

  namespace po = boost::program_options;
  po::options_description desc;
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

  fub::RunOptions(options).Print(log);
  DividerOptions(options).Print(log);

  return options;
}

fub::Polygon ReadPolygonData(std::istream& input) {
  std::string line{};
  namespace pmr = boost::container::pmr;
  pmr::vector<double> xs{};
  pmr::vector<double> ys{};
  while (std::getline(input, line)) {
    double x{}, y{};
    std::istringstream linestream(line);
    linestream >> x >> y;
    if (linestream) {
      xs.push_back(x);
      ys.push_back(y);
    }
  }
  if (xs.empty() || ys.empty()) {
    throw std::invalid_argument{"Invalid Input File:: No data found!"};
  }
  const double x0 = xs.front();
  const double x1 = xs.back();
  const double y0 = ys.front();
  const double y1 = ys.back();
  if (x0 != x1 || y0 != y1) {
    throw std::invalid_argument{
        "Invalid Input File: First and last entries are not the same point."};
  }
  return fub::Polygon(std::move(xs), std::move(ys));
}

template <typename Geometry>
struct ShockMachnumber
    : fub::amrex::cutcell::RiemannProblem<fub::PerfectGas<2>, Geometry> {
  static fub::PerfectGas<2>::Complete
  ComputePreShockState(fub::PerfectGas<2>& equation,
                       const fub::PerfectGas<2>::Complete& post_shock,
                       double M_S, const fub::Array<double, 2, 1>& normal) {
    const double g = equation.gamma;
    const double gp1 = g + 1.0;
    const double gm1 = g - 1.0;
    const double M_post = equation.Machnumber(post_shock);
    const double M2 = (M_post - M_S) * (M_post - M_S);

    const double rho_c = gp1 * M2 / (gm1 * M2 + 2.0);
    const double rho_post = post_shock.density;
    const double rho_pre = rho_post * rho_c;

    const double p_c = (2.0 * g * M2 - gm1) / gp1;
    const double p_post = post_shock.pressure;
    const double p_pre = p_post * p_c;

    const double t = rho_post / rho_pre;
    const double a_post = post_shock.speed_of_sound;
    const double u_S = M_S * a_post;
    const double u_post = equation.Velocity(post_shock).matrix().norm();
    const double u_pre = u_S * (1.0 - t) + u_post * t;

    return equation.CompleteFromPrim(rho_pre, u_pre * normal, p_pre);
  }

  ShockMachnumber(fub::PerfectGas<2> equation, Geometry geometry,
                  const fub::PerfectGas<2>::Complete& post_shock, double M_S,
                  const fub::Array<double, 2, 1>& normal)
      : fub::amrex::cutcell::RiemannProblem<fub::PerfectGas<2>, Geometry>(
            equation, std::move(geometry),
            ComputePreShockState(equation, post_shock, M_S, normal),
            post_shock) {}
};

// template <typename Logger>
// void WriteCheckpoint(Logger& log, const std::string& path,
//                      const fub::amrex::cutcell::PatchHierarchy& hierarchy) {
//   BOOST_LOG(log) << "Write Checkpoint File to '" << path << "'.\n";
//   fub::amrex::cutcell::WriteCheckpointFile(path, hierarchy);
// }

struct CheckpointOutput : fub::OutputAtFrequencyOrInterval<
                              fub::amrex::cutcell::GriddingAlgorithm> {
  CheckpointOutput(
      const fub::ProgramOptions& options)
      : OutputAtFrequencyOrInterval(options) {
    directory_ = fub::GetOptionOr(options, "directory", directory_);
    fub::SeverityLogger log = fub::GetInfoLogger();
    BOOST_LOG(log) << "CheckpointOutput:";
    OutputAtFrequencyOrInterval::Print(log);
    BOOST_LOG(log) << fmt::format(" - directory = '{}'", directory_);
  }

  void
  operator()(const fub::amrex::cutcell::GriddingAlgorithm& grid) override {
    boost::log::sources::severity_logger<boost::log::trivial::severity_level>
        log(boost::log::keywords::severity = boost::log::trivial::info);
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", grid.GetTimePoint().count());
    BOOST_LOG(log) << fmt::format("Write checkpoint to '{}'.", directory_);
    fub::amrex::cutcell::WriteCheckpointFile(directory_, grid.GetPatchHierarchy());
  }

  std::string directory_{"./Divider/"};
};

void MyMain(const std::map<std::string, pybind11::object>& vm) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);

  DividerOptions options(vm);

  const std::array<int, 2> n_cells = options.n_cells;
  const std::array<double, 2> xlower{options.x_range[0], options.y_range[0]};
  const std::array<double, 2> xupper{options.x_range[1], options.y_range[1]};
  const std::array<int, 2> periodicity{0, 0};

  fub::PerfectGas<2> equation;

  using namespace fub::amrex::cutcell;
  using State = fub::Complete<fub::PerfectGas<2>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.22;
  cons.momentum << 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> post_shock_state;
  fub::CompleteFromCons(equation, post_shock_state, cons);

  const double shock_mach_number = options.mach_number;
  const fub::Array<double, 2, 1> normal{1.0, 0.0};

  ShockMachnumber<fub::Halfspace> initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.01), post_shock_state,
      shock_mach_number, normal);

  const State& pre_shock_state = initial_data.left_;
  BOOST_LOG(log) << "Post-Shock-State:\n"
                 << "\tdensity: " << post_shock_state.density << " [kg/m^3]\n"
                 << "\tvelocity: "
                 << equation.Velocity(post_shock_state).transpose()
                 << " [m/s]\n"
                 << "\tpressure: " << post_shock_state.pressure << " [Pa]";

  BOOST_LOG(log) << "Calculated Pre-Shock-State:\n"
                 << "\tdensity: " << pre_shock_state.density << " [kg/m^3]\n"
                 << "\tvelocity: "
                 << equation.Velocity(pre_shock_state).transpose() << " [m/s]\n"
                 << "\tpressure: " << pre_shock_state.pressure << " [Pa]";

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1}},
                              &xbox, -1, periodicity.data());

  const int n_level = options.n_level;

  std::vector<fub::PolymorphicGeometry<2>> geometries;
  std::transform(options.wall_filenames.begin(), options.wall_filenames.end(),
                 std::back_inserter(geometries),
                 [](const std::string& filename) {
                   std::ifstream ifs(filename);
                   return fub::PolymorphicGeometry<2>(ReadPolygonData(ifs));
                 });

  auto embedded_boundary =
      fub::amrex::Geometry(fub::PolymorphicUnion(geometries));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  PatchHierarchyOptions hier_opts{};
  hier_opts.ngrow_eb_level_set = 9;
  hier_opts.max_number_of_levels = n_level;
  hier_opts.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, coarse_geom, n_level);

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, hier_opts), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  fub::perfect_gas::HllemMethod<2> hllem_method{equation};
  fub::FluxMethod<fub::perfect_gas::MusclHancockPrim<2>> flux_method{equation};
  fub::KbnCutCellMethod cutcell_method(flux_method, hllem_method);
  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method, 2, 0));

  // fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
  fub::NoSubcycleSolver solver(std::move(level_integrator));

  using Plotfile = PlotfileOutput<fub::PerfectGas<2>>;
  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<WriteHdf5>("HDF5");
  factory.RegisterOutput<Plotfile>("Plotfiles", equation);
  factory.RegisterOutput<CheckpointOutput>("Checkpoint");
  fub::MultipleOutputs<GriddingAlgorithm> outputs(
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
