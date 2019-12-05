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
#include <AMReX_EB_LSCore.H>

#include <boost/container/pmr/vector.hpp>
#include <boost/filesystem.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

struct DividerOptions {
  double mach_number{1.77};
  std::array<double, 2> x_range{0.0, 0.042};
  std::array<double, 2> y_range{-0.016, +0.026};
  std::array<int, 2> n_cells{200, 200};
  int n_level{1};
  std::string output_directory{"ShockTube"};

  DividerOptions() = default;

  DividerOptions(const std::map<std::string, pybind11::object>& map) {
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
                   << "\n  - output_directory = '" << output_directory << "'";
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

template <typename Logger>
void WriteCheckpoint(Logger& log, const std::string& path,
                     const fub::amrex::cutcell::PatchHierarchy& hierarchy) {
  BOOST_LOG(log) << "Write Checkpoint File to '" << path << "'.\n";
  fub::amrex::cutcell::WriteCheckpointFile(path, hierarchy);
}

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

  fub::Burke2012 burke_2012{};
  fub::IdealGasMix<2> equation(burke_2012);

  using namespace fub::amrex::cutcell;
  using State = fub::Complete<fub::IdealGasMix<2>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

  fub::FlameMasterReactor& reactor = equation.GetReactor();
  reactor.SetMoleFractions("HE:15.0,N2:0.79,O2:0.21");
  reactor.SetTemperature(293.1);
  reactor.SetPressure(16.0 * 101325.0);
  fub::Complete<fub::IdealGasMix<2>> left(equation);
  equation.CompleteFromReactor(left);

  reactor.SetMoleFractions("N2:79,O2:21");
  reactor.SetTemperature(293.1);
  reactor.SetPressure(101325.0);
  fub::Complete<fub::IdealGasMix<2>> right(equation);
  equation.CompleteFromReactor(right);

  // const double shock_mach_number = options.mach_number;
  const fub::Array<double, 2, 1> normal{1.0, 0.0};

  fub::amrex::cutcell::RiemannProblem initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.25), left, right);

  // const State& pre_shock_state = initial_data.left_;
  BOOST_LOG(log) << "Post-Shock-State:\n"
                 << "\tdensity: " << left.density << " [kg/m^3]\n"
                 << "\tpressure: " << left.pressure << " [Pa]";

  BOOST_LOG(log) << "Pre-Shock-State:\n"
                 << "\tdensity: " << right.density << " [kg/m^3]\n"
                 << "\tpressure: " << right.pressure << " [Pa]";

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1}},
                              &xbox, -1, periodicity.data());

  const int n_level = options.n_level;

  auto MakePlenum2D = [](auto... points) {
    fub::Polygon::Vector xs{};
    fub::Polygon::Vector ys{};

    (xs.push_back(std::get<0>(points)), ...);
    (ys.push_back(std::get<1>(points)), ...);

    return fub::Polygon(std::move(xs), std::move(ys));
  };

  constexpr double d_big = 0.03508;
  constexpr double d_small = 0.01248;
  constexpr double d_smallest = 0.00850;
  constexpr double r_big = d_big * 0.5;
  constexpr double r_small = d_small * 0.5;
  constexpr double r_smallest = d_smallest * 0.5;

  // constexpr double start_big = 0.0;
  constexpr double end_big = 0.29828;
  constexpr double start_small = 0.3250;
  constexpr double end_small = 1.0200;
  constexpr double start_smallest = 1.0320;
  constexpr double end_smallest = 1.0570;

  auto wall = fub::amrex::Geometry(MakePlenum2D(
      std::pair{0.0, -r_big}, std::pair{end_big, -r_big},
      std::pair{start_small, -r_small}, std::pair{end_small, -r_small},
      std::pair{start_smallest, -r_smallest},
      std::pair{end_smallest, -r_smallest}, std::pair{end_smallest, -1.0},
      std::pair{-1.0, -1.0}, std::pair{-1.0, 1.0}, std::pair{end_smallest, 1.0},
      std::pair{end_smallest, r_smallest},
      std::pair{start_smallest, r_smallest}, std::pair{end_small, r_small},
      std::pair{start_small, r_small}, std::pair{end_big, r_big},
      std::pair{0.0, r_big}, std::pair{0.0, -r_big}));

  auto shop = amrex::EB2::makeShop(amrex::EB2::makeUnion(wall));

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  PatchHierarchyOptions hier_opts{};
  hier_opts.max_number_of_levels = n_level;
  hier_opts.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, coarse_geom, n_level);

  amrex::Box box = coarse_geom.Domain();
  box.setSmall(0, box.bigEnd(0) - 1);
  BoundarySet boundary_condition{
      {IsentropicPressureBoundary{"PressureExpansion", equation, box, 101325.0,
                                  fub::Direction::X, 1},
       TransmissiveBoundary{fub::Direction::X, 0},
       TransmissiveBoundary{fub::Direction::Y, 0},
       TransmissiveBoundary{fub::Direction::Y, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, hier_opts), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<2>> signals{};
  fub::HllMethod hll_method{equation, signals};
  // fub::MusclHancockMethod flux_method(equation, hll_method);
  fub::ideal_gas::MusclHancockPrimMethod<2> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{fub::execution::simd, cutcell_method},
                          TimeIntegrator{},
                          Reconstruction{fub::execution::simd, equation}};

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method));

  fub::SubcycleFineFirstSolver solver(level_integrator);

  std::string base_name = options.output_directory;

  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<WriteHdf5>("HDF5");
  fub::MultipleOutputs<GriddingAlgorithm> outputs(
      std::move(factory),
      fub::ToMap(fub::GetOptionOr(vm, "output", pybind11::dict{})));

  using namespace std::literals::chrono_literals;
  outputs.AddOutput(std::make_unique<PlotfileOutput<fub::IdealGasMix<2>>>(
      std::vector<std::ptrdiff_t>{}, std::vector<fub::Duration>{1e-6s}, equation,
      base_name + "/Plotfiles"));

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
