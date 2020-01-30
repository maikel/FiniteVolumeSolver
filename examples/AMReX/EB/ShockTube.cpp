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

  DividerOptions() = default;

  DividerOptions(const fub::ProgramOptions& map) {
    mach_number = fub::GetOptionOr(map, "Mach_number", mach_number);
  }

  template <typename Logger> void Print(Logger& log) {
    BOOST_LOG(log) << "Divider Options:"
                   << "\n  - mach_number = " << mach_number << " [-]";
  }
};

template <typename Logger>
void WriteCheckpoint(Logger& log, const std::string& path,
                     const fub::amrex::cutcell::PatchHierarchy& hierarchy) {
  BOOST_LOG(log) << "Write Checkpoint File to '" << path << "'.\n";
  fub::amrex::cutcell::WriteCheckpointFile(path, hierarchy);
}

void MyMain(const fub::ProgramOptions& vm) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};

  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);

  DividerOptions options(vm);

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

  fub::amrex::CartesianGridGeometry grid_geometry(
      fub::GetOptions(vm, "GridGeometry"));

  PatchHierarchyOptions hierarchy_options =
      fub::GetOptions(vm, "PatchHierarchy");
  hierarchy_options.index_spaces =
      MakeIndexSpaces(shop, grid_geometry, hierarchy_options);

  IsentropicPressureBoundaryOptions boundary_options =
      fub::GetOptions(vm, "IsentropicPressureBoundary");
  BoundarySet boundary_condition{
      {IsentropicPressureBoundary{equation, boundary_options},
       TransmissiveBoundary{fub::Direction::X, 0},
       TransmissiveBoundary{fub::Direction::Y, 0},
       TransmissiveBoundary{fub::Direction::Y, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, grid_geometry, hierarchy_options), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(2)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<2>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::ideal_gas::MusclHancockPrimMethod<2> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method, 4, 2),
      fub::GodunovSplitting{});

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  std::string base_name = "ShockTube";

  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<WriteHdf5>("HDF5");
  fub::MultipleOutputs<GriddingAlgorithm> outputs(
      std::move(factory), fub::GetOptions(vm, "Output"));

  using namespace std::literals::chrono_literals;
  outputs.AddOutput(std::make_unique<PlotfileOutput<fub::IdealGasMix<2>>>(
      std::vector<std::ptrdiff_t>{}, std::vector<fub::Duration>{1e-6s},
      equation, base_name + "/Plotfiles"));

  outputs(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, fub::GetOptions(vm, "RunOptions"),
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
