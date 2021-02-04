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

#include "fub/cutcell_method/MyStabilisation.hpp"

#include "fub/AMReX/cutcell/boundary_condition/ConstantBoundary.hpp"

#include <AMReX_EB2_IF_Plane.H>

Eigen::Vector2d OrthogonalTo(const Eigen::Vector2d& x) {
  return Eigen::Vector2d{x[1], -x[0]};
}

auto Plane(const Eigen::Vector2d& p1) {
  Eigen::Vector2d p0{0.0, 0.0};
  Eigen::Vector2d norm1 = OrthogonalTo(p1 - p0).normalized();
  amrex::EB2::PlaneIF plane1({p0[0], p0[1]}, {norm1[0], norm1[1]}, false);
  return plane1;
}

void MyMain(const fub::ProgramOptions& opts) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard guard{};
  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::amrex::CartesianGridGeometry geometry =
      fub::GetOptions(opts, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  geometry.Print(log);

  fub::amrex::cutcell::PatchHierarchyOptions hier_opts =
      fub::GetOptions(opts, "PatchHierarchy");
  BOOST_LOG(log) << "PatchHierarchy:";
  hier_opts.Print(log);

  BOOST_LOG(log) << "Compute EB level set data...";
  auto embedded_boundary = Plane({-1.0, +1.0});
  auto shop = amrex::EB2::makeShop(embedded_boundary);
  hier_opts.index_spaces = MakeIndexSpaces(shop, geometry, hier_opts);

  fub::PerfectGas<2> equation{};

  Eigen::Vector2d v{0.1, -0.1};
  fub::Complete<fub::PerfectGas<2>> right = equation.CompleteFromPrim(1.0, v, 1.0);
  fub::Complete<fub::PerfectGas<2>> left = equation.CompleteFromPrim(1.0, v, 5.0);

  using namespace fub::amrex::cutcell;
  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, -1.0, 0.0}, 0.4),
                              left, right);

  using State = fub::Complete<fub::PerfectGas<2>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1}}};

  BOOST_LOG(log) << "Initialize Data...";

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, hier_opts), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  // fub::EinfeldtSignalVelocities<fub::PerfectGas<2>> signals{};
  // fub::HllMethod hll_method{equation, signals};
  fub::perfect_gas::HllemMethod<fub::PerfectGas<2>> hllem_method{equation};
  fub::MusclHancockMethod flux_method(equation, hllem_method);
  // fub::FluxMethod<fub::perfect_gas::MusclHancockPrim<2>> flux_method{equation};
  // fub::KbnCutCellMethod cutcell_method(flux_method, hllem_method);
  fub::MyCutCellMethod cutcell_method(equation, flux_method);
  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator2{},
                          Reconstruction{equation}};

  BOOST_LOG(log) << "Create Integrator Context...";

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method, 3, 0),
      fub::GodunovSplitting());

  fub::NoSubcycleSolver solver(std::move(level_integrator));

  using namespace std::literals::chrono_literals;
  using Plotfile = fub::amrex::cutcell::PlotfileOutput<fub::PerfectGas<2>>;
  using CounterOutput =
      fub::CounterOutput<GriddingAlgorithm, std::chrono::milliseconds>;
  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<Plotfile>("Plotfile", equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<fub::amrex::cutcell::WriteHdf5>("HDF5");
  factory.RegisterOutput<fub::amrex::cutcell::DebugOutput>(
      "DebugOutput",
      solver.GetGriddingAlgorithm()->GetPatchHierarchy().GetDebugStorage());
  fub::MultipleOutputs<GriddingAlgorithm> output(
      std::move(factory), fub::GetOptions(opts, "Output"));

  output(*solver.GetGriddingAlgorithm());

  fub::RunOptions run_options = fub::GetOptions(opts, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);

  BOOST_LOG(log) << "Run Simulation...";
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
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