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

std::shared_ptr<fub::amrex::cutcell::GriddingAlgorithm>
MakeGriddingAlgorithm(const fub::PerfectGas<2>& equation) {
  const std::array<int, 2> n_cells{200, 200};
  const std::array<double, 2> xlower{-1.0, -1.0};
  const std::array<double, 2> xupper{+1.0, +1.0};

  auto embedded_boundary = Plane({-1.0, +1.0});
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  amrex::RealBox xbox(xlower, xupper);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = xbox;

  const amrex::Box domain{{}, {n_cells[0] - 1, n_cells[1] - 1}};
  amrex::Geometry coarse_geom(domain, &xbox, -1, geometry.periodicity.data());

  const int n_level = 2;

  using namespace fub::amrex::cutcell;

  PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces = MakeIndexSpaces(shop, coarse_geom, n_level);

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.0;
  cons.momentum << +1.0, -1.0;
  cons.energy = equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> right;
  fub::CompleteFromCons(equation, right, cons);

  cons.energy *= 5;
  fub::Complete<fub::PerfectGas<2>> left;
  fub::CompleteFromCons(equation, left, cons);

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

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, options), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  return gridding;
}

auto MakeSolver(const fub::PerfectGas<2>& equation) {
  using namespace fub::amrex::cutcell;
  std::shared_ptr<GriddingAlgorithm> gridding = MakeGriddingAlgorithm(equation);

   fub::EinfeldtSignalVelocities<fub::PerfectGas<2>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::perfect_gas::HllemMethod<2> hllem_method{equation};
//  fub::MusclHancockMethod flux_method(equation, hllem_method);
  fub::FluxMethod<fub::perfect_gas::MusclHancockPrim<2>> flux_method{equation};
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);
  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method, 4, 2),
      fub::StrangSplitting());

  return fub::NoSubcycleSolver(std::move(level_integrator));
}

int main() {
  // fub::EnableFloatingPointExceptions();
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard amrex_scope_guard{};
  fub::InitializeLogging(MPI_COMM_WORLD);
  fub::PerfectGas<2> equation{};
  fub::NoSubcycleSolver solver = MakeSolver(equation);

  std::string base_name = "Ramp/";

  using namespace fub::amrex::cutcell;
  using namespace std::literals::chrono_literals;
  fub::MultipleOutputs<GriddingAlgorithm> output{};
  output.AddOutput(fub::MakeOutput<GriddingAlgorithm>(
      {1}, {1.0s / 180.}, PlotfileOutput(equation, base_name)));
  output.AddOutput(std::make_unique<fub::CounterOutput<GriddingAlgorithm>>(
      wall_time_reference, std::vector<std::ptrdiff_t>{2},
      std::vector<fub::Duration>{}));
  fub::RunOptions run_options{};
  run_options.final_time = 1s;
  run_options.cfl = 0.4;
  output(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}