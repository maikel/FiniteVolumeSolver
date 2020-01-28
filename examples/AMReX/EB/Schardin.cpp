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
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB_LSCore.H>

static_assert(AMREX_SPACEDIM == 2);

using Coord = Eigen::Vector2d;

Coord OrthogonalTo(const Coord& x) { return Coord{x[1], -x[0]}; }

auto Triangle(const Coord& p1, const Coord& p2, const Coord& p3) {
  Coord norm1 = OrthogonalTo(p2 - p1).normalized();
  Coord norm2 = OrthogonalTo(p3 - p2).normalized();
  Coord norm3 = OrthogonalTo(p1 - p3).normalized();
  amrex::EB2::PlaneIF plane1({p1[0], p1[1]}, {norm1[0], norm1[1]});
  amrex::EB2::PlaneIF plane2({p2[0], p2[1]}, {norm2[0], norm2[1]});
  amrex::EB2::PlaneIF plane3({p3[0], p3[1]}, {norm3[0], norm3[1]});
  return amrex::EB2::makeIntersection(plane1, plane2, plane3);
}

int main() {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard _{};
  fub::InitializeLogging(MPI_COMM_WORLD);

  const std::array<int, 2> n_cells{16 * 15, 16 * 10};
  const std::array<double, 2> xlower{0.0, 0.0};
  const std::array<double, 2> xupper{+0.15001, +0.10};
  amrex::RealBox xbox(xlower, xupper);
  const std::array<int, 2> periodicity{};

  amrex::Geometry coarse_geom(amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1}},
                              &xbox, -1, periodicity.data());

  const int n_level = 3;

  const int scratch_gcw = 8;
  const int flux_gcw = 6;

  auto embedded_boundary =
      Triangle({0.02, 0.05}, {0.05, 0.0655}, {0.05, 0.0345});
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::PerfectGas<2> equation;

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  using namespace fub::amrex::cutcell;

  PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces = fub::amrex::cutcell::MakeIndexSpaces(
      shop, coarse_geom, n_level, scratch_gcw);
  options.ngrow_eb_level_set = scratch_gcw;

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.0;
  cons.momentum << 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> right;
  fub::CompleteFromCons(equation, right, cons);

  cons.energy *= 5;
  fub::Complete<fub::PerfectGas<2>> left;
  fub::CompleteFromCons(equation, left, cons);

  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.015),
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
      TagAllOf(TagCutCells(), gradients, TagBuffer(2)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  fub::EinfeldtSignalVelocities<fub::PerfectGas<2>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::MusclHancockMethod muscl_method{equation, hll_method};
  fub::KbnCutCellMethod cutcell_method(muscl_method, hll_method);

  HyperbolicMethod method{FluxMethod{fub::execution::simd, cutcell_method},
                          TimeIntegrator{},
                          Reconstruction{fub::execution::simd, equation}};

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method, scratch_gcw, flux_gcw),
      fub::StrangSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  std::string base_name = "Schardin/";
  using namespace std::literals::chrono_literals;
  fub::MultipleOutputs<GriddingAlgorithm> output{};
  output.AddOutput(fub::MakeOutput<GriddingAlgorithm>(
      {}, {3e-4s / 180.}, PlotfileOutput(equation, base_name)));
  output.AddOutput(std::make_unique<fub::CounterOutput<GriddingAlgorithm>>(
      solver.GetContext().registry_, wall_time_reference,
      std::vector<std::ptrdiff_t>{50}, std::vector<fub::Duration>{}));

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.final_time = 3e-4s;
  run_options.cfl = 0.4;
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}
