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
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <iostream>

#include <xmmintrin.h>

using Coord = Eigen::Vector2d;

Coord OrthogonalTo(const Coord& x) { return Coord{x[1], -x[0]}; }

auto Wedge(const Coord& p1, const Coord& p2) {
  Coord p0{0.0, 0.1};
  Coord norm1 = OrthogonalTo(p1 - p0).normalized();
  Coord norm2 = OrthogonalTo(p2 - p0).normalized();
  amrex::EB2::PlaneIF plane1({p0[0], p0[1]}, {norm1[0], norm1[1]});
  amrex::EB2::PlaneIF plane2({p0[0], p0[1]}, {norm2[0], norm2[1]}, false);
  return amrex::EB2::makeComplement(
      amrex::EB2::makeIntersection(plane1, plane2));
}

int main(int argc, char** argv) {
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO |
                         _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW |
                         _MM_MASK_INVALID);

  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard _(argc, argv);

  const std::array<int, 2> n_cells{128, 128};
  const std::array<double, 2> xlower{-0.5, -1.0};
  const std::array<double, 2> xupper{+0.5, +1.0};
  amrex::RealBox xbox(xlower, xupper);
  const std::array<int, 2> periodicity{};

  amrex::Geometry coarse_geom(
      amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1}},
      &xbox, -1, periodicity.data());

  const int n_level = 1;

  auto embedded_boundary = Wedge({-1.0, +1.0}, {+1.0, +1.0});
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::PerfectGas<2> equation;

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  fub::amrex::cutcell::PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, coarse_geom, n_level);

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.0;
  cons.momentum << 0.0, 0.0; // +100.0, -100.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> right;
  fub::CompleteFromCons(equation, right, cons);

  cons.energy *= 5;
  fub::Complete<fub::PerfectGas<2>> left;
  fub::CompleteFromCons(equation, left, cons);

  fub::amrex::cutcell::RiemannProblem initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.4), left, right);

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                            TransmissiveBoundary{fub::Direction::X, 1},
                            TransmissiveBoundary{fub::Direction::Y, 0},
                            TransmissiveBoundary{fub::Direction::Y, 1}}};

  using State = fub::Complete<fub::PerfectGas<2>>;
  fub::amrex::cutcell::GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                                  std::pair{&State::density, 0.005}};

    std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, options), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  fub::MusclHancockMethod flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(std::move(flux_method));

  HyperbolicMethod method{FluxMethod{fub::execution::seq, cutcell_method},
                          TimeIntegrator{},
                          Reconstruction{fub::execution::seq, equation}};

  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      equation, fub::amrex::cutcell::IntegratorContext(gridding, method)));

  std::string base_name = "Wedge/";
  auto output = [&](const fub::amrex::cutcell::PatchHierarchy& hierarchy,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::cutcell::WritePlotFile(name, hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };
  output(solver.GetPatchHierarchy(), solver.GetCycles(), solver.GetTimePoint());

  using namespace std::literals::chrono_literals;
  fub::RunOptions run_options{};
  run_options.final_time = 1e-4s;
  run_options.output_interval = 1e-5s;
  //  run_options.output_frequency = 1;
  run_options.cfl = 0.5 * 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
