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

#include <boost/filesystem.hpp>

#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Union.H>

#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO |
                         _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW |
                         _MM_MASK_INVALID);

  static_assert(AMREX_SPACEDIM == 3);

  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  const fub::amrex::ScopeGuard _(argc, argv);

  const int r = 3;

  const std::array<int, 3> n_cells{8 * r * 12, 8 * r * 3, 8 * r * 2};
  const std::array<double, 3> xlower{0.0, -0.01, -0.02};
  const std::array<double, 3> xupper{0.23, 0.05, 0.02};
  const std::array<int, 3> periodicity{0, 0, 0};

  const int n_level = 1;

  Eigen::Matrix<double, 3, Eigen::Dynamic> centerline =
      fub::ReadPointsFromFile("centerline.txt");

  const double radius = 0.5 * 0.0125;
  auto embedded_boundary = amrex::EB2::makeComplement(amrex::EB2::makeUnion(
      amrex::EB2::CylinderIF(radius, 1.0, 0, {-0.1, radius, 0.0}, false),
      fub::ExpandTube(centerline, radius)));

  auto shop = amrex::EB2::makeShop(embedded_boundary);

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(
      amrex::Box{
          {}, {AMREX_D_DECL(n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1)}},
      &xbox, -1, periodicity.data());

  fub::PerfectGas<3> equation;

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = xbox;
  geometry.periodicity = periodicity;

  using namespace fub::amrex::cutcell;

  PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, coarse_geom, n_level);

  fub::Conservative<fub::PerfectGas<3>> cons;
  cons.density = 1.22;
  cons.momentum << 0.0, 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<3>> right;
  fub::CompleteFromCons(equation, right, cons);

  // cons.density = 3.15736;
  // cons.momentum << 1258.31, 0.0, 0.0;
  // cons.energy = 416595.0 * equation.gamma_minus_1_inv +
  //               0.5 * cons.momentum[0] * cons.momentum[0] / cons.density;
  cons.energy *= 8.0;
  fub::Complete<fub::PerfectGas<3>> left;
  fub::CompleteFromCons(equation, left, cons);

  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.04),
                              left, right);

  using State = fub::Complete<fub::PerfectGas<3>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1},
                                  TransmissiveBoundary{fub::Direction::Z, 0},
                                  TransmissiveBoundary{fub::Direction::Z, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, options), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  fub::EinfeldtSignalVelocities<fub::PerfectGas<3>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::MusclHancockMethod flux_method(equation, hll_method);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{fub::execution::seq, cutcell_method},
                          TimeIntegrator{},
                          Reconstruction{fub::execution::seq, equation}};
  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      equation, IntegratorContext(gridding, method)));
  std::string base_name = "Divider/";

  auto output = [&](const std::shared_ptr<GriddingAlgorithm>& gridding, std::ptrdiff_t cycle,
                    fub::Duration) {
    std::string name = fmt::format("{}{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    WritePlotFile(name, gridding->GetPatchHierarchy(), equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(), solver.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.final_time = 0.0005s;
  run_options.output_interval = 0.5 * 0.0000125s;
  run_options.cfl = 0.5 * 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
