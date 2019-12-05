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

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard _(argc, argv);

  const std::array<int, 3> n_cells{64, 64, 64};
  const std::array<double, 3> xlower{-0.10, -0.15, -0.15};
  const std::array<double, 3> xupper{+0.20, +0.15, +0.15};
  const std::array<int, 3> periodicity{0, 0, 0};

  const int n_level = 3;

  auto embedded_boundary = amrex::EB2::makeIntersection(
      amrex::EB2::PlaneIF({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, false),
      amrex::EB2::CylinderIF(0.015, -1.0, 0, {1e6, 0.0, 0.0}, true));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::Burke2012 mechanism;
  fub::IdealGasMix<3> equation{mechanism};
  //  fub::PerfectGas<3> equation{};

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(
      amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1}}, &xbox,
      -1, periodicity.data());

  using namespace fub::amrex::cutcell;

  PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, coarse_geom, n_level);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  equation.GetReactor().SetMoleFractions("N2:80,O2:20");
  equation.GetReactor().SetTemperature(300.0);
  equation.GetReactor().SetDensity(1.22);
  fub::Complete<fub::IdealGasMix<3>> right(equation);
  equation.CompleteFromReactor(right);

  equation.GetReactor().SetMoleFractions("N2:80,O2:20");
  equation.GetReactor().SetTemperature(500.0);
  equation.GetReactor().SetDensity(3.15);
  fub::Complete<fub::IdealGasMix<3>> left(equation);
  equation.CompleteFromReactor(left, {400.0, 0.0, 0.0});

  // fub::Conservative<fub::PerfectGas<3>> cons;
  //  cons.density = 1.0;
  //  cons.momentum << 0.0, 0.0, 0.0;
  //  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  //  fub::Complete<fub::PerfectGas<3>> right;
  //  fub::CompleteFromCons(equation, right, cons);
  //
  //  cons.energy *= 4;
  //  fub::Complete<fub::PerfectGas<3>> left;
  //  fub::CompleteFromCons(equation, left, cons);

  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04),
                              left, right);

  //  using Complete = fub::Complete<fub::PerfectGas<3>>;
  using Complete = fub::Complete<fub::IdealGasMix<3>>;
  GradientDetector gradients{equation, std::pair{&Complete::pressure, 0.05},
                             std::pair{&Complete::density, 0.05}};

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

  //  fub::EinfeldtSignalVelocities<fub::PerfectGas<3>> signals{};
  fub::EinfeldtSignalVelocities<fub::IdealGasMix<3>> signals{};
  fub::HllMethod hll_method{equation, signals};
  //  fub::MusclHancockMethod flux_method(equation, hll_method);
  fub::ideal_gas::MusclHancockPrimMethod<3> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{fub::execution::simd, cutcell_method},
                          TimeIntegrator{},
                          Reconstruction{fub::execution::simd, equation}};

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<3>, IntegratorContext(gridding, method));

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  std::string base_name = "LinearShock3d/";

  using namespace std::literals::chrono_literals;
  fub::AnyOutput<GriddingAlgorithm> output(
      {}, {0.0000125s}, [&](const GriddingAlgorithm& gridding) {
        std::ptrdiff_t cycle = gridding.GetCycles();
        std::string name = fmt::format("{}plt{:05}", base_name, cycle);
        amrex::Print() << "Start output to '" << name << "'.\n";
        WritePlotFile(name, gridding.GetPatchHierarchy(), equation);
        amrex::Print() << "Finished output to '" << name << "'.\n";
      });

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.final_time = 0.002s;
  run_options.cfl = 0.8;
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}
