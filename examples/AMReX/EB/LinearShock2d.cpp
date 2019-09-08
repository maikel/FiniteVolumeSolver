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
// OUT OF OR ILN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "fub/AMReX.hpp"
#include "fub/AMReX_CutCell.hpp"
#include "fub/Solver.hpp"

#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <iostream>

#include <xmmintrin.h>

auto Rectangle(const std::array<double, 2>& lower,
               const std::array<double, 2>& upper) {
  amrex::EB2::PlaneIF lower_x({lower[0], lower[1]}, {0, +1});
  amrex::EB2::PlaneIF lower_y({lower[0], lower[1]}, {+1, 0});
  amrex::EB2::PlaneIF upper_x({upper[0], upper[1]}, {0, -1});
  amrex::EB2::PlaneIF upper_y({upper[0], upper[1]}, {-1, 0});
  return amrex::EB2::makeIntersection(lower_x, lower_y, upper_x, upper_y);
}

int main(int argc, char** argv) {
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO |
                         _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW |
                         _MM_MASK_INVALID);

  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard _(argc, argv);

  const std::array<int, 2> n_cells{128, 128};
  const std::array<double, 2> xlower{-0.10, -0.6};
  const std::array<double, 2> xupper{+1.10, +0.6};

  const int n_level = 3;

  auto embedded_boundary =
      amrex::EB2::makeUnion(Rectangle({-1.0, +0.015}, {0.0, 1.0}),
                            Rectangle({-1.0, -1.0}, {0.0, -0.015}));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::Burke2012 mech{};
  fub::IdealGasMix<2> equation(mech);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  amrex::Geometry coarse_geom(amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1}},
                              &geometry.coordinates, -1,
                              geometry.periodicity.data());

  using namespace fub::amrex::cutcell;

  PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces = MakeIndexSpaces(shop, coarse_geom, n_level);

  //  auto hierarchy = fub::amrex::cutcell::ReadCheckpointFile(
  //      "LinearShock2d/Checkpoint", desc, geometry, options);

  fub::Conservative<fub::IdealGasMix<2>> cons;
  fub::FlameMasterReactor& reactor = equation.GetReactor();
  reactor.SetMoleFractions("O2:20,N2:80");
  reactor.SetTemperature(300.0);
  reactor.SetPressure(101325.0);
  fub::Complete<fub::IdealGasMix<2>> right{equation};
  equation.CompleteFromReactor(right);

  reactor.SetPressure(4 * 101325.0);
  fub::Complete<fub::IdealGasMix<2>> left{equation};
  equation.CompleteFromReactor(left);

  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04),
                              left, right);

  using State = fub::Complete<fub::IdealGasMix<2>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.05}};

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, options), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<2>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::ideal_gas::MusclHancockPrimMethod<2> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{fub::execution::simd, cutcell_method},
                          TimeIntegrator{},
                          Reconstruction{fub::execution::simd, equation}};

  fub::DimensionalSplitLevelIntegrator solver(
      fub::int_c<2>, IntegratorContext(gridding, method),
      fub::GodunovSplitting());

  std::string base_name = "LinearShock2d";

  auto output = [&](const std::shared_ptr<GriddingAlgorithm>& gridding,
                    std::ptrdiff_t cycle, fub::Duration, int = 0) {
    std::string name = fmt::format("{}/plt{:04}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    WritePlotFile(name, gridding->GetPatchHierarchy(), equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.final_time = 0.002s;
  run_options.output_interval = {0.0000125s};
  run_options.cfl = 0.8;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
