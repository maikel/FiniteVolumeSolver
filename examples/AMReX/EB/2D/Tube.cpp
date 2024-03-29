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

#include <fub/AMReX.hpp>
#include <fub/AMReX_CutCell.hpp>
#include <fub/Solver.hpp>

#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Rotation.H>

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard _(argc, argv);

  static constexpr int Rank = AMREX_SPACEDIM;

  const std::array<int, Rank> n_cells{128, 128};
  const std::array<double, Rank> xlower{-1.0, -1.0};
  const std::array<double, Rank> xupper{+1.0, +1.0};
  amrex::RealBox xbox(xlower, xupper);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = xbox;
  geometry.periodicity = std::array<int, 2>{1, 1};

  auto embedded_boundary = ::amrex::EB2::rotate(
      ::amrex::EB2::BoxIF({-0.5, -0.5}, {+0.5, +0.5}, true), 0.25 * M_PI, 2);
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::PerfectGas<2> equation;

  using namespace fub::amrex::cutcell;

  PatchHierarchyOptions options{};
  options.max_number_of_levels = 1;
  options.index_spaces = MakeIndexSpaces(shop, geometry, options);

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.0;
  cons.momentum << 0.0, 0.0;
  cons.energy = equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> right;
  fub::CompleteFromCons(equation, right, cons);

  cons.energy *= 5;
  fub::Complete<fub::PerfectGas<2>> left;
  fub::CompleteFromCons(equation, left, cons);

  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, +1.0, 0.0}, 0.2),
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

  fub::FluxMethod<fub::perfect_gas::MusclHancockPrim<2>> flux_method{equation};
//  fub::MusclHancockMethod flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(std::move(flux_method));

  HyperbolicMethod method{FluxMethod{cutcell_method},
                          TimeIntegrator{}, Reconstruction{equation}};

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method, 4, 2),
      fub::StrangSplitting{});

  fub::SubcycleFineFirstSolver solver(level_integrator);

  std::string base_name = "Tube/";
  using namespace std::literals::chrono_literals;
  fub::AnyOutput<GriddingAlgorithm> output(
      {1}, {}, [&](const GriddingAlgorithm& gridding) {
        std::string name =
            fmt::format("{}plt{:05}", base_name, gridding.GetCycles());
        ::amrex::Print() << "Start output to '" << name << "'.\n";
        fub::amrex::cutcell::WritePlotFile(name, gridding.GetPatchHierarchy(),
                                           equation);
        ::amrex::Print() << "Finished output to '" << name << "'.\n";
      });

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.final_time = 1e4s;
  run_options.cfl = 0.5 * 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}
