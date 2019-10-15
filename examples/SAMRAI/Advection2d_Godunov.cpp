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

#include "fub/SAMRAI.hpp"
#include "fub/Solver.hpp"

#include <fmt/format.h>
#include <iostream>

#include <xmmintrin.h>

// View Complete == struct V ( PatchDataView mass );
// Complete == struct C ( double mass );
struct CircleData {
  using Complete = fub::Complete<fub::Advection2d>;
  fub::samrai::DataDescription data_description_;
  fub::Advection2d equation_;

  void InitializeData(fub::samrai::PatchHierarchy& hierarchy, int level_number, Duration) const {
    SAMRAI::hier::PatchLevel& level = *hierarchy.GetPatchLevel(level_number);
    for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
      fub::View<Complete> states = fub::samrai::MakeView<Complete>(patch, equation_, hierarchy.GetDataDescription());
      const SAMRAI::geom::CartesianPatchGeometry& geom = GetPatchGeometry(patch);
      fub::ForEachIndex(
          fub::Box<0>(states), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
            std::array<double, 2> x = GetCellCenter(patch, i, j);
            const double norm = std::sqrt(x[0] * x[0] + x[1] * x[1]);
            if (norm < 0.25) {
              states.mass(i, j) = 3.0;
            } else {
              states.mass(i, j) = 1.0;
            }
          });
    }
  }
};

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::amrex::ScopeGuard guard(argc, argv);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM >= 2);

  const std::array<int, Dim> n_cells{AMREX_D_DECL(128, 128, 1)};
  const std::array<double, Dim> xlower{AMREX_D_DECL(-1.0, -1.0, -1.0)};
  const std::array<double, Dim> xupper{AMREX_D_DECL(+1.0, +1.0, +1.0)};

  fub::Advection2d equation{{1.0, 1.0}};

  fub::samrai::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = std::array<int, Dim>{AMREX_D_DECL(1, 1, 1)};

  fub::samrai::PatchHierarchyOptions hier_opts{};
  hier_opts.max_number_of_levels = 3;

  using State = fub::Advection2d::Complete;
  fub::samrai::GradientDetector gradient{equation,
                                        std::pair{&State::mass, 1e-3}};

  fub::samrai::BoundarySet boundary;
  using fub::samrai::TransmissiveBoundary;
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 0});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 1});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::Y, 0});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::Y, 1});

  std::shared_ptr gridding = std::make_shared<fub::samrai::GriddingAlgorithm>(
      fub::samrai::PatchHierarchy(equation, geometry, hier_opts), CircleData{equation},
      gradient, boundary);
  gridding->InitializeHierarchy(0.0);

  fub::samrai::HyperbolicMethod method{
      fub::samrai::FluxMethod(fub::execution::seq, fub::GodunovMethod{equation}),
      fub::samrai::ForwardIntegrator(fub::execution::seq),
      fub::samrai::Reconstruction(fub::execution::seq, equation)};

  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      equation, fub::samrai::IntegratorContext(gridding, method)));

  std::string base_name = "Advection_Godunov/";

  auto output = [&](const fub::samrai::PatchHierarchy& hierarchy,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}{:04}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::samrai::WritePlotFile(name, hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };
  auto print_msg = [](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(solver.GetPatchHierarchy(), 0, 0s);
  fub::RunOptions run_options{};
  run_options.final_time = 2.0s;
  run_options.output_interval = 0.1s;
  run_options.cfl = 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}
