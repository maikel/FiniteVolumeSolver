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
#include "fub/Solver.hpp"

#include <fmt/format.h>
#include <iostream>

struct CircleData {
  CircleData(const fub::ShallowWater& eq) : equation_{eq} {
    inner_.height = 1.4;
    inner_.momentum = Eigen::Array<double, 2, 1>::Zero();
    outer_.height = 1.0;
    outer_.momentum = inner_.momentum;
  }

  void InitializeData(const fub::View<fub::Complete<fub::ShallowWater>>& states,
                      const fub::amrex::PatchHierarchy& hierarchy,
                      fub::amrex::PatchHandle patch) const {
    const ::amrex::Geometry& geom = hierarchy.GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->tilebox();
    fub::CartesianCoordinates x =
        fub::amrex::GetCartesianCoordinates(geom, box);
    fub::ForEachIndex(fub::Box<0>(states), [&](auto... is) {
      const double norm = x(is...).norm();
      if (norm < 0.25) {
        Store(states, inner_, {is...});
      } else {
        Store(states, outer_, {is...});
      }
    });
  }

  fub::ShallowWater equation_;
  fub::Complete<fub::ShallowWater> inner_{};
  fub::Complete<fub::ShallowWater> outer_{};
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

  fub::ShallowWater equation{};

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = ::amrex::RealBox(xlower, xupper);
  geometry.periodicity = std::array<int, Dim>{AMREX_D_DECL(1, 1, 1)};

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 3;

  using State = fub::ShallowWater::Complete;
  fub::GradientDetector gradient{equation, std::pair(&State::height, 1e-2)};

  fub::amrex::BoundarySet boundary;
  using fub::amrex::TransmissiveBoundary;
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 0});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 1});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::Y, 0});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::Y, 1});

  fub::amrex::GriddingAlgorithm gridding(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts),
      CricleData{equation, left, right}, gradient, boundary);
  gridding.InitializeHierarchy(0.0);

  fub::HyperbolicSplitPatchIntegrator patch_integrator{equation};
  fub::HllMethod hll_method{equation, fub::ShallowWaterSignalVelocities{}};
  fub::MusclHancockMethod muscl_method{equation, hll_method};
  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethod(execution::seq, muscl_method),
      fub::amrex::ForwardIntegrator(execution::seq),
      fub::amrex::Reconstruction(execution::seq, equation)};
  fub::HyperbolicSplitSystemSolver solver(
      equation, fub::amrex::IntegratorContext(gridding, method));

  std::string base_name = "ShallowWater2d/";

  auto output = [&](const fub::amrex::PatchHierarchy& hierarchy,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::WritePlotFile(name, hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  using namespace std::literals::chrono_literals;
  output(solver.GetPatchHierarchy(), 0, 0.0s);
  fub::RunOptions run_options{};
  run_options.final_time = 1.0s;
  run_options.output_interval = 0.01s;
  run_options.cfl = 0.8;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
