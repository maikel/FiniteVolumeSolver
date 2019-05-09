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

#include "fub/equations/ShallowWater.hpp"

#include "fub/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/HyperbolicSplitSystemSolver.hpp"

#include "fub/ext/Eigen.hpp"

#include "fub/CartesianCoordinates.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"
#include "fub/flux_method/HllMethod.hpp"
#include "fub/AMReX/FluxMethod.hpp"
#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/HyperbolicSplitIntegratorContext.hpp"
#include "fub/AMReX/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/AMReX/Reconstruction.hpp"
#include "fub/AMReX/ScopeGuard.hpp"
#include "fub/tagging/GradientDetector.hpp"
#include "fub/tagging/TagBuffer.hpp"

#include "fub/RunSimulation.hpp"

#include <fmt/format.h>
#include <iostream>

struct CircleData {
  CircleData(const fub::ShallowWater& eq) : equation_{eq} {
    inner_.heigth = 1.4;
    inner_.momentum = Eigen::Array<double, 2, 1>::Zero();
    outer_.heigth = 1.0;
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

  const std::array<int, Dim> n_cells{AMREX_D_DECL(10 * 8, 10 * 8, 1)};
  const std::array<double, Dim> xlower{AMREX_D_DECL(-1.0, -1.0, -1.0)};
  const std::array<double, Dim> xupper{AMREX_D_DECL(+1.0, +1.0, +1.0)};

  fub::ShallowWater equation{};
  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = ::amrex::RealBox(xlower, xupper);
  geometry.periodicity = std::array<int, Dim>{AMREX_D_DECL(1, 1, 1)};

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 3;

  using State = fub::ShallowWater::Complete;
  fub::GradientDetector gradient{equation, std::pair(&State::heigth, 1e-2)};
  CircleData initial_data(equation);

  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(desc, geometry, hier_opts),
      fub::amrex::AdaptInitialData(initial_data, equation),
      fub::amrex::AdaptTagging(equation, gradient, fub::TagBuffer(4)));
  gridding->InitializeHierarchy(0.0);

  fub::HyperbolicSplitPatchIntegrator patch_integrator{equation};
  fub::HllMethod hll_method{equation, fub::ShallowWaterSignalVelocities{}};
  fub::MusclHancockMethod flux_method{equation, hll_method};

  const int gcw = flux_method.GetStencilWidth();
  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      fub::amrex::HyperbolicSplitIntegratorContext(std::move(gridding), gcw),
      fub::amrex::HyperbolicSplitPatchIntegrator(patch_integrator),
      fub::amrex::FluxMethod(flux_method),
      fub::amrex::Reconstruction(equation)));

  std::string base_name = "ShallowWater2d/";

  auto output = [&](const fub::amrex::PatchHierarchy& hierarchy,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::WritePlotFile(name, hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  auto print_msg = [](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(solver.GetPatchHierarchy(), 0, 0.0s);
  fub::RunOptions run_options{};
  run_options.final_time = 2.0s;
  run_options.output_interval = 0.01s;
  run_options.cfl = 0.5;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}
