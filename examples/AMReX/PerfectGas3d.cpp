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

#include "fub/CartesianCoordinates.hpp"
#include "fub/equations/PerfectGas.hpp"

#include "fub/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/HyperbolicSplitSystemSolver.hpp"
#include "fub/ext/Eigen.hpp"
#include "fub/flux_method/HllMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"
#include "fub/grid/AMReX/FluxMethod.hpp"
#include "fub/grid/AMReX/GriddingAlgorithm.hpp"
#include "fub/grid/AMReX/HyperbolicSplitIntegratorContext.hpp"
#include "fub/grid/AMReX/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/grid/AMReX/Reconstruction.hpp"
#include "fub/grid/AMReX/ScopeGuard.hpp"
#include "fub/split_method/StrangSplitting.hpp"
#include "fub/tagging/GradientDetector.hpp"
#include "fub/tagging/TagBuffer.hpp"

#include "fub/RunSimulation.hpp"

#include <fmt/format.h>
#include <iostream>

struct CircleData {
  std::shared_ptr<fub::amrex::PatchHierarchy> hierarchy;
  fub::PerfectGas<3> equation_;
  void
  InitializeData(const fub::View<fub::Complete<fub::PerfectGas<3>>>& states,
                 fub::amrex::PatchHandle patch) const {
    const ::amrex::Geometry& geom = hierarchy->GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->tilebox();
    fub::CartesianCoordinates x =
        fub::amrex::GetCartesianCoordinates(geom, box);
    fub::ForEachIndex(fub::Box<0>(states), [&](auto... is) {
      const double norm = x(is...).norm();
      fub::Conservative<fub::PerfectGas<3>> state;
      if (norm < 0.25) {
        state.energy = 8 * 101325. * equation_.gamma_minus_1_inv;
      } else {
        state.energy = 101325. * equation_.gamma_minus_1_inv;
      }
      state.density = 1.0;
      state.momentum.fill(0);
      fub::Complete<fub::PerfectGas<3>> complete;
      CompleteFromCons(equation_, complete, state);
      Store(states, complete, {is...});
    });
  }
};

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::amrex::ScopeGuard guard(argc, argv);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM == 3);

  const std::array<int, Dim> n_cells{64, 64, 64};

  const std::array<double, Dim> xlower{-1.0, -1.0, -1.0};
  const std::array<double, Dim> xupper{+1.0, +1.0, +1.0};

  fub::PerfectGas<3> equation{};

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = std::array<int, Dim>{1, 1, 1};

  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 2;

  auto hierarchy =
      std::make_shared<fub::amrex::PatchHierarchy>(desc, geometry, hier_opts);

  using Complete = fub::PerfectGas<3>::Complete;
  fub::GradientDetector gradient{std::make_pair(&Complete::density, 5e-3),
                                 std::make_pair(&Complete::pressure, 5e-2)};
  fub::TagBuffer buffer{2};

  CircleData initial_data{hierarchy, equation};
  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      hierarchy, fub::amrex::AdaptInitialData(initial_data, equation),
      fub::amrex::AdaptTagging(equation, hierarchy, gradient, buffer));
  gridding->InitializeHierarchy(0.0);

  fub::HyperbolicSplitPatchIntegrator patch_integrator{equation};
  fub::HllMethod base_method{equation,
                             fub::EinfeldtSignalVelocities<fub::PerfectGas<3>>};
  fub::MusclHancockMethod flux_method{equation, base_method};

  const int gcw = flux_method.GetStencilWidth();
  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      fub::amrex::HyperbolicSplitIntegratorContext(std::move(gridding), gcw),
      fub::amrex::HyperbolicSplitPatchIntegrator(patch_integrator),
      fub::amrex::FluxMethod(flux_method),
      fub::amrex::Reconstruction(equation)));

  std::string base_name = "AMReX/PerfectGas3d_";

  auto output = [&](auto& hierarchy, int cycle, fub::Duration) {
    std::string name = fmt::format("{}{:04}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::WritePlotFile(name, *hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  auto print_msg = [](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(hierarchy, 0, 0.0s);
  fub::RunOptions run_options{};
  run_options.final_time = 0.004s;
  run_options.output_interval = 0.00005s;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}