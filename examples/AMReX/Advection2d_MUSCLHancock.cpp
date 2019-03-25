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

#include "fub/equations/Advection.hpp"
#include "fub/CartesianCoordinates.hpp"

#include "fub/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/HyperbolicSplitSystemSolver.hpp"
#include "fub/ext/Eigen.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"
#include "fub/grid/AMReX/GriddingAlgorithm.hpp"
#include "fub/grid/AMReX/FluxMethod.hpp"
#include "fub/grid/AMReX/HyperbolicSplitIntegratorContext.hpp"
#include "fub/grid/AMReX/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/grid/AMReX/Reconstruction.hpp"
#include "fub/grid/AMReX/ScopeGuard.hpp"
#include "fub/tagging/GradientDetector.hpp"
#include "fub/tagging/TagBuffer.hpp"
#include "fub/split_method/StrangSplitting.hpp"

#include "fub/RunSimulation.hpp"

#include <fmt/format.h>
#include <iostream>

struct CircleData {
  std::shared_ptr<fub::amrex::PatchHierarchy> hierarchy;
  void InitializeData(const fub::View<fub::Complete<fub::Advection2d>>& states,
                      fub::amrex::PatchHandle patch) const {
    const ::amrex::Geometry& geom = hierarchy->GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->tilebox();
    fub::CartesianCoordinates x =
        fub::amrex::GetCartesianCoordinates(geom, box);
    fub::ForEachIndex(fub::Box<0>(states),
                      [&](std::ptrdiff_t i, std::ptrdiff_t j) {
                        const double norm = x(i, j).norm();
                        if (norm < 0.25) {
                          states.mass(i, j) = 3;
                        } else {
                          states.mass(i, j) = 1;
                        }
                      });
  }
};

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::amrex::ScopeGuard guard(argc, argv);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM == 2);

  const std::array<int, Dim> n_cells{256, 256};

  const std::array<double, Dim> xlower{-1.0, -1.0};
  const std::array<double, Dim> xupper{+1.0, +1.0};

  fub::Advection2d equation{{1.0, 1.0}};
  
  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = std::array<int, 2>{1, 1};
  
  fub::amrex::PatchHierarchyOptions options;
  options.max_number_of_levels = 2;

  auto hierarchy =
      std::make_shared<fub::amrex::PatchHierarchy>(desc, geometry, options);

  using State = fub::Advection2d::Complete;
  fub::GradientDetector gradient{std::pair{&State::mass, 1e-3}};
  fub::TagBuffer buffer{2};

  CircleData initial_data{hierarchy};
  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      hierarchy, fub::amrex::AdaptInitialData(initial_data, equation),
      fub::amrex::AdaptTagging(equation, hierarchy, gradient, buffer));
  gridding->InitializeHierarchy(0.0);
  fub::amrex::WritePlotFile("AMReX/Advection_M2_0000", *hierarchy, equation);

  fub::HyperbolicSplitPatchIntegrator patch_integrator{equation};
  fub::MusclHancockMethod flux_method{equation};

  const int gcw = flux_method.GetStencilWidth();
  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      fub::amrex::HyperbolicSplitIntegratorContext(std::move(gridding), gcw),
      fub::amrex::HyperbolicSplitPatchIntegrator(patch_integrator),
      fub::amrex::FluxMethod(flux_method),
      fub::amrex::Reconstruction(equation)));

  std::string base_name = "AMReX/Advection_M2_";

  auto output = [&](auto& hierarchy, int cycle, fub::Duration) {
    std::string name = fmt::format("{}{:04}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::WritePlotFile(name, *hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };
  auto print_msg = [](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(hierarchy, 0, 0s);
  fub::RunOptions run_options{};
  run_options.final_time = 2.0s;
  run_options.output_interval = 2.0s;
  run_options.output_frequency = 1;
  run_options.cfl = 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}