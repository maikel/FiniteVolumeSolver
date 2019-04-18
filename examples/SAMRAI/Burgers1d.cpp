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
#include "fub/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/HyperbolicSplitSystemSolver.hpp"
#include "fub/RunSimulation.hpp"
#include "fub/equations/Burgers.hpp"
#include "fub/flux_method/GodunovMethod.hpp"
#include "fub/tagging/GradientDetector.hpp"

#include "fub/grid/SAMRAI/CartesianPatchHierarchy.hpp"
#include "fub/grid/SAMRAI/GriddingAlgorithm.hpp"
#include "fub/grid/SAMRAI/HyperbolicSplitIntegratorContext.hpp"
#include "fub/grid/SAMRAI/RegisterVariables.hpp"
#include "fub/grid/SAMRAI/ScopeGuard.hpp"
#include "fub/grid/SAMRAI/ViewPatch.hpp"

#include <iostream>

struct WaveData {
  static double wave_package(double x) {
    return std::exp(-15 * x * x) * std::sin(4.0 * M_PI * x);
  }

  void InitializeData(fub::View<fub::Complete<fub::Burgers1d>> states,
                      const fub::CartesianCoordinates& coords) const {
    fub::Complete<fub::Burgers1d> q;
    fub::ForEachIndex(fub::Mapping<0>(states), [&](auto... is) {
      q.u = wave_package(coords(is...).norm());
      Store(states, q, {is...});
    });
  }
};

void print(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
           const fub::samrai::DataDescription& reg) {
  for (int levelnum = 0; levelnum < hierarchy->getNumberOfLevels();
       ++levelnum) {
    const SAMRAI::hier::PatchLevel& level = *hierarchy->getPatchLevel(levelnum);
    for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
      SAMRAI::pdat::CellData<double>& pdata =
          *static_cast<SAMRAI::pdat::CellData<double>*>(
              patch->getPatchData(reg.data_ids[0]).get());
      SAMRAI::hier::Box box = pdata.getBox();
      for (const SAMRAI::hier::Index& index : box) {
        SAMRAI::pdat::CellIndex cell(index);
        std::cout << index[0] << '\t' << pdata(cell) << '\n';
      }
    }
  }
  std::cout << "\n\n";
}

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::samrai::ScopeGuard guard(argc, argv);
  fub::Burgers1d equation;

  // Create a Hierarchy of specified extents
  const SAMRAI::tbox::Dimension dim(equation.Rank());
  const SAMRAI::hier::BlockId block_id(0);
  const SAMRAI::hier::Index lower(dim, 0);
  const SAMRAI::hier::Index upper(dim, 255);
  const SAMRAI::hier::Box coarse_domain(lower, upper, block_id);
  const fub::samrai::CoordinatesRange coordinates{{-1.0, -1.0}, {1.0, 1.0}};

  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy =
      fub::samrai::CartesianPatchHierarchy(
          coarse_domain, coordinates,
          {.max_number_of_levels = 4, .periodic_dimensions = {1, 1}});

  // Setup our gridding algorithm which manages refining and load balancing

  using Complete = fub::Burgers1d::Complete;
  fub::GradientDetector tagging{std::pair{&Complete::u, 1e-2}};
  WaveData initial_data{};
  const fub::samrai::DataDescription reg =
      fub::samrai::RegisterVariables(equation);
  fub::samrai::GriddingAlgorithm gridding(
      hierarchy, reg, fub::samrai::AdaptInitialData(initial_data, equation),
      fub::samrai::AdaptTagging(tagging, equation), {2, 2, 2});

  // Setup our Numeric Methods
  // These will register data ids with SAMRAIs variable database

  fub::HyperbolicSplitPatchIntegrator patch_integrator{equation};
  fub::GodunovMethod flux_method{equation};
  fub::HyperbolicSplitSystemSolver system_solver(
      fub::HyperbolicSplitLevelIntegrator(
          fub::samrai::HyperbolicSplitIntegratorContext(
              gridding, reg, flux_method.GetStencilWidth()),
          patch_integrator, flux_method));

  // Initialize hierarchy and reset the solvers configuration to update
  // communication schedules.
  gridding.InitializeHierarchy();
  system_solver.ResetHierarchyConfiguration();

  auto output =
      [&reg](const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
             std::ptrdiff_t /* cycle */,
             std::chrono::duration<double> /* time_point */) {
        // SAMRAI::tbox::pout << "Start Gnuplot Output.\n";
        // print(hierarchy, reg);
        // SAMRAI::tbox::pout << "Finished Gnuplot Output.\n";
      };

  auto print_msg = [](const std::string& msg) { SAMRAI::tbox::pout << msg; };

  using namespace std::literals::chrono_literals;
  fub::RunOptions options{};
  options.final_time = 2.0s;
  options.output_interval = 0.1s;
  fub::RunSimulation(system_solver, options, wall_time_reference, output,
                     print_msg);
}