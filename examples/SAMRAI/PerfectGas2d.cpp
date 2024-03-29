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
#include "fub/equations/PerfectGas.hpp"
#include "fub/flux_method/HllMethod.hpp"
#include "fub/tagging_method/GradientDetector.hpp"

#include "fub/grid/SAMRAI/CartesianPatchHierarchy.hpp"
#include "fub/grid/SAMRAI/GriddingAlgorithm.hpp"
#include "fub/grid/SAMRAI/HyperbolicSplitIntegratorContext.hpp"
#include "fub/grid/SAMRAI/RegisterVariables.hpp"
#include "fub/grid/SAMRAI/ScopeGuard.hpp"
#include "fub/grid/SAMRAI/ViewPatch.hpp"

#include <SAMRAI/appu/VisItDataWriter.h>
#include <SAMRAI/tbox/PIO.h>

#include <fmt/format.h>

#include <iostream>

struct CircleData {
  void InitializeData(fub::View<fub::Complete<fub::PerfectGas<2>>> states,
                      const fub::CartesianCoordinates& x) const {
    fub::ForEachIndex(fub::Mapping<0>(states), [&](int i, int j) {
      const double norm = x(i, j).norm();
      fub::Conservative<fub::PerfectGas<2>> state;
      if (norm < 0.25) {
        state.energy = 8 * 101325. * equation_.gamma_minus_1_inv;
      } else {
        state.energy = 101325. * equation_.gamma_minus_1_inv;
      }
      state.density = 1.0;
      state.momentum.fill(0);

      fub::Complete<fub::PerfectGas<2>> complete;
      CompleteFromCons(equation_, complete, state);
      Store(states, complete, {i, j});
    });
  }

  fub::PerfectGas<2> equation_;
};

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  // We have to initialize the SAMRAI framework with this scope guard.
  const fub::samrai::ScopeGuard guard(argc, argv);

  fub::PerfectGas<2> equation{};

  // Create a Hierarchy of specified extents
  const SAMRAI::tbox::Dimension dim(equation.Rank());
  const SAMRAI::hier::BlockId block_id(0);
  const SAMRAI::hier::Index lower(dim, 0);
  const SAMRAI::hier::Index upper(dim, 99);
  const SAMRAI::hier::Box coarse_domain(lower, upper, block_id);
  const fub::samrai::CoordinatesRange coordinates{{-1.0, -1.0, -1.0},
                                                  {1.0, 1.0, 1.0}};

  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy =
      fub::samrai::CartesianPatchHierarchy(
          coarse_domain, coordinates,
          {.max_number_of_levels = 4, .periodic_dimensions = {1, 1}});

  // Setup our gridding algorithm which manages refining and load balancing

  using Complete = fub::PerfectGas<2>::Complete;
  fub::GradientDetector tagging{std::pair{&Complete::pressure, 1e4}};
  CircleData initial_data{};
  const fub::samrai::DataDescription reg =
      fub::samrai::RegisterVariables(equation);
  fub::samrai::GriddingAlgorithm gridding(
      hierarchy, reg, fub::samrai::AdaptInitialData(initial_data, equation),
      fub::samrai::AdaptTagging(tagging, equation), {2, 2, 2});

  // Setup our Numeric Methods
  // These will register data ids with SAMRAIs variable database

  fub::HyperbolicSplitPatchIntegrator patch_integrator{equation};
  fub::HllMethod flux_method{equation,
                             fub::EinfeldtSignalVelocities<fub::PerfectGas<2>>};
  fub::HyperbolicSplitSystemSolver system_solver(
      fub::HyperbolicSplitLevelIntegrator(
          fub::samrai::HyperbolicSplitIntegratorContext(
              gridding, reg, flux_method.GetStencilWidth()),
          patch_integrator, flux_method));

  // Initialize hierarchy and reset the solvers configuration to update
  // communication schedules.
  gridding.InitializeHierarchy();
  hierarchy->recursivePrint(SAMRAI::tbox::pout, "", 2);
  system_solver.ResetHierarchyConfiguration();

  SAMRAI::appu::VisItDataWriter writer(dim, "VisItWriter",
                                       "SAMRAI/PerfectGas2d");
  writer.registerPlotQuantity("Density", "SCALAR", reg.data_ids[0]);
  writer.registerPlotQuantity("Momentum", "VECTOR", reg.data_ids[1]);
  writer.registerPlotQuantity("Energy", "SCALAR", reg.data_ids[2]);
  writer.registerPlotQuantity("Pressure", "SCALAR", reg.data_ids[3]);
  writer.registerPlotQuantity("SpeedOfSound", "SCALAR", reg.data_ids[4]);

  auto output =
      [&writer](const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                std::ptrdiff_t cycle,
                std::chrono::duration<double> time_point) {
        SAMRAI::tbox::pout << "Start VisIt Output.\n";
        writer.writePlotData(hierarchy, cycle, time_point.count());
        SAMRAI::tbox::pout << "Finished VisIt Output.\n";
      };

  using namespace std::literals::chrono_literals;
  output(hierarchy, 0, 0.0s);

  auto print_msg = [](const std::string& msg) { SAMRAI::tbox::pout << msg; };

  fub::RunOptions options{};
  options.final_time = 0.002s;
  options.output_interval = 0.0001s;
  fub::RunSimulation(system_solver, options, wall_time_reference, output,
                     print_msg);
}