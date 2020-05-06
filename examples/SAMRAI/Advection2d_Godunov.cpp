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

#include <SAMRAI/appu/VisItDataWriter.h>
#include <SAMRAI/hier/IntVector.h>

#include <fmt/format.h>
#include <iostream>

#include <xmmintrin.h>

using namespace fub;

struct CircleData {
  using Complete = fub::Complete<Advection2d>;
  samrai::DataDescription data_description_;
  Advection2d equation_;

  void InitializeData(std::shared_ptr<SAMRAI::hier::PatchLevel> patch_level,
                      const samrai::GriddingAlgorithm& grid, int level_number,
                      Duration) const {
    SAMRAI::hier::PatchLevel& level = *patch_level;
    const samrai::PatchHierarchy& hierarchy = grid.GetPatchHierarchy();
    const SAMRAI::geom::CartesianGridGeometry& geom =
        hierarchy.GetGeometry(level_number);
    span<const int> data_ids{data_description_.data_ids};
    std::vector<SAMRAI::pdat::CellData<double>*> datas(data_ids.size());
    for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
      samrai::GetPatchData(span{datas}, *patch, data_ids);
      View<Complete> states = samrai::MakeView<Complete>(
          span{datas}, equation_, samrai::AsIndexBox<2>(patch->getBox()));
      ForEachIndex(Box<0>(states), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
        std::array<double, 2> x = samrai::GetCellCenter(geom, i, j);
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

  using namespace fub::samrai;

  const ScopeGuard guard(argc, argv);
  InitializeLogging(MPI_COMM_WORLD);

  Advection2d equation{{1.0, 0.6}};
  constexpr int Dim = Advection2d::Rank();
  SAMRAI::tbox::Dimension dim(Dim);

  const std::array<int, Dim> n_cells{128, 128};
  const CoordinateRange<Dim> coordinates{{-1.0, -1.0}, {+1.0, +1.0}};
  std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geometry =
      MakeCartesianGridGeometry(n_cells, coordinates);
  geometry->initializePeriodicShift(SAMRAI::hier::IntVector(dim, 1));

  PatchHierarchyOptions hier_opts{SAMRAI::hier::IntVector(dim, 2), 3};

  using State = Advection2d::Complete;
  samrai::GradientDetector gradient{equation, std::pair{&State::mass, 1e-3}};

  GodunovMethod godunov_method{equation};

  DataDescription desc = RegisterVariables(equation);
  IntegratorContext::AuxialiaryDataDescription aux_desc =
      IntegratorContext::RegisterVariables(desc, godunov_method);

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(desc, geometry, hier_opts), CircleData{desc, equation},
      gradient, std::vector{2, 2, 2});
  gridding->InitializeHierarchy();

  samrai::HyperbolicMethod method{
      samrai::FluxMethodAdapter(execution::seq, godunov_method),
      samrai::HyperbolicTimeIntegrator(),
      samrai::CompleteFromConsCalculation(execution::seq, equation)};

  DimensionalSplitLevelIntegrator level_integrator(
      int_c<Dim>, IntegratorContext(gridding, method, aux_desc));

  NoSubcycleSolver solver(level_integrator);

  SAMRAI::appu::VisItDataWriter writer(dim, "VisItWriter", "SAMRAI/Advection");
  writer.registerPlotQuantity("Mass", "SCALAR", desc.data_ids[0]);

  auto output = fub::MakeOutput<fub::samrai::GriddingAlgorithm>(
      {}, {Duration(0.1)}, [&](const GriddingAlgorithm& grid) {
        SAMRAI::tbox::pout << "Start output to 'SAMRAI/Advection'.\n";
        writer.writePlotData(grid.GetPatchHierarchy().GetNative(),
                             grid.GetPatchHierarchy().GetCycles(), 
                             grid.GetPatchHierarchy().GetTimePoint().count());
        SAMRAI::tbox::pout << "Finished output to 'SAMRAI/Advection'.\n";
      });

  RunOptions run_options{};
  run_options.final_time = Duration(2.0);
  run_options.cfl = 0.8;
  (*output)(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, run_options, wall_time_reference, *output);
}
