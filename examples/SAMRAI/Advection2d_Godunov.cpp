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

#include <SAMRAI/hier/IntVector.h>

#include <fmt/format.h>
#include <iostream>

#include <xmmintrin.h>

// View Complete == struct V ( PatchDataView mass );
// Complete == struct C ( double mass );
struct CircleData {
  using Complete = fub::Complete<fub::Advection2d>;
  fub::samrai::DataDescription data_description_;
  fub::Advection2d equation_;

  void InitializeData(fub::samrai::PatchHierarchy& hierarchy, int level_number,
                      fub::Duration) const {
    SAMRAI::hier::PatchLevel& level = *hierarchy.GetPatchLevel(level_number);
    const SAMRAI::geom::CartesianGridGeometry& geom =
        hierarchy.GetGeometry(level_number);
    for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
      fub::View<Complete> states = fub::samrai::MakeView<Complete>(
          patch, equation_, hierarchy.GetDataIds(), patch.getBox());
      fub::ForEachIndex(
          fub::Box<0>(states), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
            std::array<double, 2> x = GetCellCenter(geom, patch, i, j);
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
  fub::InitializeLogging(MPI_COMM_WORLD);

  fub::Advection2d equation{{1.0, 0.6}};
  constexpr int Dim = fub::Advection2d::Rank();
  SAMRAI::tbox::Dimension dim(Dim);

  const std::array<int, Dim> n_cells{128, 128};
  const CoordinateRange<Dim> coordinates{{-1.0, -1.0}, {+1.0, +1.0}};
  std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geometry =
      MakeCartesianGridGeometry(n_cells, coordinates);
  geometry->initializePeriodicShift(SAMRAI::hier::IntVector(dim, 1));

  PatchHierarchyOptions hier_opts{.refine_ratio =
                                      SAMRAI::hier::IntVector(dim, 2),
                                  .max_number_of_levels = 3};

  using State = fub::Advection2d::Complete;
  GradientDetector gradient{equation, std::pair{&State::mass, 1e-3}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, hier_opts), CircleData{equation},
      gradient, std::vector{2, 2, 2});
  gridding->InitializeHierarchy(0.0);

  HyperbolicMethod method{
      FluxMethod(fub::execution::seq, fub::GodunovMethod{equation}),
      TimeIntegrator(), Reconstruction(fub::execution::seq, equation)};

  fub::DimensionalSplitLevelIntegrator solver(
      fub::int_c<Dim>, IntegratorContext(gridding, method));

  SAMRAI::appu::VisItDataWriter writer(dim, "VisItWriter", "SAMRAI/Advection");
  writer.registerPlotQuantity("mass", "SCALAR", desc.data_ids[0]);

  auto output = [&](const GriddingAlgorithm& grid) {
    SAMRAI::tbox::pout << "Start output to 'SAMRAI/Advection'.\n";
    writer.writePlotData(grid.GetHierarchy().GetNative(), grid.GetCycles(),
                         grid.GetTimePoint().count());
    SAMRAI::tbox::pout << "Finished output to 'SAMRAI/Advection'.\n";
  };

  fub::RunOptions run_options{};
  run_options.final_time = 2.0s;
  run_options.output_interval = 0.1s;
  run_options.cfl = 0.9;
  output(*solver.GetGriddingAlgorithm());
  fub::InvokeFunction<GriddingAlgorithm> out{{1}, {}, output};
  fub::RunSimulation(solver, run_options, wall_time_reference, out);
}
