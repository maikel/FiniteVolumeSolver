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

#include "fub/Solver.hpp"
#include "fub/AMReX.hpp"

#include "fub/equations/CompressibleAdvection.hpp"

struct InitialData {
  using Complete = fub::CompressibleAdvection<2>::Complete;
  void InitializeData(amrex::MultiFab& mf, const amrex::Geometry& geom) {
    fub::amrex::ForEachFab(mf, [&](const amrex::MFIter& mfi) {
      fub::CompressibleAdvection<2> equation{};
      amrex::FArrayBox& fab = mf[mfi];
      const amrex::Box& box = mfi.tilebox();
      fub::View<Complete> states = fub::amrex::MakeView<Complete>(fab, equation, box);
      fub::ForEachIndex(fub::Box<0>(states), [&](int i, int j) {
        const double x = geom.CellCenter(i, 0);
        const double y = geom.CellCenter(j, 0);
        if (x*x + y*y < 0.25 * 0.25) {
          states.density(i, j) = 42.0;
          states.momentum(i, j, 0) = 0.0;
          states.momentum(i, j, 1) = 0.0;
          states.PTdensity(i, j) = 1.0;
          states.velocity(i, j, 0) = 0.0;
          states.velocity(i, j, 1) = 0.0;
          states.PTinverse(i, j) = 42.0;
        } else {
          states.density(i, j) = 24.0;
          states.momentum(i, j, 0) = 0.0;
          states.momentum(i, j, 1) = 0.0;
          states.PTdensity(i, j) = 1.0;
          states.velocity(i, j, 0) = 0.0;
          states.velocity(i, j, 1) = 0.0;
          states.PTinverse(i, j) = 24.0;
        }
      });
    });
  }
};

int main()
{
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard guard{};
  fub::amrex::DataDescription desc{};
  desc.n_state_components = 7;
  desc.n_cons_components = 4;

  fub::amrex::CartesianGridGeometry geometry{};
  geometry.cell_dimensions = std::array<int, 2>{128, 128};
  geometry.coordinates = amrex::RealBox({-1.0, -1.0}, {+1.0, +1.0});
  geometry.periodicity = std::array<int, 2>{1, 1};

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 1;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 2, 1)};

  fub::amrex::PatchHierarchy hierarchy(desc, geometry, hier_opts);

  std::shared_ptr grid = std::make_shared<fub::amrex::GriddingAlgorithm>(std::move(hierarchy), InitialData{}, fub::amrex::TagBuffer{0});
  grid->InitializeHierarchy(0.0);

  fub::CompressibleAdvection<2> equation{};

  using namespace fub;
  CompressibleAdvectionFluxMethod<2> flux_method;
  flux_method.Pv_function_ = [](std::array<double, 2>, Duration, Direction dir) -> double {
    switch (dir) {
      case Direction::X:
        return 1.0;
      default:
        return 0.0;
    }
  };

  auto tag = fub::execution::seq;
  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethod(tag, flux_method),
      fub::amrex::ForwardIntegrator(tag),
      fub::amrex::Reconstruction(tag, equation)};

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, fub::amrex::IntegratorContext(grid, method, 3, 2),
      fub::StrangSplitting());

  fub::SubcycleFineFirstSolver solver(level_integrator);

  using namespace std::literals::chrono_literals;
  std::string base_name = "CompressibleAdvection";

  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> output{};
  output.AddOutput(fub::MakeOutput<fub::amrex::GriddingAlgorithm>(
      {}, {0.1s}, fub::amrex::PlotfileOutput(equation, base_name)));

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.final_time = 3.0s;
  run_options.cfl = 0.8;
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}