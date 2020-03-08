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
  using Complete = fub::Complete<fub::ShallowWater>;

  CircleData(const fub::ShallowWater& eq) : equation_{eq} {
    inner_.height = 1.4;
    inner_.momentum = Eigen::Array<double, 2, 1>::Zero();
    outer_.height = 1.0;
    outer_.momentum = inner_.momentum;
  }

  void InitializeData(amrex::MultiFab& data,
                      const amrex::Geometry& geom) const {
    fub::amrex::ForEachFab(data, [&](const amrex::MFIter& mfi) {
      const ::amrex::Box& box = mfi.tilebox();
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, box);
      fub::ForEachIndex(
          fub::Box<0>(states), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
            double x[AMREX_SPACEDIM] = {};
            geom.CellCenter(amrex::IntVect{AMREX_D_DECL(int(i), int(j), 0)}, x);
            const double norm = std::sqrt(x[0] * x[0] + x[1] * x[1]);
            if (norm < 0.25) {
              Store(states, inner_, {i, j});
            } else {
              Store(states, outer_, {i, j});
            }
          });
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
  fub::InitializeLogging(MPI_COMM_WORLD);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM >= 2);

  const std::array<int, Dim> n_cells{AMREX_D_DECL(128, 128, 1)};
  const std::array<double, Dim> xlower{AMREX_D_DECL(-1.0, -1.0, -1.0)};
  const std::array<double, Dim> xupper{AMREX_D_DECL(+1.0, +1.0, +1.0)};

  fub::ShallowWater equation{};

  fub::amrex::CartesianGridGeometry geometry{};
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = ::amrex::RealBox(xlower, xupper);
  geometry.periodicity = std::array<int, Dim>{AMREX_D_DECL(1, 1, 1)};

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 4;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 2, 1)};

  using State = fub::ShallowWater::Complete;
  fub::amrex::GradientDetector gradient{equation,
                                        std::pair(&State::height, 1e-2)};

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts),
      CircleData{equation},
      fub::amrex::TagAllOf(gradient, fub::amrex::TagBuffer(2)));
  gridding->InitializeHierarchy(0.0);

  fub::HllMethod hll_method{equation, fub::ShallowWaterSignalVelocities{}};
  fub::MusclHancockMethod muscl_method{equation, hll_method};
  fub::amrex::HyperbolicMethod method{fub::amrex::FluxMethod(muscl_method),
                                      fub::amrex::EulerForwardTimeIntegrator(),
                                      fub::amrex::NoReconstruction{}};

  const int scratch_ghost_cell_width = 8;
  const int flux_ghost_cell_width = 6;

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>,
      fub::amrex::IntegratorContext(gridding, method, scratch_ghost_cell_width,
                                    flux_ghost_cell_width),
      fub::StrangSplitting());

  // fub::NoSubcycleSolver solver(std::move(level_integrator));
  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  std::string base_name = "ShallowWater2d/";

  using namespace fub::amrex;
  using namespace std::literals::chrono_literals;
  fub::RunOptions run_options{};
  run_options.final_time = 1.0s;
  run_options.cfl = 0.8;

  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> output{};

  // Add a standard AMReX Plotfile output readable in VisIt
  output.AddOutput(std::make_unique<fub::amrex::PlotfileOutput<fub::ShallowWater>>(std::vector<std::ptrdiff_t>{},
          std::vector<fub::Duration>{0.1s}, equation, base_name));

  // Add an output to print the timer database
  using std::chrono::microseconds;
  output.AddOutput(
      std::make_unique<
          fub::CounterOutput<fub::amrex::GriddingAlgorithm, microseconds>>(
          wall_time_reference, std::vector<std::ptrdiff_t>{},
          std::vector<fub::Duration>{0.05s}));

  output(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}
