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

#include <cfenv>

struct RiemannProblem {
  using Equation = fub::PerfectGas<1>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  Equation equation_;
  Complete left_{equation_};
  Complete right_{equation_};

  void InitializeData(amrex::MultiFab& data, const amrex::Geometry& geom) {
    fub::amrex::ForEachFab(
        fub::execution::openmp, data, [&](amrex::MFIter& mfi) {
          fub::View<Complete> state = fub::amrex::MakeView<Complete>(
              data[mfi], equation_, mfi.tilebox());
          fub::ForEachIndex(fub::Box<0>(state),
                            [this, &state, &geom](std::ptrdiff_t i) {
                              const double x = geom.CellCenter(int(i), 0);
                              if (x < -0.04) {
                                Store(state, right_, {i});
                              } else if (x < 0.04) {
                                Store(state, left_, {i});
                              } else {
                                Store(state, right_, {i});
                              }
                            });
        });
  }
};

int main(int argc, char** argv) {
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::amrex::ScopeGuard guard(argc, argv);
  fub::InitializeLogging(MPI_COMM_WORLD);

  constexpr int Dim = AMREX_SPACEDIM;

  const std::array<int, Dim> n_cells{AMREX_D_DECL(128, 1, 1)};
  const std::array<double, Dim> xlower{AMREX_D_DECL(-1.0, -1.0, -1.0)};
  const std::array<double, Dim> xupper{AMREX_D_DECL(+1.0, +1.0, +1.0)};

  fub::PerfectGas<1> equation{};

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  using Complete = fub::PerfectGas<1>::Complete;
  fub::amrex::GradientDetector gradient{
      equation, std::make_pair(&Complete::density, 5e-3),
      std::make_pair(&Complete::pressure, 5e-2)};

  auto from_prim = [](Complete& state, const fub::PerfectGas<1>& equation) {
    state.energy = state.pressure * equation.gamma_minus_1_inv;
    state.speed_of_sound =
        std::sqrt(equation.gamma * state.pressure / state.density);
  };

  Complete left;
  left.density = 1.0;
  left.momentum = 0.0;
  left.pressure = 8.0;
  from_prim(left, equation);

  Complete right;
  right.density = 1.0;
  right.momentum = 0.0;
  right.pressure = 1.0;
  from_prim(right, equation);

  RiemannProblem initial_data{equation, left, right};

  fub::amrex::BoundarySet boundary;
  using fub::amrex::ReflectiveBoundary;
  using fub::amrex::TransmissiveBoundary;
  auto seq = fub::execution::seq;
  boundary.conditions.push_back(
      ReflectiveBoundary{seq, equation, fub::Direction::X, 0});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 1});

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 3;
  hier_opts.refine_ratio = ::amrex::IntVect{AMREX_D_DECL(2, 1, 1)};

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts), initial_data,
      gradient, boundary);
  gridding->InitializeHierarchy(0.0);

  auto tag = fub::execution::simd;

  fub::EinfeldtSignalVelocities<fub::PerfectGas<1>> signals{};
  fub::HllMethod hll_method(equation, signals);
  fub::MusclHancockMethod flux_method{equation, hll_method};
  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethod(tag, flux_method),
      fub::amrex::ForwardIntegrator(tag),
      fub::amrex::Reconstruction(tag, equation)};

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<1>, fub::amrex::IntegratorContext(gridding, method),
      fub::GodunovSplitting{});

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  std::string base_name = "PerfectGas1d/";

  using namespace fub::amrex;
  using namespace std::literals::chrono_literals;
  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> output{};
  output.AddOutput(fub::MakeOutput<GriddingAlgorithm>(
      {}, {1.0s / 180.}, fub::amrex::PlotfileOutput(equation, base_name)));
  output.AddOutput(
      std::make_unique<
          fub::CounterOutput<fub::amrex::GriddingAlgorithm>>(
          solver.GetContext().registry_, wall_time_reference,
          std::vector<std::ptrdiff_t>{}, std::vector<fub::Duration>{0.5s}));

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.cfl = 0.8;
  run_options.final_time = 1.0s;
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}
