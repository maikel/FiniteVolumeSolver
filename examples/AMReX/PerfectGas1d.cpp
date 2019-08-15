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
          fub::ForEachIndex(fub::amrex::AsIndexBox<1>(mfi.tilebox()),
                            [this, &state, &geom](auto... is) {
                              double x[AMREX_SPACEDIM] = {};
                              geom.CellCenter(amrex::IntVect{int(is)...}, x);
                              if (x[0] < -0.04) {
                                Store(state, left_, {is...});
                              } else {
                                Store(state, right_, {is...});
                              }
                            });
        });
  }
};

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::amrex::ScopeGuard guard(argc, argv);

  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO |
                         _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW |
                         _MM_MASK_INVALID);

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
  using fub::amrex::TransmissiveBoundary;
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 0});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 1});

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 5;
  hier_opts.refine_ratio = ::amrex::IntVect{AMREX_D_DECL(2, 1, 1)};

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts), initial_data,
      gradient, boundary);
  gridding->InitializeHierarchy(0.0);

  auto tag = fub::execution::seq;

  fub::EinfeldtSignalVelocities<fub::PerfectGas<1>> signals{};
  fub::HllMethod hll_method(equation, signals);
  fub::MusclHancockMethod flux_method{equation, hll_method};
  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethod(tag, flux_method),
      fub::amrex::ForwardIntegrator(tag),
      fub::amrex::Reconstruction(tag, equation)};

  fub::DimensionalSplitLevelIntegrator solver(
      fub::int_c<1>, fub::amrex::IntegratorContext(gridding, method),
      fub::GodunovSplitting{});

  std::string base_name = "PerfectGas1d/";

  using namespace fub::amrex;
  auto output = [&](const std::shared_ptr<GriddingAlgorithm>& gridding,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}plt{:05}", base_name, cycle);
    amrex::Print() << "Start output to '" << name << "'.\n";
    WritePlotFile(name, gridding->GetPatchHierarchy(), equation);
    amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.cfl = 0.8;
  run_options.final_time = 2.0s;
  run_options.max_cycles = 100;
  run_options.output_frequency = 1;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
