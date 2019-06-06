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

#include <fmt/format.h>
#include <iostream>

struct ShockData {
  using Equation = fub::PerfectGas<2>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  void InitializeData(amrex::MultiFab& data, const amrex::Geometry& geom) {
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(data, true); mfi.isValid(); ++mfi) {
      fub::View<Complete> state =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::amrex::ForEachIndex(
          mfi.tilebox(),
          [this, &state, &geom](std::ptrdiff_t i, std::ptrdiff_t j) {
            double x[2];
            geom.CellCenter({int(i), int(j)}, x);
            if (x[0] < -0.04) {
              Store(state, left_, {i, j});
            } else {
              Store(state, right_, {i, j});
            }
          });
    }
  }

  Equation equation_;
  Complete left_{};
  Complete right_{};
};

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::amrex::ScopeGuard guard(argc, argv);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM >= 2);

  fub::PerfectGas<2> equation{};

  fub::amrex::CartesianGridGeometry geometry{};
  geometry.cell_dimensions = std::array<int, Dim>{AMREX_D_DECL(128, 128, 1)};
  geometry.coordinates = amrex::RealBox({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                                        {AMREX_D_DECL(+1.0, +1.0, +1.0)});
  geometry.periodicity[1] = 1;

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 2;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 2, 1)};

  using Complete = fub::PerfectGas<2>::Complete;
  using CompleteArray = fub::PerfectGas<2>::CompleteArray;
  fub::amrex::GradientDetector gradient(
      equation, std::make_pair(&CompleteArray::density, 0.001),
      std::make_pair(&CompleteArray::pressure, 0.01));

  auto from_prim = [](Complete& state, const fub::PerfectGas<2>& equation) {
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

  fub::amrex::BoundarySet boundary;
  using fub::amrex::TransmissiveBoundary;
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 0});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 1});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::Y, 0});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::Y, 1});

  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts),
      ShockData{equation, left, right}, gradient, boundary);
  gridding->InitializeHierarchy(0.0);

  fub::HllMethod flux_method{
      equation, fub::EinfeldtSignalVelocities<fub::PerfectGas<2>>{}};
  fub::amrex::NumericalMethod method(std::move(flux_method));

  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      equation, fub::amrex::HyperbolicSplitIntegratorContext(
                    std::move(gridding), std::move(method))));

  std::string base_name = "PerfectGas2d/";

  auto output = [&](const fub::amrex::PatchHierarchy& hierarchy,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::WritePlotFile(name, hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  using namespace std::literals::chrono_literals;
  output(solver.GetPatchHierarchy(), 0, 0.0s);
  fub::RunOptions run_options{};
  run_options.final_time = 2s;
  run_options.output_interval = 0.1s;
  run_options.cfl = 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
