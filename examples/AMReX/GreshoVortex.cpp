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

#include <xmmintrin.h>

struct GreshoVortex {
  using Equation = fub::PerfectGas<2>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  void InitializeData(amrex::MultiFab& data, const amrex::Geometry& geom) {
    fub::amrex::ForEachFab(data, [&](const amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());

      auto from_prim = [](Complete& state, const Equation& equation) {
        state.energy = state.pressure * equation.gamma_minus_1_inv + 0.5 * state.momentum.matrix().squaredNorm();
        state.speed_of_sound =
            std::sqrt(equation.gamma * state.pressure / state.density);
      };

      ForEachIndex(fub::Box<0>(states),
                   [&](std::ptrdiff_t i, std::ptrdiff_t j) {
                     double xy[AMREX_SPACEDIM];
                     geom.CellCenter({AMREX_D_DECL(int(i), int(j), 0)}, xy);
                     const double x = xy[0];
                     const double y = xy[1];

                     const double r = std::sqrt(x * x + y * y);
                     const double phi = std::atan2(y, x);

                     double pr = 0 * r;
                     double uth = 0 * r;

                     if (r < 0.2) {
                       uth = 5. * r;
                       pr = 5. + 12.5 * r * r;
                     } else if (r < 0.4) {
                       uth = 2. - 5. * r;
                       pr = 9. - 4. * std::log(0.2) + 12.5 * r * r - 20. * r +
                            4. * std::log(r);
                     } else {
                       uth = 0.;
                       pr = 3. + 4. * std::log(2.);
                     }

                     const double u = -std::sin(phi) * uth;
                     const double v = std::cos(phi) * uth;

                     Complete state;

                     state.density = 1.;
                     state.momentum[0] = u;
                     state.momentum[1] = v;
                     state.pressure = pr;

                     from_prim(state, equation_);

                     Store(states, state, {i, j});
                   });
    });
  }

  Equation equation_;
};

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO |
                         _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW |
                         _MM_MASK_INVALID);

  const fub::amrex::ScopeGuard guard(argc, argv);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM >= 2);

  fub::PerfectGas<2> equation{};

  fub::amrex::CartesianGridGeometry geometry{};
  geometry.cell_dimensions = std::array<int, Dim>{AMREX_D_DECL(200, 200, 1)};
  geometry.coordinates = amrex::RealBox({AMREX_D_DECL(-0.5, -0.5, -0.5)},
                                        {AMREX_D_DECL(+0.5, +0.5, +0.5)});
  geometry.periodicity = std::array<int, Dim>{AMREX_D_DECL(1, 1, 1)};

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 1;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 2, 1)};

  using Complete = fub::PerfectGas<2>::Complete;
  auto velocity_norm = [](const Complete& state) {
    return (state.momentum / state.density).matrix().norm();
  };

  fub::amrex::GradientDetector gradient(
      equation, std::make_pair(&Complete::density, 1.0e-5),
      std::make_pair(&Complete::pressure, 1.0e-3),
      std::make_pair(velocity_norm, 1.0e-4));

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts),
      GreshoVortex{equation},
      fub::amrex::TagAllOf(gradient, fub::amrex::TagBuffer(4)));
  gridding->InitializeHierarchy(0.0);

  auto tag = fub::execution::openmp_simd;

  //  fub::EinfeldtSignalVelocities<fub::PerfectGas<2>> signals{};
  //  fub::HllMethod hll_method(equation, signals);
  //  fub::MusclHancockMethod muscl_method(equation, hll_method);
  //  fub::GodunovMethod godunov_method(equation, signals);
  fub::MusclHancockMethod muscl_method(equation);
  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethod(fub::execution::openmp, muscl_method),
      fub::amrex::ForwardIntegrator(tag),
      fub::amrex::Reconstruction(tag, equation)};

  fub::DimensionalSplitLevelIntegrator solver(
      fub::int_c<2>, fub::amrex::IntegratorContext(gridding, method),
      fub::StrangSplitting());

  std::string base_name = "GreshoVortex/";

  using namespace fub::amrex;
  auto output = [&](const std::shared_ptr<GriddingAlgorithm>& gridding,
                    std::ptrdiff_t cycle, fub::Duration, int t = 0) {
    if (t <= 0) {
      std::string name = fmt::format("{}plt{:05}", base_name, cycle);
      amrex::Print() << "Start output to '" << name << "'.\n";
      WritePlotFile(name, gridding->GetPatchHierarchy(), equation);
      amrex::Print() << "Finished output to '" << name << "'.\n";
    }
    double rho_max = 0.0;
    double rho_min = std::numeric_limits<double>::infinity();
    for (int level = 0; level < gridding->GetPatchHierarchy().GetNumberOfLevels(); ++level) {
      const ::amrex::MultiFab& mf = gridding->GetPatchHierarchy().GetPatchLevel(level).data;
      rho_max = std::max(rho_max, mf.max(0));
      rho_min = std::min(rho_min, mf.min(0));
    }
    const double rho_err = std::max(std::abs(rho_max - 1.0), std::abs(rho_min - 1.0));
    ::amrex::Print() << fmt::format("Density Max Error: {:.6E}\n", rho_err);
  };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.final_time = 3.0s;
  run_options.output_interval = {fub::Duration(1.0 / 30.0)};
  run_options.output_frequency = {0, 1};
  run_options.cfl = 0.8;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
