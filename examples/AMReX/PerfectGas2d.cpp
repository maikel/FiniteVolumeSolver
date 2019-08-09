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
    fub::amrex::ForEachFab(
        fub::execution::openmp, data, [&](const amrex::MFIter& mfi) {
          fub::View<Complete> states = fub::amrex::MakeView<Complete>(
              data[mfi], equation_, mfi.tilebox());

          auto from_prim = [](Complete& state, const Equation& equation) {
            state.energy = state.pressure * equation.gamma_minus_1_inv;
            state.speed_of_sound =
            std::sqrt(equation.gamma * state.pressure / state.density);
          };

          ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
            double xy[AMREX_SPACEDIM];
            geom.CellCenter({AMREX_D_DECL(int(i), int(j), 0)}, xy);
            const double x = xy[0];
            const double y = xy[1];

            const double r   = std::sqrt(x*x + y*y);
            const double phi = std::atan2(y,x);

            double pr  = 0*r;
            double uth = 0*r;

            if (r < 0.2) {
              uth = 5. * r;
              pr  = 5. + 12.5*r*r;
            }
            else if (r < 0.4){
              uth = 2. - 5. * r;
              pr  = 9.-4.*std::log(0.2) + 12.5*r*r - 20.*r + 4*std::log(r);
            }
            else{
              uth = 0.;
              pr  = 3. + 4.*std::log(2.);
            }

            const double u   = -std::sin(phi)*uth;
            const double v   =  std::cos(phi)*uth;

            Complete state;

            state.density     = 1.;
            state.momentum[0] = u;
            state.momentum[1] = v;
            state.pressure    = pr;

            from_prim(state, equation_);

            Store(states, state, {i,j});
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
  geometry.cell_dimensions = std::array<int, Dim>{AMREX_D_DECL(64, 64, 1)};
  geometry.coordinates = amrex::RealBox({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                                        {AMREX_D_DECL(+1.0, +1.0, +1.0)});

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 6;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 2, 1)};

  using Complete = fub::PerfectGas<2>::Complete;
  fub::amrex::GradientDetector gradient(
      equation, std::make_pair(&Complete::density, 1.E-4),
      std::make_pair(&Complete::pressure, 1.E-4),
      std::make_pair([](const Complete& state) {
        return (state.momentum / state.density).matrix().norm();
      }, 1.E-4));

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.0;
  cons.momentum << 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> right;
  fub::CompleteFromCons(equation, right, cons);

  cons.energy *= 4;
  fub::Complete<fub::PerfectGas<2>> left;
  fub::CompleteFromCons(equation, left, cons);

  fub::amrex::BoundarySet boundary;
  using fub::amrex::TransmissiveBoundary;
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 0});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 1});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::Y, 0});
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::Y, 1});

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts),
                                                                             GreshoVortex{equation}, fub::amrex::TagAllOf(gradient, fub::amrex::TagBuffer(4)), boundary);
  gridding->InitializeHierarchy(0.0);

  auto tag = fub::execution::seq;

  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethod(tag, fub::MusclHancockMethod(equation)),
      fub::amrex::ForwardIntegrator(tag),
      fub::amrex::Reconstruction(tag, equation)};

  fub::DimensionalSplitLevelIntegrator solver(fub::int_c<2>, fub::amrex::IntegratorContext(gridding, method), fub::StrangSplitting());

  std::string base_name = "PerfectGas2d/";

  using namespace fub::amrex;
  auto output = [&](const std::shared_ptr<GriddingAlgorithm>& gridding,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}plt{:05}", base_name, cycle);
    amrex::Print() << "Start output to '" << name << "'.\n";
    WritePlotFile(name, gridding->GetPatchHierarchy(), equation);
    amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(), solver.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.final_time = 4s;
  run_options.output_interval = fub::Duration(1.0 / (run_options.final_time.count() * 30.0));
  run_options.cfl = 0.8;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
