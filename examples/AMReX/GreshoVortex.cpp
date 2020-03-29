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

struct GreshoVortex {
  using Equation = fub::PerfectGas<2>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  void InitializeData(fub::amrex::PatchLevel& patch_level, const fub::amrex::GriddingAlgorithm& grid, int level, fub::Duration /*time*/) const {
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& data = patch_level.data;
    fub::amrex::ForEachFab(data, [&](const amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());

      auto from_prim = [](Complete& state, const Equation& equation) {
        state.energy = state.pressure * equation.gamma_minus_1_inv +
                       0.5 * state.momentum.matrix().squaredNorm();
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

  const fub::amrex::ScopeGuard guard(argc, argv);
  fub::InitializeLogging(MPI_COMM_WORLD);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM >= 2);

  fub::PerfectGas<2> equation{};

  fub::amrex::CartesianGridGeometry geometry{};
  geometry.cell_dimensions = std::array<int, Dim>{AMREX_D_DECL(128, 128, 1)};
  geometry.coordinates = amrex::RealBox({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                                        {AMREX_D_DECL(+1.0, +1.0, +1.0)});
  geometry.periodicity = std::array<int, Dim>{AMREX_D_DECL(1, 1, 1)};

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 4;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 2, 1)};

  using Complete = fub::PerfectGas<2>::Complete;
  fub::amrex::GradientDetector gradient(
      equation, std::pair{&Complete::density, 1.0e-4},
      std::pair{&Complete::pressure, 1.0e-4},
      std::pair{[](const Complete& state) {
                  return (state.momentum / state.density).matrix().norm();
                },
                1.0e-4});

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts),
      GreshoVortex{equation},
      fub::amrex::TagAllOf(gradient, fub::amrex::TagBuffer(2)));
  gridding->InitializeHierarchy(0.0);

  fub::EinfeldtSignalVelocities<fub::PerfectGas<2>> signals{};
  fub::HllMethod hll_method(equation, signals);
  fub::MusclHancockMethod muscl_method(equation, hll_method);
  fub::amrex::HyperbolicMethod method{fub::amrex::FluxMethod(muscl_method),
                                      fub::amrex::EulerForwardTimeIntegrator(),
                                      fub::amrex::Reconstruction(equation)};

  const int scratch_gcw = 4 * muscl_method.GetStencilWidth();
  const int flux_gcw = 3 * muscl_method.GetStencilWidth();

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>,
      fub::amrex::IntegratorContext(gridding, method, scratch_gcw, flux_gcw),
      fub::StrangSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  std::string base_name = "GreshoVortex/";

  using namespace fub::amrex;
  using namespace std::literals::chrono_literals;
  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);
  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> output{};
  output.AddOutput(fub::MakeOutput<GriddingAlgorithm>(
      {1}, {fub::Duration(1.0 / 30.0)}, [&](const GriddingAlgorithm& gridding) {
        std::ptrdiff_t cycle = gridding.GetCycles();
        fub::Duration tp = gridding.GetTimePoint();
        BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", tp.count());
        std::string name = fmt::format("{}plt{:05}", base_name, cycle);
        BOOST_LOG(log) << "Start output to '" << name << "'.";
        WritePlotFile(name, gridding.GetPatchHierarchy(), equation);
        BOOST_LOG(log) << "Finished output to '" << name << "'.";
        double rho_max = 0.0;
        double rho_min = std::numeric_limits<double>::infinity();
        for (int level = 0;
             level < gridding.GetPatchHierarchy().GetNumberOfLevels();
             ++level) {
          const ::amrex::MultiFab& mf =
              gridding.GetPatchHierarchy().GetPatchLevel(level).data;
          rho_max = std::max(rho_max, mf.max(0));
          rho_min = std::min(rho_min, mf.min(0));
        }
        const double rho_err =
            std::max(std::abs(rho_max - 1.0), std::abs(rho_min - 1.0));
        BOOST_LOG(log) << fmt::format("Density Max Error: {:.6E}", rho_err);
      }));
  output.AddOutput(
      std::make_unique<fub::CounterOutput<fub::amrex::GriddingAlgorithm,
                                          std::chrono::milliseconds>>(
          wall_time_reference, std::vector<std::ptrdiff_t>{},
          std::vector<fub::Duration>{0.01s}));

  using namespace std::literals::chrono_literals;
  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.final_time = 3.0s;
  run_options.cfl = 0.8;
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}
