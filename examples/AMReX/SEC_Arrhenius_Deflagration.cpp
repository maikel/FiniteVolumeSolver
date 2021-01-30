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

#include "fub/AMReX/cutcell/FluxMethodFactory.hpp"
#include "fub/equations/perfect_gas_mix/PlenaControl.hpp"

#include "fub/equations/PerfectGasMix.hpp"
#include "fub/flux_method/MusclHancockMethod2.hpp"

#include "fub/AMReX/DiffusionSourceTerm.hpp"
#include "fub/AMReX/boundary_condition/GenericPressureValveBoundary.hpp"
#include "fub/AMReX/boundary_condition/IsentropicPressureExpansion.hpp"
#include "fub/equations/perfect_gas_mix/ArrheniusKinetics.hpp"

#include <cmath>
#include <fmt/format.h>
#include <iostream>

#include <fenv.h>

#include <boost/log/utility/manipulators/add_value.hpp>

struct RiemannProblem {
  using Equation = fub::PerfectGasMix<1>;
  using Complete = fub::Complete<Equation>;
  using KineticState = fub::KineticState<Equation>;
  using Conservative = fub::Conservative<Equation>;

  Equation equation_;

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& /* grid */,
                      int /* level */, fub::Duration /*time*/) const {
    // const amrex::Geometry& geom =
    // grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& data = patch_level.data;
    fub::amrex::ForEachFab(
        fub::execution::openmp, data, [&](amrex::MFIter& mfi) {
          KineticState state(equation_);
          Complete complete(equation_);
          fub::Array<double, 1, 1> velocity{0.0};
          fub::View<Complete> states = fub::amrex::MakeView<Complete>(
              data[mfi], equation_, mfi.tilebox());
          fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
            // const double x = geom.CellCenter(int(i), 0);
            const double pressure = 1.0;
            state.temperature = 1.0;
            state.density = pressure / state.temperature / equation_.Rspec;
            state.mole_fractions[0] = 0.0;
            state.mole_fractions[1] = 1.0;
            // state.mole_fractions[2] = !(x < 0.4) ? 1.0 : 0.0;
            fub::euler::CompleteFromKineticState(equation_, complete, state,
                                                 velocity);
            fub::Store(states, complete, {i});
          });
        });
  }
};

struct ChangeTOpened {
  template <typename EulerEquation>
  [[nodiscard]] std::optional<fub::Duration>
  operator()(EulerEquation&, std::optional<fub::Duration>, double,
             const fub::perfect_gas_mix::gt::PlenumState&,
             const fub::amrex::GriddingAlgorithm& gridding, int) const
      noexcept {
    return gridding.GetTimePoint();
  }
};

struct IsNeverBlocked {
  template <typename EulerEquation>
  [[nodiscard]] bool
  operator()(EulerEquation&, std::optional<fub::Duration> /* t_opened */,
             double, const fub::perfect_gas_mix::gt::PlenumState&,
             const fub::amrex::GriddingAlgorithm& /* gridding */,
             int /* level */) const noexcept {
    return false;
  }
};

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard scope_guard{};

  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::PerfectGasMix<1> equation{};
  equation.n_species = 1;

  fub::amrex::CartesianGridGeometry grid_geometry(
      fub::GetOptions(options, "GridGeometry"));
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  using Complete = fub::PerfectGasMix<1>::Complete;
  fub::amrex::GradientDetector gradient{
      equation, std::make_pair(&Complete::density, 5e-3),
      std::make_pair(&Complete::pressure, 5e-2)};

  RiemannProblem initial_data{equation};

  fub::ProgramOptions hier_opts = fub::GetOptions(options, "PatchHierarchy");
  fub::amrex::PatchHierarchyOptions hierarchy_options(hier_opts);
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  fub::perfect_gas_mix::ArrheniusKinetics<1> source_term{
      equation, fub::GetOptions(options, "ArrheniusKinetics")};
  BOOST_LOG(log) << "ArrheniusKinetics:";
  source_term.options.Print(log);
  // const double eps = std::sqrt(std::numeric_limits<double>::epsilon());

  auto inflow_function =
      [prim = fub::Primitive<fub::PerfectGasMix<1>>(equation)](
          const fub::PerfectGasMix<1>& eq,
          fub::Complete<fub::PerfectGasMix<1>>& boundary_state,
          const fub::perfect_gas_mix::gt::PlenumState& compressor_state,
          double inner_pressure, fub::Duration, const amrex::MultiFab&,
          const fub::amrex::GriddingAlgorithm&, int) mutable {
        // fixed fuel concentration in deflagration mode
        const double X_inflow_left = 1.0;

        const double p_inflow_left = compressor_state.pressure;
        const double T_inflow_left = compressor_state.temperature;
        const double rho_inflow_left = p_inflow_left / T_inflow_left * eq.ooRspec;

        const double p = inner_pressure;
        const double ppv = p_inflow_left;
        const double rhopv = rho_inflow_left;
        const double Tpv = T_inflow_left;
        const double pin = p;
        const double Tin = Tpv * pow(pin / ppv, eq.gamma_minus_one_over_gamma);
        const double uin = std::sqrt(2.0 * eq.gamma_over_gamma_minus_one *
                                     std::max(0.0, Tpv - Tin));
        double rhoin = rhopv * pow(pin / ppv, eq.gamma_inv);

        FUB_ASSERT(rhoin > 0.0);
        FUB_ASSERT(pin > 0.0);
        prim.density = rhoin;
        prim.velocity[0] = uin;
        prim.pressure = pin;
        prim.species[0] = X_inflow_left;

        fub::CompleteFromPrim(eq, boundary_state, prim);
      };

  fub::ProgramOptions compressor_options = fub::GetOptions(options, "CompressorState");
  std::shared_ptr compressor_state = std::make_shared<fub::perfect_gas_mix::gt::PlenumState>();
  compressor_state->pressure = fub::GetOptionOr(compressor_options, "pressure", 1.0);
  compressor_state->temperature = fub::GetOptionOr(compressor_options, "temperature", 1.0);
  BOOST_LOG(log) << "CompressorState:";
  BOOST_LOG(log) << fmt::format("  - pressure = {}", compressor_state->pressure);
  BOOST_LOG(log) << fmt::format("  - temperature = {}", compressor_state->temperature);;

  using DeflagrationValve = fub::amrex::GenericPressureValveBoundary<
      fub::PerfectGasMix<1>, std::decay_t<decltype(inflow_function)>,
      ChangeTOpened, IsNeverBlocked>;
  DeflagrationValve valve(equation, compressor_state, inflow_function);

  fub::amrex::BoundarySet boundary;
  boundary.conditions.push_back(valve);
  using fub::amrex::IsentropicPressureExpansion;
  boundary.conditions.push_back(
      IsentropicPressureExpansion<fub::PerfectGasMix<1>>{equation, 1.0,
                                                         fub::Direction::X, 1});

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, grid_geometry, hierarchy_options),
      initial_data, gradient, boundary);
  gridding->InitializeHierarchy(0.0);

  using namespace std::literals;

  auto [flux_method, time_integrator] =
      fub::amrex::GetFluxMethod(fub::GetOptions(options, "FluxMethod"),
                                gridding->GetPatchHierarchy(), equation);
  fub::amrex::HyperbolicMethod method{flux_method, time_integrator,
                                      fub::amrex::Reconstruction(equation)};

  const int scratch_ghost_cell_width =
      fub::GetOptionOr(hier_opts, "scratch_gcw", 2);
  const int flux_ghost_cell_width =
      fub::GetOptionOr(hier_opts, "numeric_flux_gcw", 0);

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<1>,
      fub::amrex::IntegratorContext(gridding, method, scratch_ghost_cell_width,
                                    flux_ghost_cell_width),
      fub::GodunovSplitting());

  fub::amrex::DiffusionSourceTerm diffusion(equation, fub::GetOptions(options, "DiffusionSourceTerm"));
  diffusion.options_.Print(log);
  fub::SplitSystemSourceLevelIntegrator diffusive_integrator(
      std::move(level_integrator), std::move(diffusion),
      fub::StrangSplittingLumped());

  fub::SplitSystemSourceLevelIntegrator reactive_integrator(
      std::move(diffusive_integrator), std::move(source_term),
      fub::GodunovSplitting());
  // fub::StrangSplittingLumped());

  fub::SubcycleFineFirstSolver solver(std::move(reactive_integrator));

  using namespace fub::amrex;
  using namespace std::literals::chrono_literals;

  fub::OutputFactory<fub::amrex::GriddingAlgorithm> factory{};
  factory.RegisterOutput<fub::amrex::WriteHdf5>("HDF5");
  using CounterOutput = fub::CounterOutput<fub::amrex::GriddingAlgorithm,
                                           std::chrono::milliseconds>;
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> outputs(
      std::move(factory), fub::GetOptions(options, "Output"));

  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);

  outputs(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
  pybind11::scoped_interpreter interpreter{};
  std::optional<fub::ProgramOptions> opts = fub::ParseCommandLine(argc, argv);
  if (opts) {
    MyMain(*opts);
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}