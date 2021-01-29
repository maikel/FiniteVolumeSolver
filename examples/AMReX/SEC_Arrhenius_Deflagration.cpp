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

using FactoryFunction =
    std::function<fub::AnyFluxMethod<fub::amrex::IntegratorContext>(
        const fub::PerfectGasMix<1>&)>;

template <typename... Pairs> auto GetFluxMethodFactory(Pairs... ps) {
  std::map<std::string, FactoryFunction> factory;
  ((factory[ps.first] = ps.second), ...);
  return factory;
}

template <typename FluxMethod> struct MakeFlux {
  fub::AnyFluxMethod<fub::amrex::IntegratorContext>
  operator()(const fub::PerfectGasMix<1>& eq) const {
    FluxMethod flux_method{eq};
    fub::amrex::FluxMethodAdapter adapter(std::move(flux_method));
    return adapter;
  }
};

struct ChangeTOpened {
  template <typename EulerEquation>
  [[nodiscard]] std::optional<fub::Duration>
  operator()(EulerEquation&, std::optional<fub::Duration>, double,
             const fub::KineticState<EulerEquation>&,
             const fub::amrex::GriddingAlgorithm& gridding,
             int) const noexcept {
    return gridding.GetTimePoint();
  }
};

struct IsNeverBlocked {
  template <typename EulerEquation>
  [[nodiscard]] bool
  operator()(EulerEquation&, std::optional<fub::Duration> /* t_opened */,
             double, const fub::KineticState<EulerEquation>&,
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

  fub::ProgramOptions equation_options = fub::GetOptions(options, "Equation");

  fub::PerfectGasMix<1> equation{};
  equation.n_species = 1;
  equation.Rspec =
      fub::GetOptionOr(equation_options, "R_specific", equation.Rspec);
  equation.gamma = fub::GetOptionOr(equation_options, "gamma", equation.gamma);
  equation.gamma_minus_1_inv = 1.0 / (equation.gamma - 1.0);
  equation.gamma_array_ = fub::Array1d::Constant(equation.gamma);
  equation.gamma_minus_1_inv_array_ =
      fub::Array1d::Constant(equation.gamma_minus_1_inv);

  BOOST_LOG(log) << "Equation:";
  BOOST_LOG(log) << fmt::format(" - n_species = {}", equation.n_species);
  BOOST_LOG(log) << fmt::format(" - R_specific = {}", equation.Rspec);
  BOOST_LOG(log) << fmt::format(" - gamma = {}", equation.gamma);

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
          const fub::KineticState<fub::PerfectGasMix<1>>& compressor_state,
          double inner_pressure, fub::Duration, const amrex::MultiFab&,
          const fub::amrex::GriddingAlgorithm&, int) mutable {
        const double X_inflow_left =
            1.0; /* fixed fuel concentration in deflagration mode */

        const double p_inflow_left = fub::euler::Pressure(eq, compressor_state);
        const double rho_inflow_left = compressor_state.density;
        const double T_inflow_left = compressor_state.temperature;

        const double p = inner_pressure;
        const double ppv = p_inflow_left;
        const double rhopv = rho_inflow_left;
        const double Tpv = T_inflow_left;
        const double pin = p;
        const double g = eq.gamma;
        const double Gamma = (g - 1.0) / g;
        const double Gammainv = g / (g - 1.0);
        const double Tin = Tpv * pow(pin / ppv, Gamma);
        const double uin = std::sqrt(2.0 * Gammainv * std::max(0.0, Tpv - Tin));
        double rhoin = rhopv * pow(pin / ppv, 1.0 / g);

        prim.density = rhoin;
        prim.velocity[0] = uin;
        prim.pressure = pin;
        prim.species[0] = X_inflow_left;

        fub::CompleteFromPrim(eq, boundary_state, prim);
      };

  fub::KineticState<fub::PerfectGasMix<1>> compressor_state(equation);
  compressor_state.density = std::pow(2.0, 1.0 / equation.gamma);
  compressor_state.temperature = 2.0 / compressor_state.density;
  compressor_state.mole_fractions[0] = 1.0;

  using DeflagrationValve = fub::amrex::GenericPressureValveBoundary<
      fub::PerfectGasMix<1>, std::decay_t<decltype(inflow_function)>,
      ChangeTOpened, IsNeverBlocked>;
  DeflagrationValve valve(equation, compressor_state, inflow_function);

  /* fub::amrex::BoundarySet boundary;
  using fub::amrex::IsentropicPressureExpansion;
  fub::amrex::GenericPressureValveBoundary left(equation, inflow_function, {});
  boundary.conditions.push_back(left);
  boundary.conditions.push_back(
      IsentropicPressureExpansion<fub::PerfectGasMix<1>>{equation, 1.0,
                                                         fub::Direction::X, 1});
 */
  fub::amrex::BoundarySet boundary;
  boundary.conditions.push_back(valve);
  using fub::amrex::IsentropicPressureExpansion;
  // boundary.conditions.push_back(
  //    ReflectiveBoundary{seq, equation, fub::Direction::X, 1});
  boundary.conditions.push_back(
      IsentropicPressureExpansion<fub::PerfectGasMix<1>>{equation, 1.0,
                                                         fub::Direction::X, 1});

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, grid_geometry, hierarchy_options),
      initial_data, gradient, boundary);
  gridding->InitializeHierarchy(0.0);

  using namespace std::literals;
  using HLLEM = fub::perfect_gas::HllemMethod<fub::PerfectGasMix<1>>;

  using ConservativeReconstruction = fub::FluxMethod<fub::MusclHancock2<
      fub::PerfectGasMix<1>,
      fub::ConservativeGradient<
          fub::PerfectGasMix<1>,
          fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
      fub::ConservativeReconstruction<fub::PerfectGasMix<1>>, HLLEM>>;

  using PrimitiveReconstruction = fub::FluxMethod<fub::MusclHancock2<
      fub::PerfectGasMix<1>,
      fub::PrimitiveGradient<
          fub::PerfectGasMix<1>,
          fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
      fub::PrimitiveReconstruction<fub::PerfectGasMix<1>>, HLLEM>>;

  using CharacteristicsReconstruction = fub::FluxMethod<fub::MusclHancock2<
      fub::PerfectGasMix<1>,
      fub::CharacteristicsGradient<
          fub::PerfectGasMix<1>,
          fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
      fub::CharacteristicsReconstruction<fub::PerfectGasMix<1>>, HLLEM>>;

  auto flux_method_factory = GetFluxMethodFactory(
      std::pair{"NoReconstruct"s, MakeFlux<HLLEM>()},
      std::pair{"Conservative"s, MakeFlux<ConservativeReconstruction>()},
      std::pair{"Primitive"s, MakeFlux<PrimitiveReconstruction>()},
      std::pair{"Characteristics"s, MakeFlux<CharacteristicsReconstruction>()});

  std::string reconstruction =
      fub::GetOptionOr(options, "reconstruction", "HLLEM"s);
  BOOST_LOG(log) << "Reconstruction: " << reconstruction;
  auto flux_method = flux_method_factory.at(reconstruction)(equation);

  fub::amrex::HyperbolicMethod method{flux_method,
                                      fub::amrex::EulerForwardTimeIntegrator(),
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

  fub::amrex::DiffusionSourceTermOptions diffopt{};
  diffopt.mul = 0.0;
  fub::amrex::DiffusionSourceTerm diffusion(equation, diffopt);
  fub::SplitSystemSourceLevelIntegrator diffusive_integrator(
      std::move(level_integrator), std::move(diffusion),
      fub::StrangSplittingLumped());

  fub::SplitSystemSourceLevelIntegrator reactive_integrator(
      std::move(diffusive_integrator), std::move(source_term),
      fub::GodunovSplitting());
  // fub::StrangSplittingLumped());

  fub::SubcycleFineFirstSolver solver(std::move(reactive_integrator));
  // fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

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