// Copyright (c) 2021 Maikel Nadolski
// Copyright (c) 2021 Christian Zenker
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
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX_CutCell.hpp"
#include "fub/Solver.hpp"

#include "fub/AMReX/boundary_condition/GenericPressureValveBoundary.hpp"
#include "fub/equations/PerfectGasMix.hpp"
#include "fub/equations/perfect_gas_mix/ArrheniusKinetics.hpp"
#include "fub/equations/perfect_gas_mix/PlenaControl.hpp"

#include "fub/ext/CopyInputFile.hpp"

#include "fub/flux_method/MusclHancockMethod2.hpp"

#include "fub/AMReX/AxialFluxMethodAdapter.hpp"
#include "fub/AMReX/AxialTimeIntegrator.hpp"
#include "fub/AMReX/DiffusionSourceTerm.hpp"

#include "fub/AMReX/multi_block/MultiBlockBoundary2.hpp"
#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm2.hpp"
#include "fub/AMReX/multi_block/MultiBlockIntegratorContext2.hpp"

#include "fub/AMReX/cutcell/FluxMethodFactory.hpp"
#include "fub/AMReX/cutcell/boundary_condition/IsentropicPressureExpansion.hpp"
#include "fub/AMReX/cutcell/boundary_condition/MachnumberBoundary.hpp"
#include "fub/AMReX/cutcell/boundary_condition/ReflectiveBoundary2.hpp"
#include "fub/AMReX/cutcell/boundary_condition/TurbineMassflowBoundary.hpp"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <boost/filesystem.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Translation.H>

#include <cmath>
#include <iostream>

#include <range/v3/view/enumerate.hpp>
#include <xmmintrin.h>

namespace GT = fub::perfect_gas_mix::gt;

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 2;

static_assert(AMREX_SPACEDIM == 2);

static constexpr double r_tube = 0.015;

static constexpr int n_species = 2;
static constexpr int n_passive_scalars = 1;

struct TracePassiveScalarBoundary {
  using Equation = fub::PerfectGasMix<1>;
  using Complete = Equation::Complete;
  using Conservative = Equation::Conservative;

  Equation equation_;
  fub::Duration dt_;

  void FillBoundary(::amrex::MultiFab& mf,
                    const fub::amrex::GriddingAlgorithm& grid, int level) {

    dt_ = grid.GetPatchHierarchy().GetStabledt();

    using HLLEM_Lar = fub::perfect_gas::HllemMethod<Equation>;

    using CharacteristicsReconstruction = fub::FluxMethod<fub::MusclHancock2<
        Equation,
        fub::CharacteristicsGradient<
            Equation, fub::CentralDifferenceGradient<fub::MinModLimiter>>,
        fub::CharacteristicsReconstruction<Equation>, HLLEM_Lar>>;

    CharacteristicsReconstruction flux_method(equation_);
    Conservative flux(equation_);
    std::array<Complete, 4> stencil{};
    stencil.fill(Complete(equation_));

    int ngrow = 2;
    const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    const double dx = geom.CellSize(0);
    const ::amrex::IntVect lo{-2 * ngrow, 0};
    const ::amrex::IntVect hi{+ngrow, 0};
    const ::amrex::Box stencil_box{lo, hi};
    fub::amrex::ForEachFab(mf, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::Box section = mfi.growntilebox() & stencil_box;
      if (section == stencil_box) {
        fub::View<Complete> stencilv =
            fub::amrex::MakeView<Complete>(mf[mfi], equation_, stencil_box);
        for (int i = 0; i < int(stencil.size()); ++i) {
          fub::Load(stencil[i], stencilv, fub::Index<1>{i - ngrow});
        }
        flux_method.ComputeNumericFlux(flux, stencil, dt_, dx,
                                       fub::Direction::X);
        const double X_right =
            stencilv.passive_scalars(0, 0) / stencilv.density(0);
        const double X_left = X_right + dt_.count() / dx;
        stencilv.passive_scalars(-1, 0) = X_left * stencilv.density(-1);
        stencilv.passive_scalars(-2, 0) = X_left * stencilv.density(-2);
        stencilv.passive_scalars(-3, 0) = X_left * stencilv.density(-3);
        stencilv.passive_scalars(-4, 0) = X_left * stencilv.density(-4);
      } else if (!section.isEmpty()) {
        fub::View<Complete> stencilv =
            fub::amrex::MakeView<Complete>(mf[mfi], equation_, section);
        fub::ForEachIndex(fub::Box<0>(stencilv), [&](std::ptrdiff_t i) {
          stencilv.passive_scalars(i, 0) = stencilv.passive_scalars(0, 0);
        });
      }
    });
  }

  void FillBoundary(::amrex::MultiFab& mf,
                    const fub::amrex::GriddingAlgorithm& grid, int level,
                    fub::Direction dir) {
    if (dir == fub::Direction::X) {
      FillBoundary(mf, grid, level);
    }
  }
};

struct InitialDataInTube {
  using Equation = fub::PerfectGasMix<1>;
  using Complete = fub::Complete<Equation>;
  using KineticState = fub::KineticState<Equation>;
  using Conservative = fub::Conservative<Equation>;

  Equation equation_;
  fub::perfect_gas_mix::ArrheniusKineticsOptions kinetics_;
  fub::perfect_gas_mix::gt::ControlOptions control_;
  double x_0_;
  double initially_filled_x_{0.4};

  double CompressorPressureRatio(double rpm) const noexcept {
    double pratio = control_.pratiomean +
                    control_.c_0 * control_.pratiovar *
                        std::atan(M_PI * (rpm / control_.rpmmean - 1.0)) / M_PI;
    return pratio;
  }

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration /*time*/) const {
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& data = patch_level.data;
    fub::amrex::ForEachFab(fub::execution::seq, data, [&](amrex::MFIter& mfi) {
      Conservative state(equation_);
      Complete complete(equation_);
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
        const double x = geom.CellCenter(int(i), 0);
        const double rel_x = x - x_0_;

        const double pV = CompressorPressureRatio(control_.rpmmin);
        const double TV =
            1.0 + (std::pow(pV, equation_.gamma_minus_one_over_gamma) - 1.0) /
                      control_.efficiency_compressor;
        const double p0 = 2.0;
        const double p = 0.95 * p0;
        const double T = TV + kinetics_.Q * equation_.gamma_minus_one;
        const double rho = p / T; // * equation_.ooRspec;

        state.density = rho;
        state.species[0] =
            (rel_x < initially_filled_x_) ? 1.0 * rho : 0.0 * rho;
        state.passive_scalars[0] = -x * rho;
        state.energy = p * equation_.gamma_minus_one_inv;
        equation_.CompleteFromCons(complete, state);
        fub::Store(states, complete, {i});
      });
    });
  }
};

struct InitialDataInPlenum {
  using Equation = fub::PerfectGasMix<2>;
  using Complete = fub::Complete<Equation>;
  using KineticState = fub::KineticState<Equation>;
  using Conservative = fub::Conservative<Equation>;

  Equation equation_;
  fub::perfect_gas_mix::ArrheniusKineticsOptions kinetics_;

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::cutcell::GriddingAlgorithm& grid,
                      int level, fub::Duration /*time*/) const {
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& data = patch_level.data;
    fub::amrex::ForEachFab(
        fub::execution::openmp, data, [&](amrex::MFIter& mfi) {
          KineticState state(equation_);
          Complete complete(equation_);
          fub::Array<double, 2, 1> velocity{0.0, 0.0};
          fub::View<Complete> states = fub::amrex::MakeView<Complete>(
              data[mfi], equation_, mfi.tilebox());
          fub::ForEachIndex(
              fub::Box<0>(states), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
                const double x = geom.CellCenter(int(i), 0);
                // const double pressure = 1.0;
                const double p0 = 2.0;
                // const double rho0 = std::pow(p0, equation_.gamma_inv);
                // const double T0 = p0 / rho0;
                // const double p = 0.95 * p0;
                const double T = 9.0; // T0 + kinetics_.Q * equation_.gamma_minus_one;
                const double rho = p0 / T; // * equation_.ooRspec;

                state.temperature = T;
                state.density = rho;
                state.mole_fractions[0] = 0.0;
                state.mole_fractions[1] = 1.0;
                state.passive_scalars[0] = -x;
                fub::euler::CompleteFromKineticState(equation_, complete, state,
                                                     velocity);
                fub::Store(states, complete, {i, j});
              });
        });
  }
};

auto MakeTubeSolver(
    const fub::ProgramOptions& options,
    const fub::PerfectGasConstants& constants,
    const std::shared_ptr<fub::CounterRegistry>& counters,
    const std::shared_ptr<const GT::ControlState>& control_state,
    GT::ControlOptions& control_options) {
  using namespace fub::amrex;
  auto make_tube_solver = counters->get_timer("MakeTubeSolver");

  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << "==================== Tube =========================";

  CartesianGridGeometry grid_geometry(fub::GetOptions(options, "GridGeometry"));
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  PatchHierarchyOptions hierarchy_options(
      fub::GetOptions(options, "PatchHierarchy"));
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  using Eq = fub::PerfectGasMix<Tube_Rank>;
  using Complete = Eq::Complete;
  fub::PerfectGasMix<Tube_Rank> equation{constants, n_species - 1,
                                         n_passive_scalars};

  GradientDetector gradient{equation, std::make_pair(&Complete::density, 1e-2),
                            std::make_pair(&Complete::pressure, 1e-2)};

  DataDescription desc = MakeDataDescription(equation);

  ::amrex::Box refine_box{{grid_geometry.cell_dimensions[0] - 5, 0},
                          {grid_geometry.cell_dimensions[0] - 1, 0}};
  ConstantBox constant_box{refine_box};

  fub::perfect_gas_mix::ArrheniusKinetics<1> source_term{
      equation, fub::GetOptions(options, "ArrheniusKinetics")};
  BOOST_LOG(log) << "ArrheniusKinetics:";
  source_term.options.Print(log);

  const double initially_filled_x =
      fub::GetOptionOr(options, "initially_filled_x", 0.4);
  BOOST_LOG(log) << "InitialData:";
  BOOST_LOG(log) << "  - initially_filled_x = " << initially_filled_x << " [m]";
  InitialDataInTube initial_data{equation, source_term.options, control_options,
                                 grid_geometry.coordinates.lo()[0],
                                 initially_filled_x};

  FUB_ASSERT(control_state);

  //!!!!!!!!!!!!!!// begin Deflagrationvalve
  /*

  auto combustor_inflow_function =
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
        const double rho_inflow_left =
            p_inflow_left / T_inflow_left;// * eq.ooRspec;

        const double p = inner_pressure;
        const double ppv = p_inflow_left;
        const double rhopv = rho_inflow_left;
        const double Tpv = T_inflow_left;
        const double pin = p;
        const double Tin = Tpv * pow(pin / ppv,
        eq.gamma_minus_one_over_gamma); const double uin = std::sqrt(2.0 *
        eq.gamma_over_gamma_minus_one *
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

    fub::amrex::DeflagrationValve valve(equation, control_state,
    combustor_inflow_function);

  */
  //!!!!!!!!!!!!!!// End Defagrationvalve

  //!!!!!!!!!!!!!!// Begin PressureValveBoundary_Klein
  /*

  const  double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  auto inflow_function =
      [eps, source_term,
       prim = fub::Primitive<fub::PerfectGasMix<1>>(equation)](
          const fub::PerfectGasMix<1>& eq,
          fub::Complete<fub::PerfectGasMix<1>>& boundary_state,
          const fub::perfect_gas_mix::gt::PlenumState& compressor_state, double
  inner_pressure, fub::Duration t_diff, const amrex::MultiFab&, const
  fub::amrex::GriddingAlgorithm&, int) mutable { const double fuel_retardatation
  = 0.06;                // Reference: 0.06;   for icx = 256:  0.1 // const
  double tti = 0.75; // Reference: 0.75; // const double timin = 0.1; const
  double X_inflow_left = 1.0;

        const double p_inflow_left = compressor_state.pressure;
        const double T_inflow_left = compressor_state.temperature;
        const double rho_inflow_left =
            p_inflow_left / T_inflow_left;// * eq.ooRspec;

        const double p = inner_pressure;
        const double ppv = p_inflow_left;
        const double rhopv = rho_inflow_left;
        const double Tpv = T_inflow_left;
        const double pin = p;
        const double Tin = Tpv * pow(pin / ppv,
        eq.gamma_minus_one_over_gamma); const double uin = std::sqrt(2.0 *
        eq.gamma_over_gamma_minus_one *
                                     std::max(0.0, Tpv - Tin));
        double rhoin = rhopv * pow(pin / ppv, eq.gamma_inv);

        const double tign = std::max(timin, tti - t_diff.count());
        const double Xin = X_inflow_left;
        const double Tin1 =
        fub::perfect_gas_mix::TemperatureForIgnitionDelay(
            eq, source_term.options, tign, Xin, Tin, eps);

        // adjust density to match desired temperature //
        rhoin *= Tin / Tin1;

        prim.density = rhoin;
        prim.velocity[0] = uin;
        prim.pressure = pin;
        auto heaviside = [](double x) { return (x > 0); };
        prim.species[0] = std::clamp(
            Xin * heaviside(t_diff.count() - fuel_retardatation), 0.0, 1.0);

        fub::CompleteFromPrim(eq, boundary_state, prim);
      };

  fub::amrex::PressureValveBoundary_Klein valve(equation, control_state,
                                                inflow_function);

  */

  //!!!!!!!!!!!!!!// End PressureValveBoundary_Klein

  //!!!!!!!!!!!!!!// Begin SwitchDeflagrationToSECValve

  fub::ProgramOptions SEC_options =
      fub::GetOptions(options, "InflowOptionsSEC");
  double SEC_buffer = fub::GetOptionOr(SEC_options, "SEC_buffer", 0.06);
  double SEC_tti = fub::GetOptionOr(SEC_options, "SEC_tti", 1.2);
  double SEC_timin = fub::GetOptionOr(SEC_options, "SEC_timin", 0.1);
  BOOST_LOG(log) << "InflowOptionsSEC:";
  BOOST_LOG(log) << "  - SEC_buffer = " << SEC_buffer;
  BOOST_LOG(log) << "  - SEC_tti = " << SEC_tti;
  BOOST_LOG(log) << "  - SEC_timin = " << SEC_timin;

  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  auto switch_to_SEC_inflow_function =
      [eps, source_term, SEC_buffer, SEC_tti, SEC_timin,
       prim = fub::Primitive<fub::PerfectGasMix<1>>(equation)](
          const fub::PerfectGasMix<1>& eq,
          fub::Complete<fub::PerfectGasMix<1>>& boundary_state,
          const fub::perfect_gas_mix::gt::PlenumState& compressor_state,
          double inner_pressure, fub::Duration t_diff, const amrex::MultiFab&,
          const fub::amrex::GriddingAlgorithm&, int) mutable {
        if (!compressor_state.SEC_Mode) {
          // Deflagration Mode
          const double X_inflow_left = 1.0;

          const double p_inflow_left = compressor_state.pressure;
          const double T_inflow_left = compressor_state.temperature;
          const double rho_inflow_left =
              p_inflow_left / T_inflow_left; // * eq.ooRspec;

          const double p = inner_pressure;
          const double ppv = p_inflow_left;
          const double rhopv = rho_inflow_left;
          const double Tpv = T_inflow_left;
          const double pin = p;
          const double Tin =
              Tpv * pow(pin / ppv, eq.gamma_minus_one_over_gamma);
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
        } else {
          // SEC Mode
          const double fuel_retardatation =
              SEC_buffer; // Reference: 0.06;   for icx = 256:  0.1//
          const double tti = SEC_tti;     // Reference: 0.75; best 1.2 //
          const double timin = SEC_timin; // Reference: 0.1 //
          const double X_inflow_left = 1.0;

          const double p_inflow_left = compressor_state.pressure;
          const double T_inflow_left = compressor_state.temperature;
          const double rho_inflow_left =
              p_inflow_left / T_inflow_left; // * eq.ooRspec;

          const double p = inner_pressure;
          const double ppv = p_inflow_left;
          const double rhopv = rho_inflow_left;
          const double Tpv = T_inflow_left;
          const double pin = p;
          const double Tin =
              Tpv * pow(pin / ppv, eq.gamma_minus_one_over_gamma);
          const double uin = std::sqrt(2.0 * eq.gamma_over_gamma_minus_one *
                                       std::max(0.0, Tpv - Tin));
          double rhoin = rhopv * pow(pin / ppv, eq.gamma_inv);

          const double tign = std::max(timin, tti - t_diff.count());
          const double Xin = X_inflow_left;
          const double Tin1 = fub::perfect_gas_mix::TemperatureForIgnitionDelay(
              eq, source_term.options, tign, Xin, Tin, eps);

          // adjust density to match desired temperature //
          rhoin *= Tin / Tin1;

          prim.density = rhoin;
          prim.velocity[0] = uin;
          prim.pressure = pin;
          auto heaviside = [](double x) { return (x > 0); };
          prim.species[0] = std::clamp(
              Xin * heaviside(t_diff.count() - fuel_retardatation), 0.0, 1.0);

          fub::CompleteFromPrim(eq, boundary_state, prim);
        }
      };

  fub::amrex::SwitchDeflagrationToSECValve valve(equation, control_state,
                                                 switch_to_SEC_inflow_function);

  //!!!!!!!!!!!!!!// End SwitchDeflagrationToSECValve

  BoundarySet boundaries{{valve}};

  // If a checkpoint path is specified we will fill the patch hierarchy with
  // data from the checkpoint file, otherwise we will initialize the data by
  // the initial data function.
  std::shared_ptr<GriddingAlgorithm> gridding = [&] {
    std::string checkpoint =
        fub::GetOptionOr(options, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      BOOST_LOG(log) << "Initialize grid by initial condition...";
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(desc, grid_geometry, hierarchy_options), initial_data,
          TagAllOf(gradient, constant_box), boundaries);
      gridding->GetPatchHierarchy().SetCounterRegistry(counters);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
      PatchHierarchy h = ReadCheckpointFile(checkpoint, desc, grid_geometry,
                                            hierarchy_options);
      std::shared_ptr<GriddingAlgorithm> gridding =
          std::make_shared<GriddingAlgorithm>(std::move(h), initial_data,
                                              TagAllOf(gradient, constant_box),
                                              boundaries);
      return gridding;
    }
  }();

  auto [flux_method, time_integrator] =
      fub::amrex::GetFluxMethod(fub::GetOptions(options, "FluxMethod"),
                                gridding->GetPatchHierarchy(), equation);
  HyperbolicMethod method{flux_method, time_integrator,
                          Reconstruction(equation)};

  IntegratorContext context(gridding, method,
                            fub::GetOptions(options, "IntegratorContext"));
  BOOST_LOG(log) << "IntegratorContext:";
  context.GetOptions().Print(log);

  fub::Duration good_guess_dt = context.GetPatchHierarchy().GetStabledt();
  TracePassiveScalarBoundary passive_scalar_boundary{equation, good_guess_dt};
  boundaries.conditions.push_back(std::move(passive_scalar_boundary));

  context.GetGriddingAlgorithm()->GetBoundaryCondition() = boundaries;

  BOOST_LOG(log) << "==================== End Tube =========================";

  return std::pair{context, source_term};
}

auto MakePlenumSolver(
    const std::map<std::string, pybind11::object>& options,
    const fub::PerfectGasConstants& constants,
    const fub::perfect_gas_mix::ArrheniusKineticsOptions& kinetics) {
  std::shared_ptr<fub::CounterRegistry> registry =
      std::make_shared<fub::CounterRegistry>();
  auto make_plenum_solver = registry->get_timer("MakePlenumSolver");
  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << "==================== Plenum =========================";
  using namespace fub::amrex::cutcell;
  auto MakePolygon = [](auto... points) {
    fub::Polygon::Vector xs{};
    fub::Polygon::Vector ys{};

    (xs.push_back(std::get<0>(points)), ...);
    (ys.push_back(std::get<1>(points)), ...);

    return fub::Polygon(std::move(xs), std::move(ys));
  };

  std::vector<fub::PolymorphicGeometry<Plenum_Rank>> inlets{};
  std::vector<pybind11::dict> eb_dicts{};
  eb_dicts = fub::GetOptionOr(options, "InletGeometries", eb_dicts);
  for (pybind11::dict& dict : eb_dicts) {
    fub::ProgramOptions inlet_options = fub::ToMap(dict);
    double r_inlet_start = r_tube;
    double r_inlet_end = 2.0 * r_tube;
    double y_0 = 0.0;
    r_inlet_start = fub::GetOptionOr(inlet_options, "r_start", r_inlet_start);
    r_inlet_end = fub::GetOptionOr(inlet_options, "r_end", r_inlet_end);
    y_0 = fub::GetOptionOr(inlet_options, "y_0", y_0);
    static constexpr double height = 2.0;
    const double xlo = -height;
    const double xhi = 0.0;
    const double r = r_inlet_start;
    const double r2 = r_inlet_end;
    const double xdiv = xhi - 4.0 * r;
    auto polygon =
        MakePolygon(std::pair{xlo, y_0 + r}, std::pair{xdiv, y_0 + r},
                    std::pair{xhi, y_0 + r2}, std::pair{xhi, y_0 - r2},
                    std::pair{xdiv, y_0 - r}, std::pair{xlo, y_0 - r},
                    std::pair{xlo, y_0 + r});
    inlets.push_back(polygon);
  }
  fub::PolymorphicUnion<Plenum_Rank> union_of_inlets(std::move(inlets));

  auto embedded_boundary = amrex::EB2::makeIntersection(
      amrex::EB2::PlaneIF({0.0, 0.0}, {1.0, 0.0}, false),
      fub::amrex::Geometry(fub::Invert(union_of_inlets)));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::amrex::CartesianGridGeometry grid_geometry =
      fub::GetOptions(options, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  PatchHierarchyOptions hierarchy_options =
      fub::GetOptions(options, "PatchHierarchy");
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  BOOST_LOG(log) << "Compute EB level set data...";
  {
    fub::Timer timer = registry->get_timer("MakeIndexSpaces");
    hierarchy_options.index_spaces =
        MakeIndexSpaces(shop, grid_geometry, hierarchy_options);
  }

  fub::PerfectGasMix<Plenum_Rank> equation{constants, n_species - 1,
                                           n_passive_scalars};

  InitialDataInPlenum initial_data{equation, kinetics};

  //  using Complete = fub::Complete<fub::PerfectGas<Plenum_Rank>>;
  using Complete = fub::Complete<fub::PerfectGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&Complete::pressure, 0.05},
                             std::pair{&Complete::density, 0.05}};

  fub::amrex::IntegratorContextOptions context_options =
      fub::GetOptions(options, "IntegratorContext");

  const int scratch_gcw = context_options.scratch_gcw;

  const int lower_left_corner_x0 = -scratch_gcw;
  const int lower_left_corner_y0 = -scratch_gcw;
  ::amrex::IntVect lower_left_corner_lo{lower_left_corner_x0,
                                        lower_left_corner_y0};
  ::amrex::IntVect lower_left_corner_hi{lower_left_corner_x0 + scratch_gcw - 1,
                                        lower_left_corner_y0 + scratch_gcw - 1};
  ::amrex::Box lower_left_corner{lower_left_corner_lo, lower_left_corner_hi};

  const int upper_left_corner_x0 = -scratch_gcw;
  const int upper_left_corner_y0 = grid_geometry.cell_dimensions[1];
  ::amrex::IntVect upper_left_corner_lo{upper_left_corner_x0,
                                        upper_left_corner_y0};
  ::amrex::IntVect upper_left_corner_hi{upper_left_corner_x0 + scratch_gcw - 1,
                                        upper_left_corner_y0 + scratch_gcw - 1};
  ::amrex::Box upper_left_corner{upper_left_corner_lo, upper_left_corner_hi};

  const int lower_right_corner_x0 = grid_geometry.cell_dimensions[0];
  const int lower_right_corner_y0 = -scratch_gcw;
  ::amrex::IntVect lower_right_corner_lo{lower_right_corner_x0,
                                         lower_right_corner_y0};
  ::amrex::IntVect lower_right_corner_hi{
      lower_right_corner_x0 + scratch_gcw - 1,
      lower_right_corner_y0 + scratch_gcw - 1};
  ::amrex::Box lower_right_corner{lower_right_corner_lo, lower_right_corner_hi};

  const int upper_right_corner_x0 = grid_geometry.cell_dimensions[0];
  const int upper_right_corner_y0 = grid_geometry.cell_dimensions[1];
  ::amrex::IntVect upper_right_corner_lo{upper_right_corner_x0,
                                         upper_right_corner_y0};
  ::amrex::IntVect upper_right_corner_hi{
      upper_right_corner_x0 + scratch_gcw - 1,
      upper_right_corner_y0 + scratch_gcw - 1};
  ::amrex::Box upper_right_corner{upper_right_corner_lo, upper_right_corner_hi};

  BoundarySet boundary_condition{
      {ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::X, 0},
       ReflectiveBoundary2{equation, fub::Direction::X, 0, lower_left_corner},
       ReflectiveBoundary2{equation, fub::Direction::X, 0, upper_left_corner},
       ReflectiveBoundary2{equation, fub::Direction::X, 1, lower_right_corner},
       ReflectiveBoundary2{equation, fub::Direction::X, 1, upper_right_corner},
       ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::X, 1},
       ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::Y, 0},
       ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::Y,
                          1}}};

  TurbineMassflowBoundaryOptions boundary_options =
      fub::GetOptions(options, "TurbineMassflowBoundaries");
  fub::amrex::cutcell::TurbineMassflowBoundaryOptions tb_opts(boundary_options);
  tb_opts.dir = fub::Direction::X;
  tb_opts.side = 1;
  BOOST_LOG(log) << "TurbineMassflowBoundaries:";
  tb_opts.Print(log);
  fub::amrex::cutcell::TurbineMassflowBoundary<fub::PerfectGasMix<2>>
      pressure_outflow(equation, tb_opts);
  boundary_condition.conditions.push_back(std::move(pressure_outflow));

  ::amrex::RealBox xbox = grid_geometry.coordinates;
  ::amrex::Geometry coarse_geom = fub::amrex::GetCoarseGeometry(grid_geometry);

  ::amrex::RealBox inlet{{xbox.lo(0), -0.15}, {0.18, +0.15}};
  ::amrex::Box refine_box = fub::amrex::BoxWhichContains(inlet, coarse_geom);
  ConstantBox constant_refinebox{refine_box};

  std::shared_ptr gridding = [&] {
    fub::SeverityLogger log = fub::GetInfoLogger();
    std::string checkpoint =
        fub::GetOptionOr(options, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      BOOST_LOG(log) << "Initialize grid by initial condition...";
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(equation, grid_geometry, hierarchy_options),
          initial_data, TagAllOf(TagCutCells(), constant_refinebox),
          boundary_condition);
      grid->GetPatchHierarchy().SetCounterRegistry(registry);
      grid->InitializeHierarchy(0.0);
      return grid;
    } else {
      checkpoint += "/Plenum";
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), grid_geometry,
          hierarchy_options);
      h.SetCounterRegistry(registry);
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), constant_refinebox), boundary_condition);
      return grid;
    }
  }();

  auto flux_method = fub::amrex::cutcell::GetCutCellMethod(
      fub::GetOptions(options, "FluxMethod"), equation);

  HyperbolicMethod method{flux_method, TimeIntegrator{},
                          Reconstruction{equation}};

  IntegratorContext context(gridding, method, context_options);
  BOOST_LOG(log) << "IntegratorContext:";
  context.GetOptions().Print(log);

  BOOST_LOG(log) << "==================== End Plenum =========================";

  return context;
}

void MyMain(const std::map<std::string, pybind11::object>& vm);

int main(int argc, char** argv) {
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  int provided{-1};
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
  if (provided < MPI_THREAD_FUNNELED) {
    fmt::print(
        "Aborting execution. MPI could not provide a thread-safe instance.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  pybind11::scoped_interpreter interpreter{};
  std::optional<fub::ProgramOptions> opts = fub::ParseCommandLine(argc, argv);
  if (opts) {
    fub::CopyInputFile(MPI_COMM_WORLD,
                       fub::GetOptions(*opts, "InputFileOptions"), argc, argv);
    fub::InitializeLogging(MPI_COMM_WORLD,
                           fub::GetOptions(*opts, "LogOptions"));
    MyMain(*opts);
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}

void WriteCheckpoint(
    const std::string& path,
    const fub::amrex::MultiBlockGriddingAlgorithm2& grid,
    const std::shared_ptr<const fub::perfect_gas_mix::gt::ControlState>&
        control_state) {
  auto tubes = grid.GetTubes();
  for (auto&& [i, tube] : ranges::view::enumerate(tubes)) {
    std::string name = fmt::format("{}/Tube_{}", path, i);
    fub::amrex::WriteCheckpointFile(name, tube->GetPatchHierarchy());
  }
  std::string name = fmt::format("{}/Plenum", path);
  fub::amrex::cutcell::WriteCheckpointFile(
      name, grid.GetPlena()[0]->GetPatchHierarchy());
  int rank = -1;
  ::MPI_Comm_rank(amrex::ParallelDescriptor::Communicator(), &rank);
  if (rank == 0) {
    name = fmt::format("{}/ControlState", path);
    std::ofstream valve_checkpoint(name);
    boost::archive::text_oarchive oa(valve_checkpoint);
    oa << *control_state;
  }
}

void MyMain(const std::map<std::string, pybind11::object>& vm) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};

  std::vector<fub::amrex::cutcell::IntegratorContext> plenum{};
  std::vector<fub::amrex::IntegratorContext> tubes{};
  using Kinetics = fub::perfect_gas_mix::ArrheniusKinetics<Tube_Rank>;
  std::vector<Kinetics> kinetics{};
  std::vector<fub::amrex::BlockConnection> connectivity{};

  fub::ProgramOptions eq_options = fub::GetOptions(vm, "Equation");
  double Rspec = fub::GetOptionOr(eq_options, "Rspec", 1.0);
  double gamma = fub::GetOptionOr(eq_options, "gamma", 1.28);
  fub::PerfectGasConstants constants{Rspec, gamma};
  fub::SeverityLogger log = fub::GetInfoLogger();

  BOOST_LOG(log) << "Equation:";
  BOOST_LOG(log) << fmt::format("  - Rspec = {}", constants.Rspec);
  BOOST_LOG(log) << fmt::format("  - gamma = {}", constants.gamma);

  fub::PerfectGasMix<Plenum_Rank> plenum_equation{constants, n_species - 1,
                                                  n_passive_scalars};
  fub::PerfectGasMix<Tube_Rank> tube_equation{constants, n_species - 1,
                                              n_passive_scalars};
  fub::ProgramOptions control_options_map =
      fub::GetOptions(vm, "ControlOptions");
  GT::ControlOptions control_options(control_options_map);
  BOOST_LOG(log) << "ControlOptions:";
  control_options.Print(log);
  GT::Control control(plenum_equation, control_options);
  std::shared_ptr<GT::ControlState> control_state = control.GetSharedState();
  std::string checkpoint =
      fub::GetOptionOr(control_options_map, "checkpoint", std::string{});
  if (!checkpoint.empty()) {
    std::string input =
        fub::ReadAndBroadcastFile(checkpoint + "/ControlState",
                                  ::amrex::ParallelDescriptor::Communicator());
    std::istringstream ifs(input);
    boost::archive::text_iarchive ia(ifs);
    ia >> *control_state;
  }

  fub::perfect_gas_mix::ArrheniusKineticsOptions kinetic_opts =
      fub::GetOptions(vm, "ArrheniusKinetics");
  plenum.push_back(
      MakePlenumSolver(fub::GetOptions(vm, "Plenum"), constants, kinetic_opts));
  auto counter_database = plenum[0].GetCounterRegistry();

  std::vector<pybind11::dict> tube_dicts = {};
  tube_dicts = fub::GetOptionOr(vm, "Tubes", tube_dicts);
  for (pybind11::dict& dict : tube_dicts) {
    fub::ProgramOptions tube_options = fub::ToMap(dict);
    auto&& [tube, source] =
        MakeTubeSolver(tube_options, constants, counter_database, control_state,
                       control_options);
    fub::amrex::BlockConnection connection;
    connection.direction = fub::Direction::X;
    connection.side = 0;
    connection.ghost_cell_width = std::max(tube.GetOptions().scratch_gcw,
                                           plenum[0].GetOptions().scratch_gcw);
    connection.plenum.id = 0;
    connection.tube.id = tubes.size();
    connection.tube.mirror_box = tube.GetGriddingAlgorithm()
                                     ->GetPatchHierarchy()
                                     .GetGeometry(0)
                                     .Domain();
    amrex::Box plenum_mirror_box{};
    FUB_ASSERT(!plenum_mirror_box.ok());
    plenum_mirror_box =
        fub::GetOptionOr(tube_options, "plenum_mirror_box", plenum_mirror_box);
    if (!plenum_mirror_box.ok()) {
      throw std::runtime_error(
          "You need to specify plenum_mirror_box for each tube.");
    }
    connection.plenum.mirror_box = plenum_mirror_box;
    tubes.push_back(std::move(tube));
    kinetics.push_back(source);
    connectivity.push_back(connection);
  }

  fub::amrex::MultiBlockIntegratorContext2 context(
      tube_equation, plenum_equation, std::move(tubes), std::move(plenum),
      std::move(connectivity));

  // get the coarse_average_mirror_box from the TurbineMassflowBoundary, only
  // there we have an mass flow outside of the domain
  const fub::ProgramOptions plenum_options = fub::GetOptions(vm, "Plenum");
  fub::amrex::cutcell::TurbineMassflowBoundaryOptions boundary_options =
      fub::GetOptions(plenum_options, "TurbineMassflowBoundaries");
  ::amrex::Box coarse_average_mirror_box =
      boundary_options.coarse_average_mirror_box;

  GT::ControlFeedback<Plenum_Rank> feedback(plenum_equation, tube_equation,
                                            control, coarse_average_mirror_box);
  context.SetPostAdvanceHierarchyFeedback(feedback);
  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<Plenum_Rank>, std::move(context), fub::GodunovSplitting{});

  fub::amrex::MultiBlockSourceTerm<Kinetics> source_term(kinetics);

  fub::SplitSystemSourceLevelIntegrator level_integrator(
      std::move(system_solver), std::move(source_term),
      fub::StrangSplittingLumped{});

  fub::amrex::DiffusionSourceTermOptions diff_opts =
      fub::GetOptions(vm, "DiffusionSourceTerm");
  std::vector<fub::amrex::DiffusionSourceTerm<fub::PerfectGasMix<1>>> diffs(
      kinetics.size(), fub::amrex::DiffusionSourceTerm<fub::PerfectGasMix<1>>{
                           tube_equation, diff_opts});
  fub::amrex::MultiBlockSourceTerm<
      fub::amrex::DiffusionSourceTerm<fub::PerfectGasMix<1>>>
      diff_term(diffs);
  fub::SplitSystemSourceLevelIntegrator diff_integrator(
      std::move(level_integrator), std::move(diff_term),
      fub::StrangSplittingLumped{});

  fub::SubcycleFineFirstSolver solver(std::move(diff_integrator));

  using namespace fub::amrex;
  struct MakeCheckpoint
      : public fub::OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm2> {
    std::string directory_ = "SEC_Plenum_Arrhenius";
    std::shared_ptr<const fub::perfect_gas_mix::gt::ControlState>
        control_state_;
    MakeCheckpoint(
        const fub::ProgramOptions& options,
        std::shared_ptr<const fub::perfect_gas_mix::gt::ControlState> cs)
        : OutputAtFrequencyOrInterval(options), control_state_{cs} {
      directory_ = fub::GetOptionOr(options, "directory", directory_);
    }
    void operator()(const MultiBlockGriddingAlgorithm2& grid) override {
      std::string name = fmt::format("{}/{:09}", directory_, grid.GetCycles());
      fub::SeverityLogger log = fub::GetInfoLogger();
      BOOST_LOG(log) << fmt::format("Write Checkpoint to '{}'!", name);
      WriteCheckpoint(name, grid, control_state_);
    }
  };

  fub::OutputFactory<MultiBlockGriddingAlgorithm2> factory{};
  factory.RegisterOutput<MakeCheckpoint>("Checkpoint", control_state);
  factory.RegisterOutput<GT::ControlOutput>("ControlOutput", control_state);
  using CounterOutput =
      fub::CounterOutput<fub::amrex::MultiBlockGriddingAlgorithm2,
                         std::chrono::nanoseconds>;
  factory.RegisterOutput<fub::amrex::MultiWriteHdf5WithNames>(
      "HDF5", plenum_equation, tube_equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<
      MultiBlockPlotfileOutput2<fub::PerfectGasMix<1>, fub::PerfectGasMix<2>>>(
      "Plotfiles", tube_equation, plenum_equation);
  fub::MultipleOutputs<MultiBlockGriddingAlgorithm2> outputs(
      std::move(factory), fub::GetOptions(vm, "Output"));

  outputs(*solver.GetGriddingAlgorithm());
  // fub::SeverityLogger log = fub::GetInfoLogger();
  fub::RunOptions run_options = fub::GetOptions(vm, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);
  fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
}
