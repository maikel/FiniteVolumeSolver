// Copyright (c) 2020 Maikel Nadolski
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
#include "fub/flux_method/MusclHancockMethod2.hpp"

#include "fub/AMReX/AxialFluxMethodAdapter.hpp"
#include "fub/AMReX/AxialTimeIntegrator.hpp"

#include "fub/AMReX/multi_block/MultiBlockBoundary2.hpp"
#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm2.hpp"
#include "fub/AMReX/multi_block/MultiBlockIntegratorContext2.hpp"

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

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 2;

static_assert(AMREX_SPACEDIM == 2);

static constexpr double r_tube = 0.015;

struct InitialDataInTube {
  using Equation = fub::PerfectGasMix<1>;
  using Complete = fub::Complete<Equation>;
  using KineticState = fub::KineticState<Equation>;
  using Conservative = fub::Conservative<Equation>;

  Equation equation_;
  double x_0_;
  double initially_filled_x_{0.4};

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
            // const double rel_x = x - x_0_;
            const double pressure = 1.0;
            state.temperature = 1.0; // (rel_x < 0.5) ? 1.0 : 2.5;
            state.density = pressure / state.temperature / equation_.Rspec;
            state.mole_fractions[0] = 0.0;
            state.mole_fractions[1] = 1.0;
            fub::euler::CompleteFromKineticState(equation_, complete, state,
                                                 velocity);
            fub::Store(states, complete, {i});
          });
        });
  }
};

template <typename Limiter, typename GradientMethod> struct RebindLimiterTo_;

template <typename Limiter, typename GradientMethod>
using RebindLimiterTo =
    typename RebindLimiterTo_<Limiter, GradientMethod>::type;

template <typename Limiter, typename OtherLimiter>
struct RebindLimiterTo_<Limiter, fub::CentralDifferenceGradient<OtherLimiter>> {
  using type = fub::CentralDifferenceGradient<Limiter>;
};

template <typename Limiter, typename Equation, typename GradientMethod>
struct RebindLimiterTo_<Limiter,
                        fub::ConservativeGradient<Equation, GradientMethod>> {
  using type =
      fub::ConservativeGradient<Equation,
                                RebindLimiterTo<Limiter, GradientMethod>>;
};

template <typename Limiter, typename Equation, typename GradientMethod>
struct RebindLimiterTo_<Limiter,
                        fub::PrimitiveGradient<Equation, GradientMethod>> {
  using type = fub::PrimitiveGradient<Equation,
                                      RebindLimiterTo<Limiter, GradientMethod>>;
};

template <typename Limiter, typename Equation, typename GradientMethod>
struct RebindLimiterTo_<
    Limiter, fub::CharacteristicsGradient<Equation, GradientMethod>> {
  using type =
      fub::CharacteristicsGradient<Equation,
                                   RebindLimiterTo<Limiter, GradientMethod>>;
};

template <typename Limiter, typename Equation, typename GradientMethod,
          typename Reconstruction, typename BaseMethod>
struct RebindLimiterTo_<
    Limiter,
    fub::MusclHancock2<Equation, GradientMethod, Reconstruction, BaseMethod>> {
  using type =
      fub::MusclHancock2<Equation, RebindLimiterTo<Limiter, GradientMethod>,
                         Reconstruction, BaseMethod>;
};

template <typename Limiter, typename BaseMethod>
struct RebindLimiterTo_<Limiter, fub::FluxMethod<BaseMethod>> {
  using type = fub::FluxMethod<RebindLimiterTo<Limiter, BaseMethod>>;
};

template <typename BaseMethod, typename FluxMethod> struct RebindBaseMethod_;

template <typename BaseMethod, typename FluxMethod>
using RebindBaseMethod =
    typename RebindBaseMethod_<BaseMethod, FluxMethod>::type;

template <typename BaseMethod, typename FM>
struct RebindBaseMethod_<BaseMethod, fub::FluxMethod<FM>> {
  using type = fub::FluxMethod<RebindBaseMethod<BaseMethod, FM>>;
};

template <typename BaseMethod, typename Equation, typename GradientMethod,
          typename Reconstruction, typename OtherBaseMethod>
struct RebindBaseMethod_<BaseMethod,
                         fub::MusclHancock2<Equation, GradientMethod,
                                            Reconstruction, OtherBaseMethod>> {
  using type =
      fub::MusclHancock2<Equation, GradientMethod, Reconstruction, BaseMethod>;
};

template <typename Limiter, typename FluxMethod, typename BaseMethod,
          typename Equation>
auto RebindLimiter(fub::Type<FluxMethod>, fub::Type<BaseMethod>,
                   const Equation& equation) {
  return RebindBaseMethod<BaseMethod, RebindLimiterTo<Limiter, FluxMethod>>(
      equation);
}

fub::AnyFluxMethod<fub::amrex::cutcell::IntegratorContext>
GetCutCellMethod(const fub::ProgramOptions& options,
                 const fub::PerfectGasMix<2>& equation) {
  using Limiter = std::variant<fub::NoLimiter2, fub::UpwindLimiter,
                               fub::MinModLimiter, fub::VanLeerLimiter>;
  using namespace std::literals;

  const std::map<std::string, Limiter> limiters{
      std::pair{"NoLimiter"s, Limiter{fub::NoLimiter2{}}},
      std::pair{"Upwind"s, Limiter{fub::UpwindLimiter{}}},
      std::pair{"MinMod"s, Limiter{fub::MinModLimiter{}}},
      std::pair{"VanLeer"s, Limiter{fub::VanLeerLimiter{}}}};

  using HLLEM = fub::perfect_gas::HllemMethod<fub::PerfectGasMix<2>, false>;
  using HLLEM_Lar = fub::perfect_gas::HllemMethod<fub::PerfectGasMix<2>>;

  using BaseMethod = std::variant<fub::Type<HLLEM>, fub::Type<HLLEM_Lar>>;

  const std::map<std::string, BaseMethod> base_methods{
      std::pair{"HLLEM"s, BaseMethod{fub::Type<HLLEM>{}}},
      std::pair{"HLLEM_Larrouturou"s, BaseMethod{fub::Type<HLLEM_Lar>{}}}};

  using ConservativeReconstruction = fub::FluxMethod<fub::MusclHancock2<
      fub::PerfectGasMix<2>,
      fub::ConservativeGradient<
          fub::PerfectGasMix<2>,
          fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
      fub::ConservativeReconstruction<fub::PerfectGasMix<2>>, HLLEM>>;

  using PrimitiveReconstruction = fub::FluxMethod<fub::MusclHancock2<
      fub::PerfectGasMix<2>,
      fub::PrimitiveGradient<
          fub::PerfectGasMix<2>,
          fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
      fub::PrimitiveReconstruction<fub::PerfectGasMix<2>>, HLLEM>>;

  using CharacteristicsReconstruction = fub::FluxMethod<fub::MusclHancock2<
      fub::PerfectGasMix<2>,
      fub::CharacteristicsGradient<
          fub::PerfectGasMix<2>,
          fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
      fub::CharacteristicsReconstruction<fub::PerfectGasMix<2>>, HLLEM>>;

  using Reconstruction = std::variant<fub::Type<ConservativeReconstruction>,
                                      fub::Type<PrimitiveReconstruction>,
                                      fub::Type<CharacteristicsReconstruction>>;

  const std::map<std::string, Reconstruction> reconstructions{
      std::pair{"Conservative"s,
                Reconstruction{fub::Type<ConservativeReconstruction>{}}},
      std::pair{"Primitive"s,
                Reconstruction{fub::Type<PrimitiveReconstruction>{}}},
      std::pair{"Characteristics"s,
                Reconstruction{fub::Type<CharacteristicsReconstruction>{}}}};

  std::string limiter_option = fub::GetOptionOr(options, "limiter", "MinMod"s);
  std::string reconstruction_option =
      fub::GetOptionOr(options, "reconstruction", "Characteristics"s);
  std::string base_method_option =
      fub::GetOptionOr(options, "base_method", "HLLEM_Larrouturou"s);

  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << "FluxMethod:";
  BOOST_LOG(log) << " - limiter = " << limiter_option;
  BOOST_LOG(log) << " - reconstruction = " << reconstruction_option;
  BOOST_LOG(log) << " - base_method = " << base_method_option;

  Limiter limiter = limiters.at(limiter_option);
  Reconstruction reconstruction = reconstructions.at(reconstruction_option);
  BaseMethod base_method = base_methods.at(base_method_option);

  return std::visit(
      [&equation](auto limiter, auto reconstruction, auto base_method_type) {
        using ThisLimiter = fub::remove_cvref_t<decltype(limiter)>;
        auto flux_method = RebindLimiter<ThisLimiter>(
            reconstruction, base_method_type, equation);
        const auto base_method = flux_method.GetBaseMethod();
        fub::KbnCutCellMethod cutcell_method(flux_method, base_method);
        fub::amrex::cutcell::FluxMethod adapter(std::move(cutcell_method));
        fub::AnyFluxMethod<fub::amrex::cutcell::IntegratorContext> any(adapter);
        return any;
      },
      limiter, reconstruction, base_method);
}

std::pair<fub::AnyFluxMethod<fub::amrex::IntegratorContext>,
          fub::AnyTimeIntegrator<fub::amrex::IntegratorContext>>
GetFluxMethod(const fub::ProgramOptions& options,
              const fub::amrex::PatchHierarchy& hier,
              const fub::PerfectGasMix<1>& equation) {
  using Limiter = std::variant<fub::NoLimiter2, fub::UpwindLimiter,
                               fub::MinModLimiter, fub::VanLeerLimiter>;
  using namespace std::literals;

  const std::map<std::string, Limiter> limiters{
      std::pair{"NoLimiter"s, Limiter{fub::NoLimiter2{}}},
      std::pair{"Upwind"s, Limiter{fub::UpwindLimiter{}}},
      std::pair{"MinMod"s, Limiter{fub::MinModLimiter{}}},
      std::pair{"VanLeer"s, Limiter{fub::VanLeerLimiter{}}}};

  using HLLEM = fub::perfect_gas::HllemMethod<fub::PerfectGasMix<1>, false>;
  using HLLEM_Lar = fub::perfect_gas::HllemMethod<fub::PerfectGasMix<1>>;

  using BaseMethod = std::variant<fub::Type<HLLEM>, fub::Type<HLLEM_Lar>>;

  const std::map<std::string, BaseMethod> base_methods{
      std::pair{"HLLEM"s, BaseMethod{fub::Type<HLLEM>{}}},
      std::pair{"HLLEM_Larrouturou"s, BaseMethod{fub::Type<HLLEM_Lar>{}}}};

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

  using Reconstruction = std::variant<fub::Type<ConservativeReconstruction>,
                                      fub::Type<PrimitiveReconstruction>,
                                      fub::Type<CharacteristicsReconstruction>>;

  const std::map<std::string, Reconstruction> reconstructions{
      std::pair{"Conservative"s,
                Reconstruction{fub::Type<ConservativeReconstruction>{}}},
      std::pair{"Primitive"s,
                Reconstruction{fub::Type<PrimitiveReconstruction>{}}},
      std::pair{"Characteristics"s,
                Reconstruction{fub::Type<CharacteristicsReconstruction>{}}}};

  std::string limiter_option = fub::GetOptionOr(options, "limiter", "MinMod"s);
  std::string reconstruction_option =
      fub::GetOptionOr(options, "reconstruction", "Characteristics"s);
  std::string base_method_option =
      fub::GetOptionOr(options, "base_method", "HLLEM_Larrouturou"s);

  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << "FluxMethod:";
  BOOST_LOG(log) << " - limiter = " << limiter_option;
  BOOST_LOG(log) << " - reconstruction = " << reconstruction_option;
  BOOST_LOG(log) << " - base_method = " << base_method_option;

  auto area = [&]() -> std::optional<pybind11::function> {
    if (auto iter = options.find("area_variation"); iter != options.end()) {
      pybind11::object obj = iter->second;
      pybind11::function fun(obj);
      return fun;
    }
    return {};
  }();

  auto area_lambda = [area](double x) -> double {
    if (area) {
      pybind11::function py_area = *area;
      pybind11::object py_y = py_area(x);
      double y = py_y.cast<double>();
      return y;
    }
    return 1.0;
  };

  Limiter limiter = limiters.at(limiter_option);
  Reconstruction reconstruction = reconstructions.at(reconstruction_option);
  BaseMethod base_method = base_methods.at(base_method_option);

  return std::visit(
      [&equation, &hier, &area_lambda](auto limiter, auto reconstruction,
                                       auto base_method_type) {
        using ThisLimiter = fub::remove_cvref_t<decltype(limiter)>;
        auto flux_method = RebindLimiter<ThisLimiter>(
            reconstruction, base_method_type, equation);
        fub::amrex::AxialFluxMethodAdapter adapted(flux_method);
        auto pressure = adapted.SharedInterfacePressure();
        const ::amrex::Vector<::amrex::Geometry> geoms = hier.GetGeometries();
        std::vector<::amrex::Geometry> std_geoms(geoms.begin(), geoms.end());
        fub::amrex::AxialTimeIntegrator time_integrator(equation, std_geoms,
                                                        pressure, area_lambda);
        fub::AnyFluxMethod<fub::amrex::IntegratorContext> any_flux(
            std::move(adapted));
        fub::AnyTimeIntegrator<fub::amrex::IntegratorContext> any_time(
            std::move(time_integrator));
        return std::pair{any_flux, any_time};
      },
      limiter, reconstruction, base_method);
}

auto MakeTubeSolver(const fub::ProgramOptions& options,
                    const std::shared_ptr<fub::CounterRegistry>& counters) {
  using namespace fub::amrex;

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
  fub::PerfectGasMix<Tube_Rank> equation{};
  equation.n_species = 1;
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 1e-2),
                            std::make_pair(&Complete::pressure, 1e-2)};

  DataDescription desc = MakeDataDescription(equation);

  ::amrex::Box refine_box{{grid_geometry.cell_dimensions[0] - 5, 0},
                          {grid_geometry.cell_dimensions[0] - 1, 0}};
  ConstantBox constant_box{refine_box};

  const double initially_filled_x =
      fub::GetOptionOr(options, "initially_filled_x", 0.4);
  BOOST_LOG(log) << "InitialData:";
  BOOST_LOG(log) << "  - initially_filled_x = " << initially_filled_x << " [m]";
  InitialDataInTube initial_data{equation, grid_geometry.coordinates.lo()[0],
                                 initially_filled_x};

  fub::perfect_gas_mix::ArrheniusKinetics<1> source_term{
      equation, fub::GetOptions(options, "ArrheniusKinetics")};
  BOOST_LOG(log) << "ArrheniusKinetics:";
  source_term.options.Print(log);

  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());

  auto ignition_delay = [=](double Y, double T) {
    const double gm1 = equation.gamma - 1.0;
    double TT = T / (gm1 * source_term.options.Q * std::max(Y, eps));
    double Ea0 = source_term.options.EA * (1.0 / T);
    double B0 =
        source_term.options.B * exp(source_term.options.EA * (1.0 - (1.0 / T)));
    return ((TT / (B0 * Ea0)) * (1.0 + (2.0 + TT) / Ea0));
  };

  // /*
  // ==================================================================================
  // */

  auto d_ignition_delay_dT = [=](double Y, double T) {
    const double gm1 = equation.gamma - 1.0;
    double TQ = T / (gm1 * source_term.options.Q * std::max(Y, eps));
    double TE = T / source_term.options.EA;
    double A = 1.0 / (source_term.options.B * source_term.options.EA * gm1 *
                      source_term.options.Q * Y);
    double dtidT = A * source_term.options.EA *
                   ((-1.0 + 2.0 * TE) * (1.0 + TE * (2.0 + TQ)) +
                    2.0 * TE * TE * (1.0 + TQ)) *
                   exp(-source_term.options.EA * (1.0 - 1.0 / T));
    return dtidT;
  };

  // /*
  // ==================================================================================
  // */

  auto temperature_for_ignition_delay = [=](double tign, double Y, double T0) {
    double T = T0;
    double ti = ignition_delay(Y, T);
    double dti = ti - tign;
    while (std::abs(dti) > eps) {
      T -= dti / d_ignition_delay_dT(Y, T);
      ti = ignition_delay(Y, T);
      dti = ti - tign;
    }
    return T;
  };

  auto inflow_function =
      [temperature_for_ignition_delay,
       prim = fub::Primitive<fub::PerfectGasMix<1>>(equation)](
          const fub::PerfectGasMix<1>& eq,
          fub::Complete<fub::PerfectGasMix<1>>& boundary_state,
          const fub::KineticState<fub::PerfectGasMix<1>>& compressor_state,
          double inner_pressure, fub::Duration t_diff, const amrex::MultiFab&,
          const fub::amrex::GriddingAlgorithm&, int) mutable {
        const double fuel_retardatation =
            0.06;                /* Reference: 0.06;   for icx = 256:  0.1*/
        const double tti = 0.75; /* Reference: 0.75; */
        const double timin = 0.1;
        const double X_inflow_left = 1.0;

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

        const double tign = std::max(timin, tti - t_diff.count());
        const double Xin = X_inflow_left;
        const double Tin1 = temperature_for_ignition_delay(tign, Xin, Tin);

        /* adjust density to match desired temperature */
        rhoin *= Tin / Tin1;

        prim.density = rhoin;
        prim.velocity[0] = uin;
        prim.pressure = pin;
        auto heaviside = [](double x) { return (x > 0); };
        prim.species[0] = std::clamp(
            Xin * heaviside(t_diff.count() - fuel_retardatation), 0.0, 1.0);

        fub::CompleteFromPrim(eq, boundary_state, prim);
      };

  fub::KineticState<fub::PerfectGasMix<1>> compressor_state(equation);
  compressor_state.temperature = 1.0;
  compressor_state.density = 1.05 / compressor_state.temperature;
  compressor_state.mole_fractions[1] = 1.0;

  fub::amrex::PressureValveBoundary_Klein valve(equation, compressor_state,
                                                inflow_function);

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
      checkpoint = fmt::format("{}/Tube", checkpoint);
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
      GetFluxMethod(fub::GetOptions(options, "FluxMethod"),
                    gridding->GetPatchHierarchy(), equation);
  HyperbolicMethod method{flux_method, time_integrator,
                          Reconstruction(equation)};

  const int scratch_gcw = 4;
  const int flux_gcw = 2;

  IntegratorContext context(gridding, method, scratch_gcw, flux_gcw);

  BOOST_LOG(log) << "==================== End Tube =========================";

  return std::pair{context, source_term};
}

auto MakePlenumSolver(const std::map<std::string, pybind11::object>& options) {
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
  std::shared_ptr<fub::CounterRegistry> registry =
      std::make_shared<fub::CounterRegistry>();
  {
    fub::Timer timer = registry->get_timer("MakeIndexSpaces");
    hierarchy_options.index_spaces =
        MakeIndexSpaces(shop, grid_geometry, hierarchy_options);
  }

  constexpr static int n_species = 2;
  fub::PerfectGasMix<Plenum_Rank> equation{};
  equation.n_species = n_species - 1;

  fub::Complete<fub::PerfectGasMix<Plenum_Rank>> left(equation);
  fub::Complete<fub::PerfectGasMix<Plenum_Rank>> right(equation);
  {
    using namespace std::literals;
    const fub::ProgramOptions initial_options =
        fub::GetOptions(options, "InitialCondition");
    const fub::ProgramOptions left_options =
        fub::GetOptions(initial_options, "left");
    std::array<double, n_species> moles{0.0, 1.0};
    moles = fub::GetOptionOr(left_options, "moles", moles);
    double temperature = fub::GetOptionOr(left_options, "temperature", 1.0);
    double density = fub::GetOptionOr(left_options, "density", 1.0);
    std::array<double, Plenum_Rank> velocity{};
    velocity = fub::GetOptionOr(left_options, "velocity", velocity);
    fub::KineticState<fub::PerfectGasMix<Plenum_Rank>> kin(equation);
    kin.temperature = temperature;
    kin.density = density;
    kin.mole_fractions = fub::AsEigenVector(moles);
    fub::euler::CompleteFromKineticState(equation, left, kin,
                                         fub::AsEigenVector(velocity));

    const fub::ProgramOptions right_options =
        fub::GetOptions(initial_options, "right");
    moles = std::array<double, n_species>{0.0, 1.0};
    moles = fub::GetOptionOr(right_options, "moles", moles);
    temperature = fub::GetOptionOr(right_options, "temperature", 1.0);
    density = fub::GetOptionOr(right_options, "density", 1.0);
    velocity = std::array<double, Plenum_Rank>{};
    velocity = fub::GetOptionOr(right_options, "velocity", velocity);
    kin.temperature = temperature;
    kin.density = density;
    kin.mole_fractions = fub::AsEigenVector(moles);
    fub::euler::CompleteFromKineticState(equation, right, kin,
                                         fub::AsEigenVector(velocity));
  }
  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.0),
                              left, right);

  //  using Complete = fub::Complete<fub::PerfectGas<Plenum_Rank>>;
  using Complete = fub::Complete<fub::PerfectGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&Complete::pressure, 0.05},
                             std::pair{&Complete::density, 0.05}};

  const int scratch_gcw = 4;
  const int flux_gcw = 2;

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
      {TransmissiveBoundary{fub::Direction::X, 0},
       ReflectiveBoundary2{equation, fub::Direction::X, 1, lower_right_corner},
       ReflectiveBoundary2{equation, fub::Direction::X, 1, upper_right_corner},
       ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::X, 1},
       ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::Y, 0},
       ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::Y,
                          1}}};

  std::vector<pybind11::dict> dicts{};
  dicts = fub::GetOptionOr(options, "PressureOutflowBoundaries", dicts);
  for (pybind11::dict& dict : dicts) {
    fub::ProgramOptions boundary_options = fub::ToMap(dict);
    fub::amrex::IsentropicPressureBoundaryOptions pb_opts(boundary_options);
    pb_opts.direction = fub::Direction::X;
    pb_opts.side = 1;
    BOOST_LOG(log) << "PressureOutflowBoundary:";
    pb_opts.Print(log);
    fub::amrex::cutcell::IsentropicPressureExpansion<fub::PerfectGasMix<2>>
        pressure_outflow(equation, pb_opts);
    boundary_condition.conditions.push_back(std::move(pressure_outflow));
  }

  dicts.clear();
  dicts = fub::GetOptionOr(options, "TurbineMassflowBoundaries", dicts);
  for (pybind11::dict& dict : dicts) {
    fub::ProgramOptions boundary_options = fub::ToMap(dict);
    fub::amrex::cutcell::TurbineMassflowBoundaryOptions tb_opts(
        boundary_options);
    tb_opts.dir = fub::Direction::X;
    tb_opts.side = 1;
    BOOST_LOG(log) << "TurbineMassflowBoundaries:";
    tb_opts.Print(log);
    fub::amrex::cutcell::TurbineMassflowBoundary<fub::PerfectGasMix<2>>
        pressure_outflow(equation, tb_opts);
    boundary_condition.conditions.push_back(std::move(pressure_outflow));
  }

  dicts.clear();
  dicts = fub::GetOptionOr(options, "TurbineMassflowBoundaries_Jirasek", dicts);
  for (pybind11::dict& dict : dicts) {
    fub::ProgramOptions boundary_options = fub::ToMap(dict);
    fub::amrex::cutcell::TurbineMassflowBoundaryOptions tb_opts(
        boundary_options);
    tb_opts.dir = fub::Direction::X;
    tb_opts.side = 1;
    BOOST_LOG(log) << "TurbineMassflowBoundaries_Jirasek:";
    tb_opts.Print(log);
    fub::amrex::cutcell::TurbineMassflowBoundary<fub::PerfectGasMix<2>,
                                                 fub::RequireMassflow_Jirasek>
        pressure_outflow(equation, tb_opts);
    boundary_condition.conditions.push_back(std::move(pressure_outflow));
  }

  dicts.clear();
  dicts = fub::GetOptionOr(options, "MachnumberBoundaries", dicts);
  for (pybind11::dict& dict : dicts) {
    fub::ProgramOptions boundary_options = fub::ToMap(dict);
    fub::amrex::cutcell::MachnumberBoundaryOptions mb_opts(boundary_options);
    mb_opts.dir = fub::Direction::X;
    mb_opts.side = 1;
    BOOST_LOG(log) << "MachnumberBoundary:";
    mb_opts.Print(log);
    fub::amrex::cutcell::MachnumberBoundary<fub::PerfectGasMix<2>>
        mach_boundary(equation, mb_opts);
    boundary_condition.conditions.push_back(std::move(mach_boundary));
  }

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
      grid->InitializeHierarchy(0.0);
      return grid;
    } else {
      checkpoint += "/Plenum";
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), grid_geometry,
          hierarchy_options);
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), constant_refinebox), boundary_condition);
      return grid;
    }
  }();

  auto flux_method =
      GetCutCellMethod(fub::GetOptions(options, "FluxMethod"), equation);

  HyperbolicMethod method{flux_method, TimeIntegrator{},
                          Reconstruction{equation}};

  IntegratorContext context(gridding, method, scratch_gcw, flux_gcw);

  auto log_massflow = [](IntegratorContext& plenum, int level, fub::Duration,
                         std::pair<int, int>) {
    fub::SeverityLogger log =
        fub::GetLogger(boost::log::trivial::severity_level::debug);
    const fub::Duration time_point = plenum.GetTimePoint(level);
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "LogMassflow");
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", time_point.count());
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Level", level);
    const amrex::MultiFab& fluxes_x =
        plenum.GetFluxes(level, fub::Direction::X);
    const amrex::Geometry& geom = plenum.GetGeometry(level);
    const amrex::Box cells = geom.Domain();
    const amrex::Box faces_x = amrex::convert(cells, {1, 0});
    const int boundary_n = faces_x.bigEnd(0);
    const amrex::IntVect smallEnd{boundary_n, faces_x.smallEnd(1)};
    const amrex::IntVect bigEnd = faces_x.bigEnd();
    const amrex::Box right_boundary{smallEnd, bigEnd};
    // const double dy = geom.CellSize(1);
    double local_f_rho = 0.0;
    fub::amrex::ForEachFab(fluxes_x, [&](const amrex::MFIter& mfi) {
      amrex::Box local_boundary = mfi.tilebox() & right_boundary;
      if (local_boundary.ok()) {
        const amrex::FArrayBox& local_fluxes = fluxes_x[mfi];
        local_f_rho += local_fluxes.sum(local_boundary, 0);
      }
    });
    double global_f_rho = 0.0;
    ::MPI_Allreduce(&local_f_rho, &global_f_rho, 1, MPI_DOUBLE, MPI_SUM,
                    ::amrex::ParallelDescriptor::Communicator());
    global_f_rho /= right_boundary.numPts();
    BOOST_LOG(log) << fmt::format("Average F_rho = {:.12g}", global_f_rho)
                   << boost::log::add_value("average_massflow", global_f_rho);
  };
  context.SetFeedbackFunction(log_massflow);

  BOOST_LOG(log) << "==================== End Plenum =========================";

  return context;
}

void MyMain(const std::map<std::string, pybind11::object>& vm);

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

void WriteCheckpoint(const std::string& path,
                     const fub::amrex::MultiBlockGriddingAlgorithm2& grid) {
  auto tubes = grid.GetTubes();
  std::string name = fmt::format("{}/Tube", path);
  fub::amrex::WriteCheckpointFile(name, tubes[0]->GetPatchHierarchy());
  name = fmt::format("{}/Plenum", path);
  fub::amrex::cutcell::WriteCheckpointFile(
      name, grid.GetPlena()[0]->GetPatchHierarchy());
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

  plenum.push_back(MakePlenumSolver(fub::GetOptions(vm, "Plenum")));
  auto counter_database = plenum[0].GetCounterRegistry();

  std::vector<pybind11::dict> tube_dicts = {};
  tube_dicts = fub::GetOptionOr(vm, "Tubes", tube_dicts);
  for (pybind11::dict& dict : tube_dicts) {
    fub::ProgramOptions tube_options = fub::ToMap(dict);
    auto&& [tube, source] = MakeTubeSolver(tube_options, counter_database);
    fub::amrex::BlockConnection connection;
    connection.direction = fub::Direction::X;
    connection.side = 0;
    connection.ghost_cell_width = 4;
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

  fub::PerfectGasMix<Plenum_Rank> plenum_equation{};
  plenum_equation.n_species = 1;
  fub::PerfectGasMix<Tube_Rank> tube_equation{};
  tube_equation.n_species = 1;

  fub::amrex::MultiBlockIntegratorContext2 context(
      tube_equation, plenum_equation, std::move(tubes), std::move(plenum),
      std::move(connectivity));

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<Plenum_Rank>, std::move(context), fub::GodunovSplitting{});

  fub::amrex::MultiBlockSourceTerm<Kinetics> source_term(kinetics);

  fub::SplitSystemSourceLevelIntegrator level_integrator(
      std::move(system_solver), std::move(source_term),
      fub::GodunovSplitting{});

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  using namespace fub::amrex;
  struct MakeCheckpoint
      : public fub::OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm2> {
    std::string directory_ = "ConvergentNozzle";
    MakeCheckpoint(const fub::ProgramOptions& options)
        : OutputAtFrequencyOrInterval(options) {
      directory_ = fub::GetOptionOr(options, "directory", directory_);
    }
    void operator()(const MultiBlockGriddingAlgorithm2& grid) override {
      std::string name = fmt::format("{}/{:09}", directory_, grid.GetCycles());
      fub::SeverityLogger log = fub::GetInfoLogger();
      BOOST_LOG(log) << fmt::format("Write Checkpoint to '{}'!", name);
      WriteCheckpoint(name, grid);
    }
  };

  fub::OutputFactory<MultiBlockGriddingAlgorithm2> factory{};
  factory.RegisterOutput<MakeCheckpoint>("Checkpoint");
  using CounterOutput =
      fub::CounterOutput<fub::amrex::MultiBlockGriddingAlgorithm2,
                         std::chrono::milliseconds>;
  factory.RegisterOutput<fub::amrex::MultiWriteHdf52>("HDF5");
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<
      MultiBlockPlotfileOutput2<fub::PerfectGasMix<1>, fub::PerfectGasMix<2>>>(
      "Plotfiles", tube_equation, plenum_equation);
  factory.RegisterOutput<
      MultiBlockPlotfileOutput2<fub::PerfectGasMix<1>, fub::PerfectGasMix<2>>>(
      "Plotfiles", tube_equation, plenum_equation);
  fub::MultipleOutputs<MultiBlockGriddingAlgorithm2> outputs(
      std::move(factory), fub::GetOptions(vm, "Output"));

  outputs(*solver.GetGriddingAlgorithm());
  fub::SeverityLogger log = fub::GetInfoLogger();
  fub::RunOptions run_options = fub::GetOptions(vm, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);
  fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
}
