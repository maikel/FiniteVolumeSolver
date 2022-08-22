// Copyright (c) 2021 Maikel Nadolski
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

#ifndef FUB_AMREX_FLUXMETHOD_FACTORY
#define FUB_AMREX_FLUXMETHOD_FACTORY

#include "fub/flux_method/Gradient.hpp"
#include "fub/flux_method/MusclHancockMethod2.hpp"
#include "fub/flux_method/Reconstruct.hpp"

#include "fub/ext/Log.hpp"

#include "fub/equations/ideal_gas_mix/EinfeldtSignalVelocities.hpp"
#include "fub/equations/perfect_gas/HllemMethod.hpp"

#include "fub/EinfeldtSignalVelocities.hpp"
#include "fub/flux_method/HllMethod.hpp"

#include "fub/AMReX/AxialFluxMethodAdapter.hpp"
#include "fub/AMReX/AxialTimeIntegrator.hpp"
#include "fub/AMReX/FluxMethodAdapter.hpp"
#include "fub/AMReX/TimeIntegrator.hpp"
#include "fub/AMReX/IntegratorContext.hpp"

#include "fub/equations/PerfectGasMix.hpp"

#include <map>
#include <optional>
#include <string>
#include <variant>

namespace fub {
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
} // namespace fub

namespace fub::amrex {
template <typename Equation>
std::pair<fub::AnyFluxMethod<fub::amrex::IntegratorContext>,
          fub::AnyTimeIntegrator<fub::amrex::IntegratorContext>>
GetFluxMethod(const fub::ProgramOptions& options,
              const fub::amrex::PatchHierarchy& hier,
              const Equation& equation) {
  if constexpr (Equation::Rank() == 1) {
    using Limiter = std::variant<fub::NoLimiter2, fub::UpwindLimiter,
                                 fub::MinModLimiter, fub::VanLeerLimiter>;
    // using Limiter = std::variant<fub::VanLeerLimiter>;
    using namespace std::literals;

    const std::map<std::string, Limiter> limiters{
        std::pair{"NoLimiter"s, Limiter{fub::NoLimiter2{}}},
        std::pair{"Upwind"s, Limiter{fub::UpwindLimiter{}}},
        std::pair{"MinMod"s, Limiter{fub::MinModLimiter{}}},
        std::pair{"VanLeer"s, Limiter{fub::VanLeerLimiter{}}}};

    using HLLE =
        fub::HllMethod<Equation, fub::EinfeldtSignalVelocities<Equation>>;
    using HLLEM = fub::perfect_gas::HllemMethod<Equation, false>;
    using HLLEM_Lar = fub::perfect_gas::HllemMethod<Equation>;

    using BaseMethod =
        std::variant<fub::Type<HLLE>, fub::Type<HLLEM>, fub::Type<HLLEM_Lar>>;
    // std::variant<fub::Type<HLLEM_Lar>>;

    const std::map<std::string, BaseMethod> base_methods{
        std::pair{"HLLE"s, BaseMethod{fub::Type<HLLE>{}}},
        std::pair{"HLLEM"s, BaseMethod{fub::Type<HLLEM>{}}},
        std::pair{"HLLEM_Larrouturou"s, BaseMethod{fub::Type<HLLEM_Lar>{}}}};

    using ConservativeReconstruction = fub::FluxMethod<fub::MusclHancock2<
        Equation,
        fub::ConservativeGradient<
            Equation, fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
        fub::ConservativeReconstruction<Equation>, HLLEM>>;

    using PrimitiveReconstruction = fub::FluxMethod<fub::MusclHancock2<
        Equation,
        fub::PrimitiveGradient<
            Equation, fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
        fub::PrimitiveReconstruction<Equation>, HLLEM>>;

    using CharacteristicsReconstruction = fub::FluxMethod<fub::MusclHancock2<
        Equation,
        fub::CharacteristicsGradient<
            Equation, fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
        fub::CharacteristicsReconstruction<Equation>, HLLEM>>;

    using Reconstruction =
        std::variant<fub::Type<ConservativeReconstruction>,
                     fub::Type<PrimitiveReconstruction>,
                     fub::Type<CharacteristicsReconstruction>>;

    const std::map<std::string, Reconstruction> reconstructions{
        std::pair{"Conservative"s,
                  Reconstruction{fub::Type<ConservativeReconstruction>{}}},
        std::pair{"Primitive"s,
                  Reconstruction{fub::Type<PrimitiveReconstruction>{}}},
        std::pair{"Characteristics"s,
                  Reconstruction{fub::Type<CharacteristicsReconstruction>{}}}};

    std::string limiter_option =
        fub::GetOptionOr(options, "limiter", "MinMod"s);
    std::string reconstruction_option =
        fub::GetOptionOr(options, "reconstruction", "Characteristics"s);
    std::string base_method_option =
        fub::GetOptionOr(options, "base_method", "HLLEM_Larrouturou"s);

    fub::SeverityLogger log = fub::GetInfoLogger();
    BOOST_LOG(log) << "FluxMethod:";
    BOOST_LOG(log) << " - limiter = " << limiter_option;
    BOOST_LOG(log) << " - reconstruction = " << reconstruction_option;
    BOOST_LOG(log) << " - base_method = " << base_method_option;
    if (auto iter = options.find("area_variation"); iter != options.end()) {
      BOOST_LOG(log) << " - area_variation = <Function>";
    } else {
      BOOST_LOG(log) << " - area_variation = None";
    }

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
        [&equation, &hier, &area_lambda,
         &area](auto limiter, auto reconstruction, auto base_method_type) {
          using ThisLimiter = remove_cvref_t<decltype(limiter)>;
          auto flux_method = RebindLimiter<ThisLimiter>(
              reconstruction, base_method_type, equation);
          if (area) {
            AxialFluxMethodAdapter adapted(flux_method);
            auto pressure = adapted.SharedInterfacePressure();
            const ::amrex::Vector<::amrex::Geometry> geoms =
                hier.GetGeometries();
            std::vector<::amrex::Geometry> std_geoms(geoms.begin(),
                                                     geoms.end());
            AxialTimeIntegrator time_integrator(equation, std_geoms, pressure,
                                                area_lambda);
            AnyFluxMethod<IntegratorContext> any_flux(std::move(adapted));
            AnyTimeIntegrator<IntegratorContext> any_time(
                std::move(time_integrator));
            return std::pair{any_flux, any_time};
          } else {
            FluxMethodAdapter adapted(flux_method);
            AnyFluxMethod<IntegratorContext> any_flux(std::move(adapted));
            AnyTimeIntegrator<IntegratorContext> any_time(
                EulerForwardTimeIntegrator());
            return std::pair{any_flux, any_time};
          }
        },
        limiter, reconstruction, base_method);
  } else {
  }
}

extern template std::pair<fub::AnyFluxMethod<fub::amrex::IntegratorContext>,
                          fub::AnyTimeIntegrator<fub::amrex::IntegratorContext>>
GetFluxMethod<PerfectGasMix<1>>(const fub::ProgramOptions& options,
                                const fub::amrex::PatchHierarchy& hier,
                                const PerfectGasMix<1>& equation);
} // namespace fub::amrex

#endif