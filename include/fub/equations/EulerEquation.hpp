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

#ifndef FUB_EQUATIONS_EULER_EQUATION_HPP
#define FUB_EQUATIONS_EULER_EQUATION_HPP

#include "fub/core/type_traits.hpp"

namespace fub::euler {
template <int Dim>
double KineticEnergy(double density,
                     const Eigen::Array<double, Dim, 1>& momentum) noexcept {
  return 0.5 * momentum.matrix().squaredNorm() / density;
}

template <int Dim, int N, int O, int MR, int MC>
Array1d KineticEnergy(
    const Array1d& density,
    const Eigen::Array<double, Dim, N, O, MR, MC>& momentum) noexcept {
  Array1d squaredMomentum = momentum.matrix().colwise().squaredNorm();
  return 0.5 * squaredMomentum / density;
}

template <int Dim, int N, int O, int MR, int MC>
Array1d KineticEnergy(const Array1d& density,
                      const Eigen::Array<double, Dim, N, O, MR, MC>& momentum,
                      MaskArray mask) noexcept {
  mask = mask && (density > 0.0);
  Array1d squaredMomentum = momentum.matrix().colwise().squaredNorm();
  Array1d safe_density = density;
  safe_density = mask.select(density, 1.0);
  FUB_ASSERT((safe_density > 0.0).all());
  return 0.5 * squaredMomentum / safe_density;
}

inline constexpr struct MaSqFn {
  template <typename Equation>
  constexpr auto operator()(Equation&& eq) const noexcept {
    if constexpr (is_tag_invocable<MaSqFn, Equation>::value) {
      return fub::meta::tag_invoke(*this, std::forward<Equation>(eq));
    } else {
      return 1.0;
    }
  }
} MaSq;

inline constexpr struct GammaFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<GammaFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const
      noexcept(is_nothrow_tag_invocable<GammaFn, Equation, State>::value)
          -> tag_invoke_result_t<GammaFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state));
  }
} Gamma;

inline constexpr struct DensityFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<DensityFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const
      noexcept(is_nothrow_tag_invocable<DensityFn, Equation, State>::value)
          -> tag_invoke_result_t<DensityFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state));
  }
} Density;

inline constexpr struct MomentumFn {
  template <typename Equation, typename State, typename... Indices,
            typename = std::enable_if_t<is_tag_invocable<
                MomentumFn, Equation, State, Indices...>::value>>
  constexpr auto
  operator()(Equation&& eq, State&& state, Indices... d) const noexcept(
      is_nothrow_tag_invocable<MomentumFn, Equation, State, Indices...>::value)
      -> tag_invoke_result_t<MomentumFn, Equation, State, Indices...> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state), d...);
  }
} Momentum;

inline constexpr struct VelocityFn {
  template <typename Equation, typename State, typename... Indices,
            typename = std::enable_if_t<is_tag_invocable<
                VelocityFn, Equation, State, Indices...>::value>>
  constexpr auto
  operator()(Equation&& eq, State&& state, Indices... d) const noexcept(
      is_nothrow_tag_invocable<VelocityFn, Equation, State, Indices...>::value)
      -> tag_invoke_result_t<VelocityFn, Equation, State, Indices...> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state), d...);
  }
} Velocity;

inline constexpr struct SetVelocityFn {
  template <typename Equation, typename State, typename... Indices,
            typename = std::enable_if_t<is_tag_invocable<
                SetVelocityFn, Equation, State, Indices...>::value>>
  constexpr auto operator()(Equation&& eq, State&& state, Indices... d) const
      noexcept(is_nothrow_tag_invocable<SetVelocityFn, Equation, State,
                                        Indices...>::value)
          -> tag_invoke_result_t<SetVelocityFn, Equation, State, Indices...> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state), d...);
  }
} SetVelocity;

inline constexpr struct SpeedOfSoundFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<SpeedOfSoundFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const
      noexcept(is_nothrow_tag_invocable<SpeedOfSoundFn, Equation, State>::value)
          -> tag_invoke_result_t<SpeedOfSoundFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state));
  }
} SpeedOfSound;

inline constexpr struct EnergyFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<EnergyFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const
      noexcept(is_nothrow_tag_invocable<EnergyFn, Equation, State>::value)
          -> tag_invoke_result_t<EnergyFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state));
  }
} Energy;

inline constexpr struct PressureFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<PressureFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const
      noexcept(is_nothrow_tag_invocable<PressureFn, Equation, State>::value)
          -> tag_invoke_result_t<PressureFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state));
  }
} Pressure;

inline constexpr struct MachnumberFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<MachnumberFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const
      noexcept(is_nothrow_tag_invocable<MachnumberFn, Equation, State>::value)
          -> tag_invoke_result_t<MachnumberFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state));
  }
} Machnumber;

inline constexpr struct TemperatureFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<TemperatureFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const
      noexcept(is_nothrow_tag_invocable<TemperatureFn, Equation, State>::value)
          -> tag_invoke_result_t<TemperatureFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state));
  }
} Temperature;

inline constexpr struct InternalEnergyFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<InternalEnergyFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const noexcept(
      is_nothrow_tag_invocable<InternalEnergyFn, Equation, State>::value)
      -> tag_invoke_result_t<InternalEnergyFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state));
  }
} InternalEnergy;

inline constexpr struct TotalEnthalpyFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<TotalEnthalpyFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const noexcept(
      is_nothrow_tag_invocable<TotalEnthalpyFn, Equation, State>::value)
      -> tag_invoke_result_t<TotalEnthalpyFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state));
  }
} TotalEnthalpy;

inline constexpr struct MoleFractionsFn {
  template <typename Equation, typename Dest, typename State,
            typename = std::enable_if_t<is_tag_invocable<
                MoleFractionsFn, Equation, Dest, State>::value>>
  constexpr auto operator()(Equation&& eq, Dest&& dest, State&& state) const
      noexcept(
          is_nothrow_tag_invocable<MoleFractionsFn, Equation, State>::value)
          -> tag_invoke_result_t<MoleFractionsFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<Dest>(dest),
                                 std::forward<State>(state));
  }
} MoleFractions;

inline constexpr struct KineticStateFromCompleteFn {
  template <typename Equation, typename Dest, typename State,
            typename = std::enable_if_t<is_tag_invocable<
                KineticStateFromCompleteFn, Equation, Dest, State>::value>>
  constexpr auto operator()(Equation&& eq, Dest&& dest, State&& state) const
      noexcept(is_nothrow_tag_invocable<KineticStateFromCompleteFn, Equation,
                                        Dest, State>::value)
          -> tag_invoke_result_t<KineticStateFromCompleteFn, Equation, Dest,
                                 State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<Dest>(dest),
                                 std::forward<State>(state));
  }
} KineticStateFromComplete;

inline constexpr struct CompleteFromKineticStateFn {
  template <
      typename Equation, typename Dest, typename State, typename Velocity,
      typename = std::enable_if_t<is_tag_invocable<
          CompleteFromKineticStateFn, Equation, Dest, State, Velocity>::value>>
  constexpr auto operator()(Equation&& eq, Dest&& dest, State&& state,
                            Velocity&& velocity) const
      noexcept(is_nothrow_tag_invocable<CompleteFromKineticStateFn, Equation,
                                        Dest, State, Velocity>::value)
          -> tag_invoke_result_t<CompleteFromKineticStateFn, Equation, Dest,
                                 State, Velocity> {
    return fub::meta::tag_invoke(
        *this, std::forward<Equation>(eq), std::forward<Dest>(dest),
        std::forward<State>(state), std::forward<Velocity>(velocity));
  }

  template <typename Equation, typename Dest, typename State,
            typename = std::enable_if_t<is_tag_invocable<
                CompleteFromKineticStateFn, Equation, Dest, State>::value>>
  constexpr auto operator()(Equation&& eq, Dest&& dest, State&& state) const
      noexcept(is_nothrow_tag_invocable<CompleteFromKineticStateFn, Equation,
                                        Dest, State>::value)
          -> tag_invoke_result_t<CompleteFromKineticStateFn, Equation, Dest,
                                 State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<Dest>(dest),
                                 std::forward<State>(state));
  }
} CompleteFromKineticState;

inline constexpr struct SetIsentropicPressureFn {
  template <typename... Args, typename = std::enable_if_t<is_tag_invocable<
                                  SetIsentropicPressureFn, Args...>::value>>
  constexpr auto operator()(Args&&... args) const noexcept(
      is_nothrow_tag_invocable<SetIsentropicPressureFn, Args...>::value)
      -> tag_invoke_result_t<SetIsentropicPressureFn, Args...> {
    return fub::meta::tag_invoke(*this, std::forward<Args>(args)...);
  }
} SetIsentropicPressure;

inline constexpr struct IsentropicExpansionWithoutDissipationFn {
  template <typename Equation,
            typename = std::enable_if_t<
                is_tag_invocable<IsentropicExpansionWithoutDissipationFn,
                                 Equation, Complete<std::decay_t<Equation>>&,
                                 const Complete<std::decay_t<Equation>>&,
                                 double, double>::value ||
                is_tag_invocable<SetIsentropicPressureFn, Equation,
                                 Complete<std::decay_t<Equation>>&,
                                 const Complete<std::decay_t<Equation>>&,
                                 double>::value>>
  constexpr void operator()(Equation&& eq,
                            Complete<std::decay_t<Equation>>& dest,
                            const Complete<std::decay_t<Equation>>& src,
                            double pressure_dest, double efficiency) const
      noexcept(
          is_nothrow_tag_invocable<IsentropicExpansionWithoutDissipationFn,
                                   Equation, Complete<std::decay_t<Equation>>&,
                                   const Complete<std::decay_t<Equation>>&,
                                   double, double>::value) {
    if constexpr (is_tag_invocable<IsentropicExpansionWithoutDissipationFn,
                                   Equation, Complete<std::decay_t<Equation>>&,
                                   const Complete<std::decay_t<Equation>>&,
                                   double, double>::value) {
      fub::meta::tag_invoke(*this, std::forward<Equation>(eq), dest, src,
                            pressure_dest, efficiency);
    } else {
      const auto old_velocity = Velocity(eq, src);
      const double rhoE_kin = KineticEnergy(src.density, src.momentum);
      constexpr int N = std::decay_t<Equation>::Rank();
      dest = src;
      dest.momentum = Array<double, N, 1>::Zero();
      dest.energy = src.energy - rhoE_kin;
      const double h_before = (dest.energy + dest.pressure) / dest.density;
      SetIsentropicPressure(eq, dest, dest, pressure_dest);
      const double h_after = (dest.energy + dest.pressure) / dest.density;
      const double h_diff = h_before - h_after;
      const double e_kin_new = 2.0 * efficiency * std::abs(h_diff) +
                               old_velocity.matrix().squaredNorm();
      FUB_ASSERT(e_kin_new >= 0);
      const int sign = (h_diff >= 0) - (h_diff < 0);
      dest.momentum[0] = dest.density * (sign * std::sqrt(e_kin_new));
      dest.energy = dest.energy + KineticEnergy(dest.density, dest.momentum);
    }
  }
} IsentropicExpansionWithoutDissipation;

inline constexpr struct SpecificGasConstantFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<is_tag_invocable<
                SpecificGasConstantFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const noexcept(
      is_nothrow_tag_invocable<SpecificGasConstantFn, Equation, State>::value)
      -> tag_invoke_result_t<SpecificGasConstantFn, Equation, State> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state));
  }
} SpecificGasConstant;

inline constexpr struct SpeciesFn {
  template <typename Equation, typename State, typename... Indices,
            typename = std::enable_if_t<is_tag_invocable<
                SpeciesFn, Equation, State, Indices...>::value>>
  constexpr auto
  operator()(Equation&& eq, State&& state, Indices... d) const noexcept(
      is_nothrow_tag_invocable<SpeciesFn, Equation, State, Indices...>::value)
      -> tag_invoke_result_t<SpeciesFn, Equation, State, Indices...> {
    return fub::meta::tag_invoke(*this, std::forward<Equation>(eq),
                                 std::forward<State>(state), d...);
  }
} Species;

namespace meta {

template <typename State>
using species_t = decltype(std::declval<State>().species);

}

template <typename State>
struct state_with_species : is_detected<meta::species_t, State> {};

namespace meta {

template <typename State>
using passive_scalars_t = decltype(std::declval<State>().passive_scalars);

}

template <typename State>
struct state_with_passive_scalars
    : is_detected<meta::passive_scalars_t, State> {};

} // namespace fub::euler

#endif