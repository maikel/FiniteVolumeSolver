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

inline constexpr struct GammaFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<GammaFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const
      noexcept(is_nothrow_tag_invocable<GammaFn, Equation, State>::value)
          -> tag_invoke_result_t<GammaFn, Equation, State> {
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
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
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
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
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
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
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
                           std::forward<State>(state), d...);
  }
} Velocity;

inline constexpr struct SpeedOfSoundFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<SpeedOfSoundFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const
      noexcept(is_nothrow_tag_invocable<SpeedOfSoundFn, Equation, State>::value)
          -> tag_invoke_result_t<SpeedOfSoundFn, Equation, State> {
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
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
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
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
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
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
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
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
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
                           std::forward<State>(state));
  }
} Temperature;

inline constexpr struct InternalEnergyFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<InternalEnergyFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const
      noexcept(is_nothrow_tag_invocable<InternalEnergyFn, Equation, State>::value)
          -> tag_invoke_result_t<InternalEnergyFn, Equation, State> {
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
                           std::forward<State>(state));
  }
} InternalEnergy;

inline constexpr struct MoleFractionsFn {
  template <typename Equation, typename Dest, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<MoleFractionsFn, Equation, Dest, State>::value>>
  constexpr auto operator()(Equation&& eq, Dest&& dest, State&& state) const
      noexcept(is_nothrow_tag_invocable<MoleFractionsFn, Equation, State>::value)
          -> tag_invoke_result_t<MoleFractionsFn, Equation, State> {
    return fub::tag_invoke(*this, std::forward<Equation>(eq), std::forward<Dest>(dest),
                           std::forward<State>(state));
  }
} MoleFractions;


inline constexpr struct KineticStateFromCompleteFn {
  template <typename Equation, typename Dest, typename State,
            typename = std::enable_if_t<
                is_tag_invocable<KineticStateFromCompleteFn, Equation, Dest, State>::value>>
  constexpr auto operator()(Equation&& eq, Dest&& dest, State&& state) const
      noexcept(is_nothrow_tag_invocable<KineticStateFromCompleteFn, Equation, Dest, State>::value)
          -> tag_invoke_result_t<KineticStateFromCompleteFn, Equation, Dest, State> {
    return fub::tag_invoke(*this, std::forward<Equation>(eq), std::forward<Dest>(dest),
                           std::forward<State>(state));
  }
} KineticStateFromComplete;

inline constexpr struct CompleteFromKineticStateFn {
  template <typename Equation, typename Dest, typename State, typename Velocity,
            typename = std::enable_if_t<
                is_tag_invocable<CompleteFromKineticStateFn, Equation, Dest, State, Velocity>::value>>
  constexpr auto operator()(Equation&& eq, Dest&& dest, State&& state, Velocity&& velocity) const
      noexcept(is_nothrow_tag_invocable<CompleteFromKineticStateFn, Equation, Dest, State, Velocity>::value)
          -> tag_invoke_result_t<CompleteFromKineticStateFn, Equation, Dest, State, Velocity> {
    return fub::tag_invoke(*this, std::forward<Equation>(eq), std::forward<Dest>(dest),
                           std::forward<State>(state), std::forward<Velocity>(velocity));
  }
} CompleteFromKineticState;

inline constexpr struct SetIsentropicPressureFn {
  template <typename... Args,
            typename = std::enable_if_t<
                is_tag_invocable<SetIsentropicPressureFn, Args...>::value>>
  constexpr auto operator()(Args&&... args) const
      noexcept(is_nothrow_tag_invocable<SetIsentropicPressureFn, Args...>::value)
          -> tag_invoke_result_t<SetIsentropicPressureFn, Args...> {
    return fub::tag_invoke(*this, std::forward<Args>(args)...);
  }
} SetIsentropicPressure;

inline constexpr struct SpecificGasConstantFn {
  template <typename Equation, typename State,
            typename = std::enable_if_t<is_tag_invocable<
                SpecificGasConstantFn, Equation, State>::value>>
  constexpr auto operator()(Equation&& eq, State&& state) const noexcept(
      is_nothrow_tag_invocable<SpecificGasConstantFn, Equation, State>::value)
      -> tag_invoke_result_t<SpecificGasConstantFn, Equation, State> {
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
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
    return fub::tag_invoke(*this, std::forward<Equation>(eq),
                           std::forward<State>(state), d...);
  }
} Species;

template <typename Equation, typename State>
struct state_with_species : is_invocable<decltype(fub::euler::Species),
                                         const Equation&, const State&> {};

} // namespace fub::euler

#endif