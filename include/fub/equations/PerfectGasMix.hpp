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

#ifndef FUB_EQUATIONS_PERFECT_GAS_MIX_HPP
#define FUB_EQUATIONS_PERFECT_GAS_MIX_HPP

#include "fub/State.hpp"
#include "fub/StateArray.hpp"

#include "fub/CompleteFromCons.hpp"
#include "fub/Equation.hpp"
#include "fub/ext/Eigen.hpp"

#include <array>

namespace fub {
template <int Rank> struct PerfectGasMix;

/// This is a template class for constructing conservative states for the
/// perfect gas equations.
template <typename Density, typename Momentum, typename Energy,
          typename Species>
struct PerfectGasMixConservative { // Zenker: renamed PerfectGasConservativeMix --> PerfectGasMixConservative
  Density density;
  Momentum momentum;
  Energy energy;
  Species species;
};

template <int Rank> // Zenker: renamed PerfectGasConsMixShape --> PerfectGasMixConsShape
using PerfectGasMixConsShape =
    PerfectGasMixConservative<ScalarDepth, VectorDepth<Rank>, ScalarDepth,
                              VectorDepth<-1>>;

namespace meta {
template <int R> struct Rank<PerfectGasMixConsShape<R>> : int_constant<R> {};
} // namespace meta

// We "register" the conservative state with our framework.
// This enables us to name and iterate over all member variables in a given
// conservative state.
template <typename... Xs> struct StateTraits<PerfectGasMixConservative<Xs...>> {
  static constexpr auto names =
      std::make_tuple("Density", "Momentum", "Energy", "Species");

  static constexpr auto pointers_to_member =
      std::make_tuple(&PerfectGasMixConservative<Xs...>::density,
                      &PerfectGasMixConservative<Xs...>::momentum,
                      &PerfectGasMixConservative<Xs...>::energy,
                      &PerfectGasMixConservative<Xs...>::species);

  template <int Rank> using Depths = PerfectGasMixConsShape<Rank>;

  template <int Rank> using Equation = PerfectGasMix<Rank>;
};

template <typename Density, typename Velocity, typename Pressure,
          typename Species>
struct PerfectGasMixPrimitive {
  Density density;
  Velocity velocity;
  Pressure pressure;
  Species species;
};

template <int Rank>
using PerfectGasMixPrimShape =
    PerfectGasMixPrimitive<ScalarDepth, VectorDepth<Rank>, ScalarDepth,
                           VectorDepth<-1>>;

namespace meta {
template <int R> struct Rank<PerfectGasMixPrimShape<R>> : int_constant<R> {};
} // namespace meta

template <typename... Xs> struct StateTraits<PerfectGasMixPrimitive<Xs...>> {
  static constexpr auto names =
      std::make_tuple("Density", "Velocity", "Pressure", "Species");

  static constexpr auto pointers_to_member =
      std::make_tuple(&PerfectGasMixPrimitive<Xs...>::density,
                      &PerfectGasMixPrimitive<Xs...>::velocity,
                      &PerfectGasMixPrimitive<Xs...>::pressure,
                      &PerfectGasMixPrimitive<Xs...>::species);

  template <int Rank> using Depths = PerfectGasMixPrimShape<Rank>;

  template <int Rank> using Equation = PerfectGasMix<Rank>;
};

template <typename Minus, typename Zero, typename Plus, typename Species>
struct PerfectGasMixCharacteristics {
  Minus minus;
  Zero zero;
  Plus plus;
  Species species;
};

template <int Rank>
using PerfectGasMixCharShape =
    PerfectGasMixCharacteristics<ScalarDepth, VectorDepth<Rank>, ScalarDepth,
                                 VectorDepth<-1>>;

namespace meta {
template <int R> struct Rank<PerfectGasMixCharShape<R>> : int_constant<R> {};
} // namespace meta

template <typename... Xs>
struct StateTraits<PerfectGasMixCharacteristics<Xs...>> {
  static constexpr auto names =
      std::make_tuple("Minus", "Zero", "Plus", "Species");

  static constexpr auto pointers_to_member =
      std::make_tuple(&PerfectGasMixCharacteristics<Xs...>::minus,
                      &PerfectGasMixCharacteristics<Xs...>::zero,
                      &PerfectGasMixCharacteristics<Xs...>::plus,
                      &PerfectGasMixCharacteristics<Xs...>::species);

  template <int Rank> using Depths = PerfectGasMixCharShape<Rank>;

  template <int Rank> using Equation = PerfectGasMix<Rank>;
};

template <typename Density, typename Momentum, typename Energy,
          typename Species, typename Pressure, typename SpeedOfSound>
struct PerfectGasMixComplete
    : PerfectGasMixConservative<Density, Momentum, Energy, Species> {
  Pressure pressure;
  SpeedOfSound speed_of_sound;
};

template <int Rank>
using PerfectGasMixCompleteShape =
    PerfectGasMixComplete<ScalarDepth, VectorDepth<Rank>, ScalarDepth,
                          VectorDepth<-1>, ScalarDepth, ScalarDepth>;

namespace meta {
template <int R>
struct Rank<PerfectGasMixCompleteShape<R>> : int_constant<R> {};
} // namespace meta

// We "register" the complete state with our framework.
// This enables us to name and iterate over all member variables in a given
// compete state.
template <typename... Xs> struct StateTraits<PerfectGasMixComplete<Xs...>> {
  static constexpr auto names = std::make_tuple("Density", "Momentum", "Energy",
                                                "Species"
                                                "Pressure",
                                                "SpeedOfSound");
  static constexpr auto pointers_to_member = std::make_tuple(
      &PerfectGasMixComplete<Xs...>::density, &PerfectGasMixComplete<Xs...>::momentum,
      &PerfectGasMixComplete<Xs...>::energy, &PerfectGasMixComplete<Xs...>::species,
      &PerfectGasMixComplete<Xs...>::pressure,
      &PerfectGasMixComplete<Xs...>::speed_of_sound);

  template <int Rank> using Depths = PerfectGasMixCompleteShape<Rank>;

  template <int Rank> using Equation = PerfectGasMix<Rank>;
};

template <int N> struct PerfectGasMix {
  using ConservativeDepths = PerfectGasMixConsShape<N>;
  using CompleteDepths = PerfectGasMixCompleteShape<N>;
  using PrimitiveDepths = PerfectGasMixPrimShape<N>;
  using CharacteristicsDepths = PerfectGasMixCharShape<N>;

  using Conservative = ::fub::Conservative<PerfectGasMix<N>>;
  using Complete = ::fub::Complete<PerfectGasMix<N>>;
  using ConservativeArray = ::fub::ConservativeArray<PerfectGasMix<N>>;
  using CompleteArray = ::fub::CompleteArray<PerfectGasMix<N>>;

  static constexpr int Rank() noexcept { return N; }

  void Flux(Conservative& flux, const Complete& state,
            [[maybe_unused]] Direction dir = Direction::X) const noexcept;

  void Flux(ConservativeArray& flux, const CompleteArray& state,
            [[maybe_unused]] Direction dir) const noexcept;

  void Flux(ConservativeArray& flux, const CompleteArray& state, MaskArray mask,
            [[maybe_unused]] Direction dir) const noexcept;

  void
  CompleteFromCons(Complete& complete,
                   const ConservativeBase<PerfectGasMix>& cons) const noexcept;

  void CompleteFromCons(
      CompleteArray& complete,
      const ConservativeArrayBase<PerfectGasMix>& cons) const noexcept;

  void CompleteFromCons(CompleteArray& complete,
                        const ConservativeArrayBase<PerfectGasMix>& cons,
                        MaskArray mask) const noexcept;

  Complete CompleteFromPrim(double density, const Array<double, N, 1>& u,
                            double pressure, const Array<double, -1, 1>& species) const noexcept;

  // CompleteArray CompleteFromPrim(Array1d density, const Array<double, N>& u,
  //                                Array1d pressure) const noexcept;

  // CompleteArray CompleteFromPrim(Array1d density, const Array<double, N>& u,
  //                                Array1d pressure,
  //                                const MaskArray& mask) const noexcept;

  Array<double, N, 1> Velocity(const Complete& q) const noexcept;
  Array<double, N> Velocity(const CompleteArray& q) const noexcept;

  double Machnumber(const Complete& q) const noexcept;

  double Temperature(const Complete& q) const noexcept;
  Array1d Temperature(const CompleteArray& q) const noexcept;

  int n_species{0};

  double Rspec{287.058};

  double gamma{1.4};
  double gamma_minus_1_inv{1.0 / (gamma - 1.0)};

  Array1d gamma_array_{Array1d::Constant(gamma)};
  Array1d gamma_minus_1_inv_array_{Array1d::Constant(gamma_minus_1_inv)};
};

// We define this class only for dimensions 1 to 3.
// The definitions will be found in its source file PerfectGasMix.cpp
extern template struct PerfectGasMix<1>;
extern template struct PerfectGasMix<2>;
extern template struct PerfectGasMix<3>;

namespace detail {
template <typename State, int Dim>
struct DepthsImpl<State, PerfectGasMix<Dim>> {
  constexpr ToConcreteDepths<typename StateTraits<State>::template Depths<Dim>>
  operator()(const PerfectGasMix<Dim>& equation) const noexcept {
    ToConcreteDepths<typename StateTraits<State>::template Depths<Dim>> depths{};
    depths.species = equation.n_species;
    return depths;
  }
};

}

template <int Rank>
void CompleteFromPrim(const PerfectGasMix<Rank>& equation,
                      Complete<PerfectGasMix<Rank>>& complete,
                      const Primitive<PerfectGasMix<Rank>>& prim) {
  complete =
      equation.CompleteFromPrim(prim.density, prim.velocity, prim.pressure, prim.species);
}

template <int Rank>
void PrimFromComplete(const PerfectGasMix<Rank>& equation,
                      Primitive<PerfectGasMix<Rank>>& prim,
                      const Complete<PerfectGasMix<Rank>>& complete) {
  prim.density = complete.density;
  prim.pressure = complete.pressure;
  prim.velocity = complete.momentum / complete.density;
  prim.species = complete.species / complete.density;
}

/// @{
/// \brief Defines how to rotate a given state of the euler equations.
///
/// This function is needed when computing the reference state in the boundary
/// flux of the cut-cell stabilizations.
void Rotate(Conservative<PerfectGasMix<2>>& rotated,
            const Conservative<PerfectGasMix<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation,
            const PerfectGasMix<2>&);

void Rotate(Complete<PerfectGasMix<2>>& rotated,
            const Complete<PerfectGasMix<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation,
            const PerfectGasMix<2>&);

void Rotate(Conservative<PerfectGasMix<3>>& rotated,
            const Conservative<PerfectGasMix<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation,
            const PerfectGasMix<3>&);

void Rotate(Complete<PerfectGasMix<3>>& rotated,
            const Complete<PerfectGasMix<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation,
            const PerfectGasMix<3>&);
/// @}

void Reflect(Conservative<PerfectGasMix<1>>& reflected,
             const Conservative<PerfectGasMix<1>>& state,
             const Eigen::Matrix<double, 1, 1>& normal,
             const PerfectGasMix<1>& gas);

void Reflect(Conservative<PerfectGasMix<2>>& reflected,
             const Conservative<PerfectGasMix<2>>& state,
             const Eigen::Vector2d& normal, const PerfectGasMix<2>& gas);

void Reflect(Conservative<PerfectGasMix<3>>& reflected,
             const Conservative<PerfectGasMix<3>>& state,
             const Eigen::Vector3d& normal, const PerfectGasMix<3>& gas);

void Reflect(Complete<PerfectGasMix<1>>& reflected,
             const Complete<PerfectGasMix<1>>& state,
             const Eigen::Matrix<double, 1, 1>& normal,
             const PerfectGasMix<1>& gas);

void Reflect(Complete<PerfectGasMix<2>>& reflected,
             const Complete<PerfectGasMix<2>>& state,
             const Eigen::Vector2d& normal, const PerfectGasMix<2>& gas);

void Reflect(Complete<PerfectGasMix<3>>& reflected,
             const Complete<PerfectGasMix<3>>& state,
             const Eigen::Vector3d& normal, const PerfectGasMix<3>& gas);

} // namespace fub

#endif
