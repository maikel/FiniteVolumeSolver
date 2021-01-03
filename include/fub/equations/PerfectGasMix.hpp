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

#include "fub/equations/EulerEquation.hpp"

#include <array>

namespace fub {

template <int Rank> struct PerfectGasMix;

/// This is a template class for constructing conservative states for the
/// perfect gas equations.
template <typename Density, typename Momentum, typename Energy,
          typename Species>
struct PerfectGasMixConservative {
  Density density;
  Momentum momentum;
  Energy energy;
  Species species;
};

template <int Rank>
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

template <typename Density, typename Temperature, typename MoleFractions>
struct PerfectGasMixKineticState {
  Density density;
  Temperature temperature;
  MoleFractions mole_fractions;
};

using PerfectGasMixKineticStateShape =
    PerfectGasMixKineticState<ScalarDepth, ScalarDepth, VectorDepth<-1>>;

namespace meta {
template <> struct Rank<PerfectGasMixKineticStateShape> : int_constant<1> {};
} // namespace meta

template <typename... Xs> struct StateTraits<PerfectGasMixKineticState<Xs...>> {
  static constexpr auto names =
      std::make_tuple("Density", "Temperature", "MoleFractions");

  static constexpr auto pointers_to_member =
      std::make_tuple(&PerfectGasMixKineticState<Xs...>::density,
                      &PerfectGasMixKineticState<Xs...>::temperature,
                      &PerfectGasMixKineticState<Xs...>::mole_fractions);

  template <int Rank> using Depths = PerfectGasMixKineticStateShape;

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
  static constexpr auto names = std::make_tuple(
      "Density", "Momentum", "Energy", "Species", "Pressure", "SpeedOfSound");
  static constexpr auto pointers_to_member =
      std::make_tuple(&PerfectGasMixComplete<Xs...>::density,
                      &PerfectGasMixComplete<Xs...>::momentum,
                      &PerfectGasMixComplete<Xs...>::energy,
                      &PerfectGasMixComplete<Xs...>::species,
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
  using KineticStateDepths = PerfectGasMixKineticStateShape;

  using Conservative = ::fub::Conservative<PerfectGasMix<N>>;
  using Complete = ::fub::Complete<PerfectGasMix<N>>;
  using ConservativeArray = ::fub::ConservativeArray<PerfectGasMix<N>>;
  using CompleteArray = ::fub::CompleteArray<PerfectGasMix<N>>;
  using KineticState = ::fub::KineticState<PerfectGasMix<N>>;

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
                            double pressure,
                            const Array<double, -1, 1>& species) const noexcept;

  // CompleteArray CompleteFromPrim(Array1d density, const Array<double, N>& u,
  //                                Array1d pressure,
  //                                const MaskArray& mask) const noexcept;

  Array<double, N, 1> Velocity(const Complete& q) const noexcept;
  Array<double, N> Velocity(const CompleteArray& q) const noexcept;

  double Machnumber(const Complete& q) const noexcept;

  double Temperature(const Complete& q) const noexcept;
  Array1d Temperature(const CompleteArray& q) const noexcept;

  int n_species{0};

  // Rspec = cp - cv
  // gamma = cp / cv
  // cp = cv gamma
  // Rspec = gamma cv - cv
  // Rspec = (gamma - 1) cv
  // Rpsec /(gamma - 1) = cv
  double Rspec{1.};
  double gamma{1.28};
  double gamma_minus_1_inv{1.0 / (gamma - 1.0)};

  Array1d gamma_array_{Array1d::Constant(gamma)};
  Array1d gamma_minus_1_inv_array_{Array1d::Constant(gamma_minus_1_inv)};

private:
  template <typename State>
  friend constexpr auto tag_invoke(tag_t<Depths>, const PerfectGasMix& eq,
                                   Type<State>) noexcept {
    using Depths = typename State::Traits::template Depths<N>;
    ToConcreteDepths<Depths> depths{};
    if constexpr (std::is_same_v<Depths, PerfectGasMixKineticStateShape>) {
      depths.mole_fractions = eq.n_species + 1;
    } else {
      depths.species = eq.n_species;
    }
    return depths;
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  friend auto
  tag_invoke(tag_t<euler::Gamma>, const PerfectGasMix& eq,
             const PerfectGasMixConservative<Density, Momentum, Energy,
                                             Species>&) noexcept {
    if constexpr (std::is_same_v<double, Density>) {
      return eq.gamma;
    } else {
      return Array1d::Constant(eq.gamma);
    }
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  friend const Density&
  tag_invoke(tag_t<euler::Density>, const PerfectGasMix&,
             const PerfectGasMixConservative<Density, Momentum, Energy,
                                             Species>& q) noexcept {
    return q.density;
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  friend auto
  tag_invoke(tag_t<euler::InternalEnergy>, const PerfectGasMix& eq,
             const PerfectGasMixConservative<Density, Momentum, Energy,
                                             Species>& q) noexcept {
    return (q.energy - euler::KineticEnergy(q.density, q.momentum)) / q.density;
  }

  template <typename Density, typename Temperature, typename MoleFractions>
  friend auto
  tag_invoke(tag_t<euler::InternalEnergy>, const PerfectGasMix& eq,
             const PerfectGasMixKineticState<Density, Temperature,
                                             MoleFractions>& q) noexcept {
    // Rspec = cp - cv
    // gamma = cp / cv
    // cp = cv gamma
    // Rspec = gamma cv - cv
    // Rspec = (gamma - 1) cv
    // Rpsec /(gamma - 1) = cv
    const double cv = eq.Rspec * eq.gamma_minus_1_inv;
    return cv * q.temperature;
  }

  friend void tag_invoke(tag_t<euler::CompleteFromKineticState>,
                         const PerfectGasMix& eq, Complete& q,
                         const KineticState& kin,
                         const Array<double, N, 1>& u) noexcept {
    q.density = kin.density;
    q.momentum = q.density * u;
    const double e = euler::InternalEnergy(eq, kin);
    q.energy = q.density * (e + 0.5 * u.matrix().squaredNorm());
    const double RT = eq.Rspec * kin.temperature;
    q.pressure = RT * q.density;
    q.speed_of_sound = std::sqrt(eq.gamma * RT);
    const double sum = kin.mole_fractions.sum();
    for (int i = 0; i < eq.n_species; i++) {
      q.species[i] = q.density * kin.mole_fractions[i] / sum;
    }
  }

  friend void tag_invoke(tag_t<euler::KineticStateFromComplete>,
                         const PerfectGasMix& eq, KineticState& kin,
                         const Complete& q) noexcept {
    kin.density = q.density;
    kin.temperature = euler::Temperature(eq, q);
    double sum = 0.0;
    for (int i = 0; i < eq.n_species; ++i) {
      kin.mole_fractions[i] = q.species[i] / q.density;
      sum += kin.mole_fractions[i];
    }
    kin.mole_fractions[eq.n_species] = std::max(0.0, 1.0 - sum);
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  friend const Momentum&
  tag_invoke(tag_t<euler::Momentum>, const PerfectGasMix&,
             const PerfectGasMixConservative<Density, Momentum, Energy,
                                             Species>& q) noexcept {
    return q.momentum;
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  friend decltype(auto) tag_invoke(
      tag_t<euler::Momentum>, const PerfectGasMix&,
      const PerfectGasMixConservative<Density, Momentum, Energy, Species>& q,
      int d) noexcept {
    if constexpr (std::is_same_v<double, Density>) {
      return q.momentum[d];
    } else {
      return q.row(d);
    }
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  friend Momentum
  tag_invoke(tag_t<euler::Velocity>, const PerfectGasMix&,
             const PerfectGasMixConservative<Density, Momentum, Energy,
                                             Species>& q) noexcept {
    if constexpr (std::is_same_v<double, Density>) {
      return q.momentum / q.density;
    } else {
      Array<double, N> v;
      for (int i = 0; i < N; ++i) {
        v.row(i) = q.momentum.row(i) / q.density;
      }
      return v;
    }
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  friend auto tag_invoke(
      tag_t<euler::Velocity>, const PerfectGasMix&,
      const PerfectGasMixConservative<Density, Momentum, Energy, Species>& q,
      int d) noexcept {
    if constexpr (std::is_same_v<double, Density>) {
      return q.momentum[d] / q.density;
    } else {
      return q.row(d) / q.density;
    }
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  friend const Energy&
  tag_invoke(tag_t<euler::Energy>, const PerfectGasMix&,
             const PerfectGasMixConservative<Density, Momentum, Energy,
                                             Species>& q) noexcept {
    return q.energy;
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  friend const Species&
  tag_invoke(tag_t<euler::Species>, const PerfectGasMix&,
             const PerfectGasMixConservative<Density, Momentum, Energy,
                                             Species>& q) noexcept {
    return q.species;
  }

  template <typename Density, typename Velocity, typename Pressure,
            typename Species>
  friend const Species&
  tag_invoke(tag_t<euler::Species>, const PerfectGasMix&,
             const PerfectGasMixPrimitive<Density, Velocity, Pressure, Species>&
                 q) noexcept {
    return q.species;
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  friend decltype(auto) tag_invoke(
      tag_t<euler::Species>, const PerfectGasMix&,
      const PerfectGasMixConservative<Density, Momentum, Energy, Species>& q,
      int d) noexcept {
    if constexpr (std::is_same_v<double, Density>) {
      return q.species[d];
    } else {
      return q.species.row(d);
    }
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species, typename Pressure, typename SpeedOfSound>
  friend const Pressure&
  tag_invoke(tag_t<euler::Pressure>, const PerfectGasMix&,
             const PerfectGasMixComplete<Density, Momentum, Energy, Species,
                                         Pressure, SpeedOfSound>& q) noexcept {
    return q.pressure;
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species, typename Pressure, typename SpeedOfSound>
  friend const SpeedOfSound&
  tag_invoke(tag_t<euler::SpeedOfSound>, const PerfectGasMix&,
             const PerfectGasMixComplete<Density, Momentum, Energy, Species,
                                         Pressure, SpeedOfSound>& q) noexcept {
    return q.speed_of_sound;
  }

  template <typename Density, typename Momentum, typename Energy,
            typename Species, typename Pressure, typename SpeedOfSound>
  friend auto
  tag_invoke(tag_t<euler::SetIsentropicPressure>, const PerfectGasMix& eq,
             PerfectGasMixComplete<Density, Momentum, Energy, Species, Pressure,
                                   SpeedOfSound>& q,
             const PerfectGasMixComplete<Density, Momentum, Energy, Species,
                                         Pressure, SpeedOfSound>& q0,
             Pressure p_new) noexcept {
    const double rho_new =
        std::pow(p_new / q0.pressure, 1 / eq.gamma) * q0.density;
    const Array<double, N, 1> u0 = euler::Velocity(eq, q0);
    Array<double, N, 1> u_new = u0;
    u_new[0] = u0[0] +
        2.0 * std::sqrt(eq.gamma * q0.pressure / q0.density) *
            eq.gamma_minus_1_inv -
        2.0 * std::sqrt(eq.gamma * p_new / rho_new) * eq.gamma_minus_1_inv;
    const Array<double, N, 1> rhou_new = rho_new * u_new;
    const double rhoE_new =
        p_new * eq.gamma_minus_1_inv + euler::KineticEnergy(rho_new, rhou_new);

    q.density = rho_new;
    q.momentum = rhou_new;
    q.energy = rhoE_new;
    q.species = q0.species;
    q.pressure = p_new;
    q.speed_of_sound = std::sqrt(eq.gamma * p_new / rho_new);
  }
};

// We define this class only for dimensions 1 to 3.
// The definitions will be found in its source file PerfectGasMix.cpp
extern template struct PerfectGasMix<1>;
extern template struct PerfectGasMix<2>;
extern template struct PerfectGasMix<3>;

template <int Rank>
void CompleteFromPrim(const PerfectGasMix<Rank>& equation,
                      Complete<PerfectGasMix<Rank>>& complete,
                      const Primitive<PerfectGasMix<Rank>>& prim) {
  complete = equation.CompleteFromPrim(prim.density, prim.velocity,
                                       prim.pressure, prim.species);
}

template <int Rank>
void CompleteFromPrim(const PerfectGasMix<Rank>& equation,
                      CompleteArray<PerfectGasMix<Rank>>& complete,
                      const PrimitiveArray<PerfectGasMix<Rank>>& prim) {
  complete.density = prim.density;
  for (int i = 0; i < Rank; ++i) {
    complete.momentum.row(i) = complete.density * prim.velocity.row(i);
  }
  complete.pressure = prim.pressure;
  for (int s = 0; s < equation.n_species; ++s) {
    complete.species.row(s) = prim.density * prim.species.row(s);
  }
  const Array1d e_kin =
      euler::KineticEnergy(complete.density, complete.momentum);
  complete.energy = e_kin + complete.pressure * equation.gamma_minus_1_inv;
  complete.speed_of_sound =
      (equation.gamma * complete.pressure / complete.density).sqrt();
}

template <int Rank>
void PrimFromComplete(const PerfectGasMix<Rank>& equation,
                      Primitive<PerfectGasMix<Rank>>& prim,
                      const Complete<PerfectGasMix<Rank>>& complete) {
  prim.density = complete.density;
  prim.pressure = complete.pressure;
  prim.velocity = complete.momentum / complete.density;
  for (int i = 0; i < equation.n_species; ++i) {
    prim.species[i] = complete.species[i] / complete.density;
  }
}

template <int Rank>
void PrimFromComplete(const PerfectGasMix<Rank>& equation,
                      PrimitiveArray<PerfectGasMix<Rank>>& prim,
                      const CompleteArray<PerfectGasMix<Rank>>& complete) {
  prim.density = complete.density;
  prim.pressure = complete.pressure;
  for (int i = 0; i < Rank; ++i) {
    prim.velocity.row(i) = complete.momentum.row(i) / complete.density;
  }
  for (int i = 0; i < equation.n_species; ++i) {
    prim.species.row(i) = complete.species.row(i) / complete.density;
  }
}

template <int Rank>
void tag_invoke(tag_t<euler::CompleteFromKineticState>,
                const PerfectGasMix<Rank>& eq,
                CompleteArray<PerfectGasMix<Rank>>& q,
                const KineticStateArray<PerfectGasMix<Rank>>& kin,
                const Array<double, Rank>& u) noexcept {
  q.density = kin.density;
  for (int i = 0; i < Rank; ++i) {
    q.momentum.row(i) = q.density * u.row(i);
  }
  const Array1d e = euler::InternalEnergy(eq, kin);
  const Array1d rhoE_kin = euler::KineticEnergy(q.density, q.momentum);
  q.energy = q.density * e + rhoE_kin;
  const Array1d RT = eq.Rspec * kin.temperature;
  q.pressure = RT * q.density;
  q.speed_of_sound = Eigen::sqrt(eq.gamma * RT);
  const Array1d sum = kin.mole_fractions.colwise().sum();
  for (int i = 0; i < eq.n_species; i++) {
    q.species.row(i) = q.density * kin.mole_fractions.row(i) / sum;
  }
}

template <int Rank>
void tag_invoke(tag_t<euler::CompleteFromKineticState>,
                const PerfectGasMix<Rank>& eq, Complete<PerfectGasMix<Rank>>& q,
                const KineticState<PerfectGasMix<Rank>>& kin) noexcept {
  Array<double, Rank, 1> u = Array<double, Rank, 1>::Zero();
  euler::CompleteFromKineticState(eq, q, kin, u);
}

template <int Rank>
void tag_invoke(tag_t<euler::KineticStateFromComplete>,
                const PerfectGasMix<Rank>& eq,
                KineticStateArray<PerfectGasMix<Rank>>& kin,
                const CompleteArray<PerfectGasMix<Rank>>& q) noexcept {
  kin.density = q.density;
  kin.temperature = euler::Temperature(eq, q);
  kin.mole_fractions.setZero();
  for (int i = 0; i < eq.n_species; ++i) {
    kin.mole_fractions.row(i) = (q.species.row(i) / q.density).max(0.0);
  }
  Array1d sum = kin.mole_fractions.colwise().sum();
  kin.mole_fractions.row(eq.n_species) = (1.0 - sum).max(0.0);
  sum = kin.mole_fractions.colwise().sum();
  for (int i = 0; i < eq.n_species + 1; ++i) {
    kin.mole_fractions.row(i) /= sum;
  }
  FUB_ASSERT(kin.mole_fractions.colwise().sum().isApproxToConstant(1.0));
}

template <int Rank>
double tag_invoke(tag_t<euler::Temperature>, const PerfectGasMix<Rank>& eq,
                  const Complete<PerfectGasMix<Rank>>& q) noexcept {
  return eq.Temperature(q);
}

template <int Rank>
Array1d tag_invoke(tag_t<euler::Temperature>, const PerfectGasMix<Rank>& eq,
                   const CompleteArray<PerfectGasMix<Rank>>& q) noexcept {
  return eq.Temperature(q);
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
