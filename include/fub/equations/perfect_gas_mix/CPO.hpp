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

#ifndef FUB_EQUATINOS_PERFECT_GAS_MIX_CPOS_HPP
#define FUB_EQUATINOS_PERFECT_GAS_MIX_CPOS_HPP

namespace fub {

template <int Rank, typename State>
constexpr auto tag_invoke(tag_t<Depths>, const PerfectGasMix<Rank>& eq,
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
constexpr auto
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
constexpr const Density&
tag_invoke(tag_t<euler::Density>, const PerfectGasMix&,
           const PerfectGasMixConservative<Density, Momentum, Energy, Species>&
               q) noexcept {
  return q.density;
}

template <typename Density, typename Momentum, typename Energy,
          typename Species>
constexpr auto
tag_invoke(tag_t<euler::InternalEnergy>, const PerfectGasMix& eq,
           const PerfectGasMixConservative<Density, Momentum, Energy, Species>&
               q) noexcept {
  return (q.energy - euler::KineticEnergy(q.density, q.momentum)) / q.density;
}

template <typename Density, typename Temperature, typename MoleFractions>
constexpr auto
tag_invoke(tag_t<euler::InternalEnergy>, const PerfectGasMix& eq,
           const PerfectGasMixKineticState<Density, Temperature, MoleFractions>&
               q) noexcept {
  // Rspec = cp - cv
  // gamma = cp / cv
  // cp = cv gamma
  // Rspec = gamma cv - cv
  // Rspec = (gamma - 1) cv
  // Rpsec /(gamma - 1) = cv
  const double cv = eq.Rspec * eq.gamma_minus_1_inv;
  return cv * q.temperature;
}

constexpr void tag_invoke(tag_t<euler::CompleteFromKineticState>,
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

constexpr void tag_invoke(tag_t<euler::KineticStateFromComplete>,
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
constexpr const Momentum&
tag_invoke(tag_t<euler::Momentum>, const PerfectGasMix&,
           const PerfectGasMixConservative<Density, Momentum, Energy, Species>&
               q) noexcept {
  return q.momentum;
}

template <typename Density, typename Momentum, typename Energy,
          typename Species>
constexpr decltype(auto) tag_invoke(
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
constexpr Momentum
tag_invoke(tag_t<euler::Velocity>, const PerfectGasMix&,
           const PerfectGasMixConservative<Density, Momentum, Energy, Species>&
               q) noexcept {
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
constexpr auto tag_invoke(
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
constexpr const Energy&
tag_invoke(tag_t<euler::Energy>, const PerfectGasMix&,
           const PerfectGasMixConservative<Density, Momentum, Energy, Species>&
               q) noexcept {
  return q.energy;
}

template <typename Density, typename Momentum, typename Energy,
          typename Species>
constexpr const Species&
tag_invoke(tag_t<euler::Species>, const PerfectGasMix&,
           const PerfectGasMixConservative<Density, Momentum, Energy, Species>&
               q) noexcept {
  return q.species;
}

template <typename Density, typename Velocity, typename Pressure,
          typename Species>
constexpr const Species&
tag_invoke(tag_t<euler::Species>, const PerfectGasMix&,
           const PerfectGasMixPrimitive<Density, Velocity, Pressure, Species>&
               q) noexcept {
  return q.species;
}

template <typename Density, typename Momentum, typename Energy,
          typename Species>
constexpr decltype(auto) tag_invoke(
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
constexpr const Pressure&
tag_invoke(tag_t<euler::Pressure>, const PerfectGasMix&,
           const PerfectGasMixComplete<Density, Momentum, Energy, Species,
                                       Pressure, SpeedOfSound>& q) noexcept {
  return q.pressure;
}

template <typename Density, typename Momentum, typename Energy,
          typename Species, typename Pressure, typename SpeedOfSound>
constexpr const SpeedOfSound&
tag_invoke(tag_t<euler::SpeedOfSound>, const PerfectGasMix&,
           const PerfectGasMixComplete<Density, Momentum, Energy, Species,
                                       Pressure, SpeedOfSound>& q) noexcept {
  return q.speed_of_sound;
}

template <typename Density, typename Momentum, typename Energy,
          typename Species, typename Pressure, typename SpeedOfSound>
constexpr auto
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

} // namespace fub

#endif