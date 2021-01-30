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

#include "fub/equations/PerfectGasMix.hpp"

namespace fub {

template <int Dim>
void PerfectGasMix<Dim>::Flux(Conservative& flux, const Complete& state,
                              Direction dir) const noexcept {
  FUB_ASSERT(state.density > 0);
  const int d0 = static_cast<int>(dir);
  const double velocity = state.momentum[d0] / state.density;
  flux.density = state.momentum[d0];
  for (int d = 0; d < Dim; ++d) {
    flux.momentum[d] = velocity * state.momentum[d];
  }
  flux.momentum[d0] += state.pressure;
  flux.energy = velocity * (state.energy + state.pressure);
  for (int s = 0; s < n_species; ++s) {
    flux.species[s] = velocity * state.species[s];
  }
}

template <int Dim>
void PerfectGasMix<Dim>::Flux(ConservativeArray& flux,
                              const CompleteArray& state,
                              Direction dir) const noexcept {
  const int d0 = static_cast<int>(dir);
  const Array1d velocity = state.momentum.row(d0) / state.density;
  flux.density = state.momentum.row(d0);
  for (int d = 0; d < Dim; ++d) {
    flux.momentum.row(d) = velocity * state.momentum.row(d);
  }
  flux.momentum.row(d0) += state.pressure;
  flux.energy = velocity * (state.energy + state.pressure);
  for (int s = 0; s < flux.species.rows(); ++s) {
    flux.species.row(s) = velocity * state.species.row(s);
  }
}

template <int Dim>
void PerfectGasMix<Dim>::Flux(ConservativeArray& flux,
                              const CompleteArray& state, MaskArray mask,
                              Direction dir) const noexcept {
  const int d0 = static_cast<int>(dir);
  mask = mask && (state.density > 0.0);
  const Array1d density = mask.select(state.density, 1.0);
  FUB_ASSERT((density > 0.0).all());
  const Array1d velocity = state.momentum.row(d0) / density;
  flux.density = state.momentum.row(d0);
  for (int d = 0; d < Dim; ++d) {
    flux.momentum.row(d) = velocity * state.momentum.row(d);
  }
  flux.momentum.row(d0) += state.pressure;
  flux.energy = velocity * (state.energy + state.pressure);
  for (int s = 0; s < flux.species.rows(); ++s) {
    flux.species.row(s) = velocity * state.species.row(s);
  }
}

template <int Dim>
void PerfectGasMix<Dim>::CompleteFromCons(
    Complete& complete,
    const ConservativeBase<PerfectGasMix>& cons) const noexcept {
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  for (int s = 0; s < complete.species.size(); ++s) {
    complete.species[s] = cons.species[s];
  }
  const double e_kin = euler::KineticEnergy(cons.density, cons.momentum);
  FUB_ASSERT(e_kin < cons.energy);
  const double e_int = cons.energy - e_kin;
  complete.pressure = e_int / gamma_minus_one_inv;
  FUB_ASSERT(complete.pressure > 0.0 && complete.density > 0.0);
  complete.speed_of_sound =
      std::sqrt(gamma * complete.pressure / complete.density);
}

template <int Dim>
void PerfectGasMix<Dim>::CompleteFromCons(
    CompleteArray& complete,
    const ConservativeArrayBase<PerfectGasMix>& cons) const noexcept {
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  for (int s = 0; s < complete.species.rows(); ++s) {
    complete.species.row(s) = cons.species.row(s);
  }
  const Array1d e_kin = euler::KineticEnergy(cons.density, cons.momentum);
  const Array1d e_int = cons.energy - e_kin;
  complete.pressure = e_int * gamma_minus_one;
  complete.speed_of_sound =
      (gamma_array_ * complete.pressure / complete.density).sqrt();
}

template <int Dim>
void PerfectGasMix<Dim>::CompleteFromCons(
    CompleteArray& complete, const ConservativeArrayBase<PerfectGasMix>& cons,
    MaskArray mask) const noexcept {
  Array1d zero = Array1d::Zero();
  mask = mask && (cons.density > 0.0);
  complete.density = mask.select(cons.density, zero);
  for (int d = 0; d < Dim; ++d) {
    complete.momentum.row(d) = mask.select(cons.momentum.row(d), zero);
  }
  complete.energy = mask.select(cons.energy, zero);
  for (int s = 0; s < complete.species.rows(); ++s) {
    complete.species.row(s) = mask.select(cons.species.row(s), zero);
  }
  const Array1d e_kin = euler::KineticEnergy(cons.density, cons.momentum, mask);
  const Array1d e_int = cons.energy - e_kin;
  complete.pressure = mask.select(e_int / gamma_minus_one_inv, zero);
  Array1d safe_density = mask.select(complete.density, 1.0);
  FUB_ASSERT((safe_density > 0.0).all());
  complete.speed_of_sound = mask.select(
      (gamma_array_ * complete.pressure / safe_density).sqrt(), zero);
}

template <int Dim>
Complete<PerfectGasMix<Dim>> PerfectGasMix<Dim>::CompleteFromPrim(
    double rho, const Array<double, Dim, 1>& v, double p,
    const Array<double, -1, 1>& species) const noexcept {
  Complete q{*this};
  q.density = rho;
  q.momentum = rho * v;
  q.pressure = p;
  for (int s = 0; s < n_species; ++s) {
    q.species[s] = q.density * species[s];
  }
  const double e_kin = euler::KineticEnergy(q.density, q.momentum);
  q.energy = e_kin + p * gamma_minus_one_inv;
  q.speed_of_sound = std::sqrt(gamma * q.pressure / q.density);
  return q;
}

template <int Dim>
Array<double, Dim, 1>
PerfectGasMix<Dim>::Velocity(const Complete& complete) const noexcept {
  Array<double, Dim, 1> u = complete.momentum / complete.density;
  return u;
}

template <int Dim>
Array<double, Dim>
PerfectGasMix<Dim>::Velocity(const CompleteArray& complete) const noexcept {
  Array<double, Dim> u = complete.momentum;
  for (int i = 0; i < Dim; ++i) {
    u.row(i) /= complete.density;
  }
  return u;
}

template <int Dim>
double PerfectGasMix<Dim>::Machnumber(const Complete& complete) const noexcept {
  double u = Velocity(complete).matrix().norm();
  return u / complete.speed_of_sound;
}

template <int Dim>
double
PerfectGasMix<Dim>::Temperature(const Complete& complete) const noexcept {
  return complete.pressure / (complete.density * Rspec);
}

template <int Dim>
Array1d
PerfectGasMix<Dim>::Temperature(const CompleteArray& complete) const noexcept {
  return complete.pressure / (complete.density * Rspec);
}

template struct PerfectGasMix<1>;
template struct PerfectGasMix<2>;
template struct PerfectGasMix<3>;

void Rotate(Conservative<PerfectGasMix<2>>& rotated,
            const Conservative<PerfectGasMix<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation,
            const PerfectGasMix<2>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
  for (int s = 0; s < rotated.species.size(); ++s) {
    rotated.species[s] = state.species[s];
  }
}

void Rotate(Complete<PerfectGasMix<2>>& rotated,
            const Complete<PerfectGasMix<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation,
            const PerfectGasMix<2>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.pressure = state.pressure;
  rotated.speed_of_sound = state.speed_of_sound;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
  for (int s = 0; s < rotated.species.size(); ++s) {
    rotated.species[s] = state.species[s];
  }
}

void Rotate(Conservative<PerfectGasMix<3>>& rotated,
            const Conservative<PerfectGasMix<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation,
            const PerfectGasMix<3>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
  for (int s = 0; s < rotated.species.size(); ++s) {
    rotated.species[s] = state.species[s];
  }
}

void Rotate(Complete<PerfectGasMix<3>>& rotated,
            const Complete<PerfectGasMix<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation,
            const PerfectGasMix<3>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.pressure = state.pressure;
  rotated.speed_of_sound = state.speed_of_sound;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
  for (int s = 0; s < rotated.species.size(); ++s) {
    rotated.species[s] = state.species[s];
  }
}

void Reflect(Complete<PerfectGasMix<1>>& reflected,
             const Complete<PerfectGasMix<1>>& state,
             const Eigen::Matrix<double, 1, 1>& normal,
             const PerfectGasMix<1>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
  for (int s = 0; s < reflected.species.size(); ++s) {
    reflected.species[s] = state.species[s];
  }
}

void Reflect(Complete<PerfectGasMix<2>>& reflected,
             const Complete<PerfectGasMix<2>>& state,
             const Eigen::Vector2d& normal, const PerfectGasMix<2>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
  for (int s = 0; s < reflected.species.size(); ++s) {
    reflected.species[s] = state.species[s];
  }
}

void Reflect(Complete<PerfectGasMix<3>>& reflected,
             const Complete<PerfectGasMix<3>>& state,
             const Eigen::Vector3d& normal, const PerfectGasMix<3>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
  for (int s = 0; s < reflected.species.size(); ++s) {
    reflected.species[s] = state.species[s];
  }
}

void Reflect(Conservative<PerfectGasMix<1>>& reflected,
             const Conservative<PerfectGasMix<1>>& state,
             const Eigen::Matrix<double, 1, 1>& normal,
             const PerfectGasMix<1>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
  for (int s = 0; s < reflected.species.size(); ++s) {
    reflected.species[s] = state.species[s];
  }
}

void Reflect(Conservative<PerfectGasMix<2>>& reflected,
             const Conservative<PerfectGasMix<2>>& state,
             const Eigen::Vector2d& normal, const PerfectGasMix<2>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
  for (int s = 0; s < reflected.species.size(); ++s) {
    reflected.species[s] = state.species[s];
  }
}

void Reflect(Conservative<PerfectGasMix<3>>& reflected,
             const Conservative<PerfectGasMix<3>>& state,
             const Eigen::Vector3d& normal, const PerfectGasMix<3>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
  for (int s = 0; s < reflected.species.size(); ++s) {
    reflected.species[s] = state.species[s];
  }
}

} // namespace fub
