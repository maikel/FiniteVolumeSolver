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

#include "fub/equations/PerfectGas.hpp"

namespace fub {
namespace {
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
} // namespace

template <int Dim>
void PerfectGas<Dim>::Flux(Conservative& flux, const Complete& state,
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
}

template <int Dim>
void PerfectGas<Dim>::Flux(ConservativeArray& flux, const CompleteArray& state,
                           Direction dir) const noexcept {
  const int d0 = static_cast<int>(dir);
  const Array1d velocity = state.momentum.row(d0) / state.density;
  flux.density = state.momentum.row(d0);
  for (int d = 0; d < Dim; ++d) {
    flux.momentum.row(d) = velocity * state.momentum.row(d);
  }
  flux.momentum.row(d0) += state.pressure;
  flux.energy = velocity * (state.energy + state.pressure);
}

template <int Dim>
void PerfectGas<Dim>::Flux(ConservativeArray& flux, const CompleteArray& state, MaskArray mask,
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
}


template <int Dim>
void PerfectGas<Dim>::CompleteFromCons(
    Complete& complete, const ConservativeBase<PerfectGas>& cons) const
    noexcept {
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  const double e_kin = KineticEnergy(cons.density, cons.momentum);
  FUB_ASSERT(e_kin < cons.energy);
  const double e_int = cons.energy - e_kin;
  complete.pressure = e_int / gamma_minus_1_inv;
  FUB_ASSERT(complete.pressure > 0.0 && complete.density > 0.0);
  complete.speed_of_sound =
      std::sqrt(gamma * complete.pressure / complete.density);
}

template <int Dim>
void PerfectGas<Dim>::CompleteFromCons(
    CompleteArray& complete,
    const ConservativeArrayBase<PerfectGas>& cons) const noexcept {
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  const Array1d e_kin = KineticEnergy(cons.density, cons.momentum);
  const Array1d e_int = cons.energy - e_kin;
  complete.pressure = e_int / gamma_minus_1_inv;
  complete.speed_of_sound =
      (gamma_array_ * complete.pressure / complete.density).sqrt();
}

template <int Dim>
void PerfectGas<Dim>::CompleteFromCons(
    CompleteArray& complete, const ConservativeArrayBase<PerfectGas>& cons,
    MaskArray mask) const noexcept {
  Array1d zero = Array1d::Zero();
  mask = mask && (cons.density > 0.0);
  complete.density = mask.select(cons.density, zero);
  for (int d = 0; d < Dim; ++d) {
    complete.momentum.row(d) = mask.select(cons.momentum.row(d), zero);
  }
  complete.energy = mask.select(cons.energy, zero);
  const Array1d e_kin = KineticEnergy(cons.density, cons.momentum, mask);
  const Array1d e_int = cons.energy - e_kin;
  complete.pressure = mask.select(e_int / gamma_minus_1_inv, zero);
  Array1d safe_density = mask.select(complete.density, 1.0);
  FUB_ASSERT((safe_density > 0.0).all());
  complete.speed_of_sound = mask.select(
      (gamma_array_ * complete.pressure / safe_density).sqrt(), zero);
}

template <int Dim>
Complete<PerfectGas<Dim>>
PerfectGas<Dim>::CompleteFromPrim(double rho, const Array<double, Dim, 1>& v,
                                  double p) const noexcept {
  Complete q{};
  q.density = rho;
  q.momentum = rho * v;
  q.pressure = p;
  const double e_kin = KineticEnergy(q.density, q.momentum);
  q.energy = e_kin + p * gamma_minus_1_inv;
  q.speed_of_sound = std::sqrt(gamma * q.pressure / q.density);
  return q;
}

template <int Dim>
Array<double, Dim, 1> PerfectGas<Dim>::Velocity(const Complete& complete) const
    noexcept {
  Array<double, Dim, 1> u = complete.momentum / complete.density;
  return u;
}

template <int Dim>
double PerfectGas<Dim>::Machnumber(const Complete& complete) const noexcept {
  double u = Velocity(complete).matrix().norm();
  return u / complete.speed_of_sound;
}

template struct PerfectGas<1>;
template struct PerfectGas<2>;
template struct PerfectGas<3>;

void Rotate(Conservative<PerfectGas<2>>& rotated,
            const Conservative<PerfectGas<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation, const PerfectGas<2>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Rotate(Complete<PerfectGas<2>>& rotated,
            const Complete<PerfectGas<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation, const PerfectGas<2>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.pressure = state.pressure;
  rotated.speed_of_sound = state.speed_of_sound;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Rotate(Conservative<PerfectGas<3>>& rotated,
            const Conservative<PerfectGas<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation, const PerfectGas<3>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Rotate(Complete<PerfectGas<3>>& rotated,
            const Complete<PerfectGas<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation, const PerfectGas<3>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.pressure = state.pressure;
  rotated.speed_of_sound = state.speed_of_sound;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Reflect(Complete<PerfectGas<1>>& reflected,
             const Complete<PerfectGas<1>>& state,
             const Eigen::Matrix<double, 1, 1>& normal, const PerfectGas<1>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
}

void Reflect(Complete<PerfectGas<2>>& reflected,
             const Complete<PerfectGas<2>>& state,
             const Eigen::Vector2d& normal, const PerfectGas<2>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
}

void Reflect(Complete<PerfectGas<3>>& reflected,
             const Complete<PerfectGas<3>>& state,
             const Eigen::Vector3d& normal, const PerfectGas<3>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
}

} // namespace fub
