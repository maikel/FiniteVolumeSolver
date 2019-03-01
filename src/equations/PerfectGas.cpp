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

#include <numeric>

namespace fub {

namespace {
template <int Rank>
void Flux(const PerfectGas<Rank>& equation, Cons<PerfectGas<Rank>>& flux,
          const Complete<PerfectGas<Rank>>& state, Direction dir) noexcept {
  if constexpr (Rank == 1) {
    const double velocity = state.momentum / state.density;
    flux.density = state.momentum;
    flux.momentum = velocity * state.momentum + state.pressure;
    flux.energy = velocity * (state.energy + state.pressure);
  } else {
    const int d0 = static_cast<int>(dir);
    const double velocity = state.momentum[d0] / state.density;
    flux.density = state.momentum[d0];
    for (int d = 0; d < Rank; ++d) {
      flux.momentum[d] = velocity * state.momentum[d];
    }
    flux.momentum[d0] += state.pressure;
    flux.energy = velocity * (state.energy + state.pressure);
  }
}

template <int Rank>
void Reconstruct(const PerfectGas<Rank>& equation,
                 Complete<PerfectGas<Rank>>& complete,
                 const Cons<PerfectGas<Rank>>& cons) noexcept {
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  if constexpr (Rank == 1) {
    complete.pressure =
        equation.gamma *
        (cons.energy - 0.5 * cons.momentum * cons.momentum / cons.density);
  } else {
    const double E_kin =
        0.5 *
        std::inner_product(&cons.momentum[0], &cons.momentum[0] + Rank,
                           &cons.momentum[0], 0.0) /
        cons.density;
    const double rho_e_internal = complete.energy - E_kin;
    complete.pressure = (equation.gamma - 1) * rho_e_internal;
  }
  complete.speed_of_sound =
      std::sqrt(equation.gamma * complete.pressure / complete.density);
}
} // namespace

void PerfectGas<1>::Flux(Cons& flux, const Complete& state, Direction dir) const
    noexcept {
  ::fub::Flux(*this, flux, state, dir);
}

void PerfectGas<2>::Flux(Cons& flux, const Complete& state, Direction dir) const
    noexcept {
  ::fub::Flux(*this, flux, state, dir);
}

void PerfectGas<3>::Flux(Cons& flux, const Complete& state, Direction dir) const
    noexcept {
  ::fub::Flux(*this, flux, state, dir);
}

void PerfectGas<1>::Reconstruct(Complete& complete, const Cons& cons) const
    noexcept {
  ::fub::Reconstruct(*this, complete, cons);
}

void PerfectGas<2>::Reconstruct(Complete& complete, const Cons& cons) const
    noexcept {
  ::fub::Reconstruct(*this, complete, cons);
}

void PerfectGas<3>::Reconstruct(Complete& complete, const Cons& cons) const
    noexcept {
  ::fub::Reconstruct(*this, complete, cons);
}

namespace {
template <int Rank>
std::array<double, 2>
ComputeHlleSignalVelocities(const Complete<PerfectGas<Rank>>& left,
                            const Complete<PerfectGas<Rank>>& right,
                            Direction dir) noexcept {
  const double rhoL = left.density;
  const double rhoR = right.density;
  double rhoUL;
  double rhoUR;
  if constexpr (Rank == 1) {
    rhoUL = left.momentum;
    rhoUR = right.momentum;
  } else {
    rhoUL = left.momentum[int(dir)];
    rhoUR = right.momentum[int(dir)];
  }
  const double aL = left.speed_of_sound;
  const double aR = right.speed_of_sound;
  const double sqRhoL = std::sqrt(rhoL);
  const double sqRhoR = std::sqrt(rhoR);
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  const double roeU = (sqRhoL * uL + sqRhoR * uR) / (sqRhoL + sqRhoR);
  const double roeA = std::sqrt(
      (sqRhoL * aL * aL + sqRhoR * aR * aR) / (sqRhoL + sqRhoR) +
      0.5 * (sqRhoL * sqRhoR) / ((sqRhoL + sqRhoR) * (sqRhoL + sqRhoR)) *
          (uR - uL) * (uR - uL));
  const double sL1 = uL - aL;
  const double sL2 = roeU - 0.5 * roeA;
  const double sR1 = roeU + 0.5 * roeA;
  const double sR2 = uR + aR;
  return {std::min(sL1, sL2), std::max(sR1, sR2)};
}
} // namespace

std::array<double, 2> EinfeldtSignalVelocities<PerfectGas<1>>::ComputeSignals(
    const Complete& left, const Complete& right, Direction dir) const {
  return ComputeHlleSignalVelocities(left, right, dir);
}

std::array<double, 2> EinfeldtSignalVelocities<PerfectGas<2>>::ComputeSignals(
    const Complete& left, const Complete& right, Direction dir) const {
  return ComputeHlleSignalVelocities(left, right, dir);
}

std::array<double, 2> EinfeldtSignalVelocities<PerfectGas<3>>::ComputeSignals(
    const Complete& left, const Complete& right, Direction dir) const {
  return ComputeHlleSignalVelocities(left, right, dir);
}

} // namespace fub