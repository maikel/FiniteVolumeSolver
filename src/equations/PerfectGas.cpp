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

template <int Dim>
void PerfectGas<Dim>::Flux(Conservative& flux, const Complete& state,
                           Direction dir) const noexcept {
  if constexpr (Dim == 1) {
    const double velocity = state.momentum / state.density;
    flux.density = state.momentum;
    flux.momentum = velocity * state.momentum + state.pressure;
    flux.energy = velocity * (state.energy + state.pressure);
  } else {
    const int d0 = static_cast<int>(dir);
    const double velocity = state.momentum[d0] / state.density;
    flux.density = state.momentum[d0];
    for (int d = 0; d < Dim; ++d) {
      flux.momentum[d] = velocity * state.momentum[d];
    }
    flux.momentum[d0] += state.pressure;
    flux.energy = velocity * (state.energy + state.pressure);
  }
}

template struct PerfectGas<1>;
template struct PerfectGas<2>;
template struct PerfectGas<3>;

template <int Dim>
void CompleteFromConsImpl<PerfectGas<Dim>>::apply(
    const PerfectGas<Dim>& equation, Complete<PerfectGas<Dim>>& complete,
    const Conservative<PerfectGas<Dim>>& cons) {
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  if constexpr (Dim == 1) {
    complete.pressure =
        (equation.gamma - 1.0) *
        (cons.energy - 0.5 * cons.momentum * cons.momentum / cons.density);
  } else {
    const double E_kin =
        0.5 * cons.momentum.matrix().squaredNorm() / cons.density;
    const double rho_e_internal = cons.energy - E_kin;
    complete.pressure = (equation.gamma - 1.0) * rho_e_internal;
  }
  complete.speed_of_sound =
      std::sqrt(equation.gamma * complete.pressure / complete.density);
}

template <int Dim>
void CompleteFromConsImpl<PerfectGas<Dim>>::apply(
    const PerfectGas<Dim>& equation, Complete<PerfectGas<Dim>>& complete,
    const Complete<PerfectGas<Dim>>& cons) {
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  if constexpr (Dim == 1) {
    complete.pressure =
        (equation.gamma - 1.0) *
        (cons.energy - 0.5 * cons.momentum * cons.momentum / cons.density);
  } else {
    const double E_kin =
        0.5 * cons.momentum.matrix().squaredNorm() / cons.density;
    const double rho_e_internal = cons.energy - E_kin;
    complete.pressure = (equation.gamma - 1.0) * rho_e_internal;
  }
  complete.speed_of_sound =
      std::sqrt(equation.gamma * complete.pressure / complete.density);
}

template struct CompleteFromConsImpl<PerfectGas<1>>;
template struct CompleteFromConsImpl<PerfectGas<2>>;
template struct CompleteFromConsImpl<PerfectGas<3>>;

template <int Dim>
std::array<double, 2> EinfeldtSignalVelocitiesImpl<PerfectGas<Dim>>::apply(
    const PerfectGas<Dim>&, const Complete& left, const Complete& right,
    Direction dir) noexcept {
  const double rhoL = left.density;
  const double rhoR = right.density;
  double rhoUL;
  double rhoUR;
  if constexpr (Dim == 1) {
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

template struct EinfeldtSignalVelocitiesImpl<PerfectGas<1>>;
template struct EinfeldtSignalVelocitiesImpl<PerfectGas<2>>;
template struct EinfeldtSignalVelocitiesImpl<PerfectGas<3>>;

} // namespace fub