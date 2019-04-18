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

#include "fub/ideal_gas/flux_method/HlleGodunovMethod.hpp"

namespace fub {
namespace ideal_gas {

HlleGodunovMethod::HlleGodunovMethod(int n_speces) {
  flux_left.species.resize(kChunkSize, n_speces);
  flux_right.species.resize(kChunkSize, n_speces);
}

Array2d HlleGodunovMethod::ComputeEinfeldtSignalVelocities(
    const State& left, const State& right) noexcept {
  const Array1d sqRhoL = left.density.sqrt();
  const Array1d sqRhoR = right.density.sqrt();
  const Array1d uL = left.momentum.col(0) / left.density;
  const Array1d uR = right.momentum.col(0) / right.density;
  const Array1d aL2 = left.speed_of_sound * left.speed_of_sound;
  const Array1d aR2 = right.speed_of_sound * right.speed_of_sound;
  const Array1d roeU = (sqRhoL * uL + sqRhoR * uR) / (sqRhoL + sqRhoR);
  const Array1d roeA =
      0.5 * ((sqRhoL * aL2 + sqRhoR * aR2) / (sqRhoL + sqRhoR) +
             0.5 * (sqRhoL * sqRhoR) / ((sqRhoL + sqRhoR) * (sqRhoL + sqRhoR)) *
                 (uR - uL) * (uR - uL))
                .sqrt();
  Array2d signals;
  signals.col(0) = (uL - left.speed_of_sound).cwiseMin(roeU - roeA);
  signals.col(1) = (uR + right.speed_of_sound).cwiseMax(roeU + roeA);
  return signals;
}

void HlleGodunovMethod::ComputeNumericFlux(Conservative& flux,
                                           const State& left,
                                           const State& right) noexcept {
  const Array2d signals = ComputeEinfeldtSignalVelocities(left, right);
  const Array1d sL = signals.col(0);
  const Array1d sR = signals.col(1);
  IdealGasEquation::Flux(flux_left, left);
  IdealGasEquation::Flux(flux_right, right);
  ForEachVariable(
      [&](auto flux, auto qL, auto qR, auto fL, auto fR) {
        flux = (sL * fR - sR * fL + sL * sR * (qL - qR)) / (sR - sL);
      },
      flux, left, right, flux_left, flux_right);
}

void HlleGodunovMethod::ComputeNumericFluxesOnSpans(
    ConsSpan<double> fluxes, StateSpan<const double> states) noexcept {
  ForEachRow(
      Extents<0>(states),
      [&](auto fluxes, auto states) {
        FUB_ASSERT(states.extent(0) == fluxes.extent(0) + 1);
        for (int i = 0; i + 1 < Extents<0>(states).extent(0); i += kChunkSize) {
          Load(left, states, i);
          Load(right, states, i + 1);
          ComputeNumericFlux(flux, left, right);
          Store(fluxes, flux, i);
        }
      },
      fluxes, states);
}

} // namespace ideal_gas
} // namespace fub