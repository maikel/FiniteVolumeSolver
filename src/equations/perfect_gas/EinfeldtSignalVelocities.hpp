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

#include "fub/equations/perfect_gas/EinfeldtSignalVelocities.hpp"

namespace fub {

template <int Dim>
std::array<double, 2> EinfeldtSignalVelocities<PerfectGas<Dim>>::operator()(
    const PerfectGas<Dim>&, const Complete& left, const Complete& right,
    Direction dir) const noexcept {
  FUB_ASSERT(left.density > 0.0);
  FUB_ASSERT(right.density > 0.0);
  const double rhoL = left.density;
  const double rhoR = right.density;
  double rhoUL;
  double rhoUR;
  rhoUL = left.momentum[int(dir)];
  rhoUR = right.momentum[int(dir)];
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

template <int Dim>
std::array<Array1d, 2> EinfeldtSignalVelocities<PerfectGas<Dim>>::operator()(
    const PerfectGas<Dim>&, const CompleteArray& left,
    const CompleteArray& right, Direction dir) const noexcept {
  const Array1d rhoL = left.density;
  const Array1d rhoR = right.density;
  const Array1d rhoUL = left.momentum.row(int(dir));
  const Array1d rhoUR = right.momentum.row(int(dir));
  const Array1d aL = left.speed_of_sound;
  const Array1d aR = right.speed_of_sound;
  const Array1d sqRhoL = rhoL.sqrt();
  const Array1d sqRhoR = rhoR.sqrt();
  const Array1d uL = rhoUL / rhoL;
  const Array1d uR = rhoUR / rhoR;
  const Array1d roeU = (sqRhoL * uL + sqRhoR * uR) / (sqRhoL + sqRhoR);
  const Array1d roeA =
      ((sqRhoL * aL * aL + sqRhoR * aR * aR) / (sqRhoL + sqRhoR) +
       0.5 * (sqRhoL * sqRhoR) / ((sqRhoL + sqRhoR) * (sqRhoL + sqRhoR)) *
           (uR - uL) * (uR - uL))
          .sqrt();
  const Array1d sL1 = uL - aL;
  const Array1d sL2 = roeU - 0.5 * roeA;
  const Array1d sR1 = roeU + 0.5 * roeA;
  const Array1d sR2 = uR + aR;
  return {sL1.min(sL2), sR1.max(sR2)};
}

template <int Dim>
std::array<Array1d, 2> EinfeldtSignalVelocities<PerfectGas<Dim>>::operator()(
    const PerfectGas<Dim>&, const CompleteArray& left,
    const CompleteArray& right, const MaskArray& mask, Direction dir) const
    noexcept {
  const Array1d rhoL = left.density;
  const Array1d rhoR = right.density;
  const Array1d zero = Array1d::Zero();
  const Array1d one = Array1d::Constant(1.0);
  const Array1d rhoUL = mask.select(left.momentum.row(int(dir)), zero);
  const Array1d rhoUR = mask.select(right.momentum.row(int(dir)), zero);
  const Array1d aL = mask.select(left.speed_of_sound, zero);
  const Array1d aR = mask.select(right.speed_of_sound, zero);
  const Array1d rhoLs = mask.select(rhoL, one);
  const Array1d rhoRs = mask.select(rhoR, one);
  const Array1d sqRhoL = rhoLs.sqrt();
  const Array1d sqRhoR = rhoRs.sqrt();
  const Array1d sqRho = sqRhoL + sqRhoR;
  FUB_ASSERT((rhoLs > 0.0).all());
  FUB_ASSERT((rhoRs > 0.0).all());
  FUB_ASSERT((sqRho > 0.0).all());
  const Array1d uL = rhoUL / rhoLs;
  const Array1d uR = rhoUR / rhoRs;
  const Array1d sqRhoL_over_sqRho = sqRhoL / sqRho;
  const Array1d sqRhoR_over_sqRho = sqRhoR / sqRho;
  const Array1d roeU = sqRhoL_over_sqRho * uL + sqRhoR_over_sqRho * uR;
  const Array1d roeA =
      (sqRhoL_over_sqRho * aL * aL + sqRhoR_over_sqRho * aR * aR +
       Array1d::Constant(0.5) * sqRhoL_over_sqRho * sqRhoR_over_sqRho *
           (uR - uL) * (uR - uL))
          .sqrt();
  const Array1d sL1 = uL - aL;
  const Array1d sL2 = roeU - Array1d::Constant(0.5) * roeA;
  const Array1d sR1 = roeU + Array1d::Constant(0.5) * roeA;
  const Array1d sR2 = uR + aR;
  return {sL1.min(sL2), sR1.max(sR2)};
}

} // namespace fub