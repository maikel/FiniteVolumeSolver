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

#include "fub/EinfeldtSignalVelocities.hpp"
namespace fub {

template <typename Equation>
std::array<double, 2> EinfeldtSignalVelocities<Equation>::
operator()(const Equation& equation, const Complete& left,
           const Complete& right, Direction dir) const noexcept {
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
  const double oosqRhoSum = 1.0 / (sqRhoL + sqRhoR);
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  auto roeAvg = [=](double xL, double xR) {
    return (xL * sqRhoL + xR * sqRhoR) * oosqRhoSum;
  };
  const double roeU = roeAvg(uL, uR);
  const double rhoEL = left.energy;
  const double rhoER = right.energy;
  const double pL = left.pressure;
  const double pR = right.pressure;
  const double hL = (rhoEL + pL) / rhoL;
  const double hR = (rhoER + pR) / rhoR;
  const double roeH = roeAvg(hL, hR);
  const double roeA2 = (equation.gamma - 1.0) * (roeH - 0.5 * roeU * roeU);
  const double roeA = std::sqrt(roeA2);
  const double sL1 = uL - aL;
  const double sL2 = roeU - roeA;
  const double sR1 = roeU + roeA;
  const double sR2 = uR + aR;
  return {std::min(sL1, sL2), std::max(sR1, sR2)};
}

template <typename Equation>
std::array<Array1d, 2> EinfeldtSignalVelocities<Equation>::
operator()(const Equation& equation, const CompleteArray& left,
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
  const Array1d oosqRhoSum = Array1d::Constant(1.0) / (sqRhoL + sqRhoR);
  auto roeAvg = [=](Array1d xL, Array1d xR) {
    return (xL * sqRhoL + xR * sqRhoR) * oosqRhoSum;
  };
  const Array1d roeU = roeAvg(uL, uR);
  const Array1d rhoEL = left.energy;
  const Array1d rhoER = right.energy;
  const Array1d pL = left.pressure;
  const Array1d pR = right.pressure;
  const Array1d hL = (rhoEL + pL) / rhoL;
  const Array1d hR = (rhoER + pR) / rhoR;
  const Array1d roeH = roeAvg(hL, hR);
  const Array1d roeA2 = (equation.gamma_array_ - Array1d::Constant(1.0)) *
                        (roeH - Array1d::Constant(0.5) * roeU * roeU);
  const Array1d roeA = roeA2.sqrt();
  const Array1d sL1 = uL - aL;
  const Array1d sL2 = roeU - roeA;
  const Array1d sR1 = roeU + roeA;
  const Array1d sR2 = uR + aR;
  return {sL1.min(sL2), sR1.max(sR2)};
}

template <typename Equation>
std::array<Array1d, 2> EinfeldtSignalVelocities<Equation>::
operator()(const Equation& equation, const CompleteArray& left,
           const CompleteArray& right, const MaskArray& mask,
           Direction dir) const noexcept {
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
  auto roeAvg = [=](Array1d xL, Array1d xR) {
    return (xL * sqRhoL_over_sqRho + xR * sqRhoR_over_sqRho);
  };
  const Array1d roeU = roeAvg(uL, uR);
  const Array1d rhoEL = mask.select(left.energy, zero);
  const Array1d rhoER = mask.select(right.energy, zero);
  const Array1d pL = mask.select(left.pressure, zero);
  const Array1d pR = mask.select(right.pressure, zero);
  const Array1d hL = (rhoEL + pL) / rhoL;
  const Array1d hR = (rhoER + pR) / rhoR;
  const Array1d roeH = roeAvg(hL, hR);
  const Array1d roeA2 = (equation.gamma_array_ - Array1d::Constant(1.0)) *
                        (roeH - Array1d::Constant(0.5) * roeU * roeU);
  const Array1d roeA = roeA2.sqrt();
  const Array1d sL1 = uL - aL;
  const Array1d sL2 = roeU - roeA;
  const Array1d sR1 = roeU + roeA;
  const Array1d sR2 = uR + aR;
  return {sL1.min(sL2), sR1.max(sR2)};
}

} // namespace fub