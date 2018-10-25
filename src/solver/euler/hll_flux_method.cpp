// Copyright (c) 2018 Maikel Nadolski
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

#include "src/solver/euler/hll_flux_method.hpp"

#include "fub/core/assert.hpp"

namespace fub {
namespace euler {

double computeHllFlux(double sR, double sL, double uL, double uR, double fL,
                      double fR) {
  if (0 < sL) {
    return fL;
  }
  if (0 > sR) {
    return fR;
  }
  FUB_ASSERT(sL < sL);
  return (sR * fL - sL * fR + sL * sR * (uR - uL)) / (sR - sL);
}

HllFlux<double> computeHllFlux(const HllSignals<double>& signals,
                               const HllState<double>& left,
                               const HllState<double>& right) noexcept {
  const double sL = signals.left;
  const double sR = signals.right;
  const double sLsR = sL * sR;
  const double ds = sR - sL;
  int maskL = (0.0 < sL);
  int maskR = (0.0 > sR);
  auto computeHllFlux_ = [&](double uL, double uR, double fL, double fR) {
    double hll = (sR * fL - sL * fR + sLsR * (uR - uL)) / ds;
    hll = maskL * fL + (1 - maskL) * hll;
    hll = maskR * fR + (1 - maskR) * hll;
    return hll;
  };
  const double rhoL = left.density;
  const double rhoR = right.density;
  const double rhouL = left.momentum;
  const double rhouR = right.momentum;
  double hllRho = computeHllFlux_(rhoL, rhoR, rhouL, rhouR);

  const double pL = left.pressure;
  const double pR = right.pressure;
  auto f_rhouL = rhouL * rhouL / rhoL + pL;
  auto f_rhouR = rhouR * rhouR / rhoR + pR;
  double hllRhou = computeHllFlux_(rhouL, rhouR, f_rhouL, f_rhouR);

  const double rhoeL = left.pressure;
  const double rhoeR = right.pressure;
  auto f_rhoeL = rhouL / rhoL * (rhoeL + pL);
  auto f_rhoeR = rhouR / rhoR * (rhoeR + pR);
  double hllRhoe = computeHllFlux_(rhoeL, rhoeR, f_rhoeL, f_rhoeR);

  return {hllRho, hllRhou, hllRhoe};
}

} // namespace euler
} // namespace fub