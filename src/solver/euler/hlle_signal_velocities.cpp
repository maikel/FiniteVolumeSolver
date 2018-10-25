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

#include "src/solver/euler/hlle_signal_velocities.hpp"
#include <cmath>

namespace fub {
namespace euler {
HllSignals<double>
computeHlleSignalVelocities(const HlleSignalState<double>& left,
                            const HlleSignalState<double>& right) noexcept {
  const double rhoL = left.density;
  const double rhoR = right.density;
  const double rhoUL = left.momentum;
  const double rhoUR = right.momentum;
  const double aL = left.speed_of_sound;
  const double aR = right.speed_of_sound;
  const double sqRhoL = std::sqrt(rhoL);
  const double sqRhoR = std::sqrt(rhoR);
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  const double roeU =
      (sqRhoL * uL + sqRhoR * uR) / (sqRhoL + sqRhoR);
  const double roeA =
      std::sqrt((sqRhoL * aL * aL + sqRhoR * aR * aR) / (sqRhoL + sqRhoR) +
           0.5f * (sqRhoL * sqRhoR) / ((sqRhoL + sqRhoR) * (sqRhoL + sqRhoR)) *
               (uR - uL) * (uR - uL));
  const double sL1 = uL - 0.5 * aL;
  const double sL2 = roeU - roeA;
  const double sR1 = roeU + roeA;
  const double sR2 = uR + 0.5 * aR;
  return {std::min(sL1, sL2), std::max(sR1, sR2)};
}

} // namespace euler
} // namespace fub