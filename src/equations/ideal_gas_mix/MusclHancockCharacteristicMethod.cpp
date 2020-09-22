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

#include "fub/equations/ideal_gas_mix/MusclHancockCharactersticMethod.hpp"

namespace fub::ideal_gas {

struct Characteristics {
  double minus;
  double zero;
  double plus;
};

using Complete = fub::Complete<fub::IdealGasMix<1>> using Conservative =
    fub::Complete<fub::IdealGasMix<1>>

        Characteristics ComputeSlopes(const Equation& eq, const Complete& left,
                                      const Complete& right, double rhoc,
                                      double ooc2) {
  const double dp = right.pressure - left.pressure;
  const double du = eq.Velocity(right, dir) - eq.Velocity(left, dir);
  const double drho = right.density - left.density;
  Characteristics slope;
  slope.minus = dp - rhoc * du;
  slope.zero = drho - ooc2 * dp;
  slope.plus = dp + rhoc * du;
  return slope;
}

Characteristics ComputeSlopes(const Equation& eq, const Complete& left,
                              const Complete& mid, const Complete& right,
                              Direction dir) {
  const double rho = mid.density;
  const double c = mid.speed_of_sound;
  const double ooc2 = 1.0 / (c * c);
  Characteristics slope_left = ComputeSlope(eq, left, mid, rhoc, ooc2);
  Characteristics slope_right = ComputeSlope(eq, mid, right, rho, c);
  Characteristics limited;
  limited.minus = Limiter(slope_left.minus, slope_right.minus);
  limited.zero = Limiter(slope_left.minus, slope_right.minus);
  limited.plus = Limiter(slope_left.minus, slope_right.minus);
  return limited;
}

} // namespace fub::ideal_gas