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

#include "fub/equations/ShallowWater.hpp"
#include "fub/NewtonIteration.hpp"

namespace fub {

void ShallowWater::Flux(Cons& flux, const State& state, Direction dir) const
    noexcept {
  const int dir_v = int(dir);
  const Scalar velocity = state.momentum.col(dir_v);
  flux.heigth = state.momentum.col(dir_v);
  flux.momentum.col(0) = velocity * state.momentum.col(0);
  flux.momentum.col(1) = velocity * state.momentum.col(1);
  flux.momentum.col(dir_v) += gravity_ * state.heigth;
}

void ShallowWater::Reconstruct(State& state, const Cons& cons) const noexcept {
  state.heigth = cons.heigth;
  state.momentum = cons.momentum;
}

Array2d ShallowWater::ComputeSignals(const State& left, const State& right,
                                     Direction) {
  const Array1d uL = left.momentum.col(0) / left.heigth;
  const Array1d uR = right.momentum.col(0) / right.heigth;
  const Array1d sqrt_g_hL = (gravity_ * left.heigth).sqrt();
  const Array1d sqrt_g_hR = (gravity_ * right.heigth).sqrt();
  Array2d signals;
  signals.col(0) = uL - sqrt_g_hL;
  signals.col(1) = uR - sqrt_g_hR;
  return signals;
}

void ShallowWater::SolveRiemannProblem(State& state, const State& left,
                                       const State& right, Direction) {
  const Array1d uL = left.momentum.col(0) / left.heigth;
  const Array1d hL = left.heigth;
  const Array1d g = gravity_;
  const Array1d sqrt_g_hL = (g * hL).sqrt();
  auto f_left = [&](const Array1d& h) {
    const auto rarefaction = uL + 2.0 * (sqrt_g_hL - (g * h).sqrt());
    const auto shock = uL - (h - hL) * (0.5 * g * (1.0 / h + 1.0 / hL)).sqrt();
    Array1d result = (h < hL).select(rarefaction, shock);
    return result;
  };
  auto df_left = [&](const Array1d& h) {
    const auto rarefaction = -g / (g * h).sqrt();
    const auto shock = 0.25 * g * (h - hL) /
                           (h * h * (0.5 * g * (1.0 / h + 1.0 / hL)).sqrt()) -
                       (0.5 * g * (1.0 / h + 1.0 / hL)).sqrt();
    Array1d result = (h < hL).select(rarefaction, shock);
    return result;
  };
  const Array1d uR = right.momentum.col(0) / right.heigth;
  const Array1d hR = right.heigth;
  const Array1d sqrt_g_hR = (g * hR).sqrt();
  auto f_right = [&](const Array1d& h) {
    const auto rarefaction = uR - 2.0 * (sqrt_g_hR - (g * h).sqrt());
    const auto shock = uR + (h - hR) * (0.5 * g * (1.0 / h + 1.0 / hR)).sqrt();
    Array1d result = (h < hR).select(rarefaction, shock);
    return result;
  };
  auto df_right = [&](const Array1d& h) {
    const auto rarefaction = g / (g * h).sqrt();
    const auto shock =
        (0.5 * g * (1.0 / h + 1.0 / hR)).sqrt() -
        0.25 * g * (h - hR) / (h * h * (0.5 * g * (1.0 / h + 1.0 / hR)).sqrt());
    Array1d result = (h < hR).select(rarefaction, shock);
    return result;
  };
  auto f = [&](const Array1d& h) {
    Array1d result = f_right(h) - f_left(h);
    return result;
  };
  auto df = [&](const Array1d& h) {
    Array1d result = df_right(h) - df_left(h);
    return result;
  };
  Array1d h_star = NewtonIteration(f, df, 0.5 * (left.heigth + right.heigth));
  Array1d hu_star = f_left(h_star);
  const Array1d sL = uL - sqrt_g_hL;
  const Array1d sR = uR + sqrt_g_hR;
  h_star = (sL > 0.0).select(left.heigth, h_star);
  hu_star = (sL > 0.0).select(left.momentum.col(0), hu_star);
  h_star = (sR < 0.0).select(right.heigth, h_star);
  hu_star = (sR < 0.0).select(right.momentum.col(0), hu_star);
  state.heigth = h_star;
  state.momentum.col(0) = hu_star;
}

} // namespace fub