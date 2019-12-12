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

void ShallowWater::Flux(Conservative& flux, const Complete& state,
                        Direction dir) const noexcept {
  const int d = int(dir);
  const double velocity = state.momentum[d];
  flux.height = state.momentum[d];
  flux.momentum[0] = velocity * state.momentum[0];
  flux.momentum[1] = velocity * state.momentum[1];
  flux.momentum[d] += 0.5 * gravity_ * state.height * state.height;
}

void ShallowWater::Flux(ConservativeArray& flux, const CompleteArray& state,
                        Direction dir) const noexcept {
  const int d = int(dir);
  const Array1d velocity = state.momentum.row(d);
  flux.height = state.momentum.row(d);
  flux.momentum.row(0) = velocity * state.momentum.row(0);
  flux.momentum.row(1) = velocity * state.momentum.row(1);
  flux.momentum.row(d) += 0.5 * gravity_ * state.height * state.height;
}

ShallowWater::Complete ExactRiemannSolver<ShallowWater>::ComputeMiddleState(
    const Complete& left, const Complete& right, Direction dir) {
  int d0 = static_cast<int>(dir);
  const double uL = left.momentum[d0] / left.height;
  const double hL = left.height;
  const double g = equation_.gravity_;
  const double sqrt_g_hL = std::sqrt(g * hL);
  auto f_left = [&](double h) {
    if (h < hL) {
      return uL + 2.0 * (sqrt_g_hL - std::sqrt(g * h));
    } else {
      return uL - (h - hL) * std::sqrt(0.5 * g * (1.0 / h + 1.0 / hL));
    }
  };
  auto df_left = [&](double h) {
    if (h < hL) {
      return -g / std::sqrt(g * h);
    } else {
      return 0.25 * g * (h - hL) /
                 (h * h * std::sqrt(0.5 * g * (1.0 / h + 1.0 / hL))) -
             std::sqrt(0.5 * g * (1.0 / h + 1.0 / hL));
    }
  };
  const double uR = right.momentum[d0] / right.height;
  const double hR = right.height;
  const double sqrt_g_hR = std::sqrt(g * hR);
  auto f_right = [&](double h) {
    if (h < hR) {
      return uR - 2.0 * (sqrt_g_hR - std::sqrt(g * h));
    } else {
      return uR + (h - hR) * std::sqrt(0.5 * g * (1.0 / h + 1.0 / hR));
    }
  };
  auto df_right = [&](double h) {
    if (h < hR) {
      return g / std::sqrt(g * h);
    } else {
      return std::sqrt(0.5 * g * (1.0 / h + 1.0 / hR)) -
             0.25 * g * (h - hR) /
                 (h * h * std::sqrt(0.5 * g * (1.0 / h + 1.0 / hR)));
    }
  };
  auto f = [&](double h) { return f_right(h) - f_left(h); };
  auto df = [&](double h) { return df_right(h) - df_left(h); };

  const double hM = NewtonIteration(f, df, 0.5 * (left.height + right.height));
  // FUB_ASSERT(f_left(hM) == f_right(hM));
  const double uM = f_left(hM);

  Complete middle;
  middle.height = hM;
  middle.momentum[d0] = hM * uM;
  const int d1 = (d0 + 1) % 2;
  middle.momentum[d1] = 0;
  return middle;
}

std::array<double, 2> ExactRiemannSolver<ShallowWater>::ComputeSignals(
    const Complete& left, const Complete& right, Direction dir) {
  const Complete middle = ComputeMiddleState(left, right, dir);
  std::array<double, 2> signals;
  int d0 = static_cast<int>(dir);
  // Compute Left fastest signal
  if (middle.height <= left.height) {
    FUB_ASSERT(0 < left.height);
    const double uL = left.momentum[d0] / left.height;
    const double sqrt_g_hL = std::sqrt(equation_.gravity_ * left.height);
    signals[0] = uL - sqrt_g_hL;
  } else {
    signals[0] = (middle.momentum[d0] - left.momentum[d0]) /
                 (middle.height - left.height);
  }
  // Compute right fastest signal
  if (middle.height <= right.height) {
    FUB_ASSERT(0 < right.height);
    const double uR = right.momentum[d0] / right.height;
    const double sqrt_g_hR = std::sqrt(equation_.gravity_ * right.height);
    signals[1] = uR + sqrt_g_hR;
  } else {
    signals[1] = (right.momentum[d0] - middle.momentum[d0]) /
                 (right.height - middle.height);
  }
  return signals;
}

std::array<Array1d, 2> ExactRiemannSolver<ShallowWater>::ComputeSignals(
    const CompleteArray& left, const CompleteArray& right, Direction dir) {
  Complete left_i{};
  Complete right_i{};
  std::array<Array1d, 2> signals;
  for (int i = 0; i < kDefaultChunkSize; ++i) {
    left_i.height = left.height(i);
    left_i.momentum = left.momentum.col(i);
    right_i.height = right.height(i);
    right_i.momentum = right.momentum.col(i);
    std::array<double, 2> signals_i = ComputeSignals(left_i, right_i, dir);
    signals[0][i] = signals_i[0];
    signals[1][i] = signals_i[1];
  }
  return signals;
}

void ExactRiemannSolver<ShallowWater>::SolveRiemannProblem(
    Complete& state, const Complete& left, const Complete& right,
    Direction dir) {
  FUB_ASSERT(0 < left.height);
  FUB_ASSERT(0 < right.height);
  Complete middle = ComputeMiddleState(left, right, dir);
  FUB_ASSERT(0 < middle.height);
  int d0 = static_cast<int>(dir);
  int d1 = (d0 + 1) % 2;
  const double g = equation_.gravity_;
  const double uL = left.momentum[d0] / left.height;
  const double vL = left.momentum[d1] / left.height;
  const double hL = left.height;
  const double uR = right.momentum[d0] / right.height;
  const double vR = right.momentum[d1] / right.height;
  const double hR = right.height;
  const double hM = middle.height;
  const double uM = middle.momentum[d0] / middle.height;
  const double sqrt_g_hM = std::sqrt(g * hM);
  const double sqrt_g_hL = std::sqrt(g * hL);
  const double sqrt_g_hR = std::sqrt(g * hR);

  // (II) Determine Structure of Solution

  // Left side is interesting
  if (0 <= uM) {
    middle.momentum[d1] = hM * vL;
    // left side is Shock
    if (hM > hL) {
      const double sL = (hM * uM - hL * uL) / (hM - hL);
      if (0 < sL) {
        state = left;
      } else {
        state = middle;
      }
      // left side is rarefaction
    } else {
      if (0 < uL - sqrt_g_hL) {
        state = left;
      } else if (uM - sqrt_g_hM < 0) {
        state = middle;
      } else {
        const double a = uL + 2.0 * sqrt_g_hL;
        const double h_rarefaction = a * a / (9.0 * g);
        const double u_rarefaction =
            uL + 2.0 * (sqrt_g_hL - std::sqrt(g * h_rarefaction));
        state.height = h_rarefaction;
        state.momentum[d0] = h_rarefaction * u_rarefaction;
        state.momentum[d1] = h_rarefaction * vL;
      }
    }
    // Right side is interesting
  } else {
    middle.momentum[d1] = hM * vR;
    // right side is shock
    if (hM > hR) {
      const double sR = (hR * uR - hM * uM) / (hR - hM);
      if (sR < 0) {
        state = right;
      } else {
        state = middle;
      }
      // right side is rarefaction
    } else {
      if (uR + sqrt_g_hR < 0) {
        state = right;
      } else if (0 < uM + sqrt_g_hM) {
        state = middle;
      } else {
        const double a = uR - 2.0 * sqrt_g_hR;
        const double h_rarefaction = a * a / (9.0 * g);
        const double u_rarefaction =
            uR - 2.0 * (sqrt_g_hR - std::sqrt(g * h_rarefaction));
        state.height = h_rarefaction;
        state.momentum[d0] = h_rarefaction * u_rarefaction;
        state.momentum[d1] = h_rarefaction * vR;
      }
    }
  }
}

void ExactRiemannSolver<ShallowWater>::SolveRiemannProblem(
    CompleteArray& state, const CompleteArray& left, const CompleteArray& right,
    Direction dir) {
  Complete left_i{};
  Complete right_i{};
  Complete solution_i{};
  std::array<Array1d, 2> signals;
  for (int i = 0; i < kDefaultChunkSize; ++i) {
    left_i.height = left.height(i);
    left_i.momentum = left.momentum.col(i);
    right_i.height = right.height(i);
    right_i.momentum = right.momentum.col(i);
    SolveRiemannProblem(solution_i, left_i, right_i, dir);
    state.height(i) = solution_i.height;
    state.momentum.col(i) = solution_i.momentum;
  }
}

std::array<double, 2> ShallowWaterSignalVelocities::
operator()(const ShallowWater& equation, const Complete& left,
           const Complete& right, Direction dir) {
  const int d = static_cast<int>(dir);
  const double vL = left.momentum[d] / left.height;
  const double sqrt_g_hL = std::sqrt(equation.gravity_ * left.height);
  const double vR = right.momentum[d] / right.height;
  const double sqrt_g_hR = std::sqrt(equation.gravity_ * right.height);
  const double sL = std::min({vL - sqrt_g_hL, vR - sqrt_g_hR, 0.0});
  const double sR = std::max({vL + sqrt_g_hL, vR + sqrt_g_hR, 0.0});
  return {sL, sR};
}

std::array<Array1d, 2> ShallowWaterSignalVelocities::
operator()(const ShallowWater& equation, const CompleteArray& left,
           const CompleteArray& right, Direction dir) {
  const int d = static_cast<int>(dir);
  const Array1d vL = left.momentum.row(d) / left.height;
  const Array1d sqrt_g_hL = (equation.gravity_ * left.height).sqrt();
  const Array1d vR = right.momentum.row(d) / right.height;
  const Array1d sqrt_g_hR = (equation.gravity_ * right.height).sqrt();
  const Array1d sL = (vL - sqrt_g_hL).min(vR - sqrt_g_hR).min(Array1d::Zero());
  const Array1d sR = (vL + sqrt_g_hL).max(vR + sqrt_g_hR).max(Array1d::Zero());
  return {sL, sR};
}

} // namespace fub
