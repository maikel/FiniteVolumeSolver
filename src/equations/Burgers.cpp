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

#include "fub/equations/Burgers.hpp"

namespace fub {

Burgers1d::Conservative Burgers1d::Flux(Complete state, Direction) const
    noexcept {
  Conservative flux{};
  flux.u = 0.5 * state.u * state.u;
  return flux;
}

Burgers1d::ConservativeArray Burgers1d::Flux(CompleteArray state,
                                             Direction) const noexcept {
  ConservativeArray flux{};
  flux.u = 0.5 * state.u * state.u;
  return flux;
}

void ExactRiemannSolver<Burgers1d>::SolveRiemannProblem(
    Complete& riemann_solution, const Complete& left, const Complete& right,
    Direction) const {
  if (right.u < left.u) {
    const double s = 0.5 * (left.u + right.u);
    riemann_solution.u = s < 0 ? right.u : left.u;
  } else {
    if (left.u > 0) {
      riemann_solution.u = left.u;
    } else if (right.u < 0) {
      riemann_solution.u = right.u;
    } else {
      riemann_solution.u = 0.0;
    }
  }
}

void ExactRiemannSolver<Burgers1d>::SolveRiemannProblem(
    CompleteArray& riemann_solution, const CompleteArray& left,
    const CompleteArray& right, Direction) const {
  const Array1d s = 0.5 * (left.u + right.u);
  const Array1d by_s = (s < 0.0).select(right.u, left.u);
  riemann_solution.u = (left.u > 0.0).select(left.u, Array1d::Zero());
  riemann_solution.u = (right.u < 0.0).select(right.u, riemann_solution.u);
  riemann_solution.u = (right.u < left.u).select(by_s, riemann_solution.u);
}

std::array<double, 1> ExactRiemannSolver<Burgers1d>::ComputeSignals(
    const Complete& left, const Complete& right, Direction) const {
  if (right.u < left.u) {
    const double s = 0.5 * (left.u + right.u);
    return {s};
  }
  return {std::max(std::abs(left.u), std::abs(right.u))};
}

std::array<Array1d, 1> ExactRiemannSolver<Burgers1d>::ComputeSignals(
    const CompleteArray& left, const CompleteArray& right, Direction) const {
  const Array1d s = 0.5 * (left.u + right.u);
  const Array1d max = left.u.abs().max(right.u.abs());
  return {(right.u < left.u).select(s, max)};
}

} // namespace fub