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

#ifndef FUB_EQUATIONS_PERFECT_GAS_EXACT_RIEMANN_SOLVER_HPP
#define FUB_EQUATIONS_PERFECT_GAS_EXACT_RIEMANN_SOLVER_HPP

#include "fub/equations/PerfectGas.hpp"
#include "fub/ExactRiemannSolver.hpp"

namespace fub {

template <int Dim> class ExactRiemannSolver<PerfectGas<Dim>> {
public:
  using Complete = typename PerfectGas<Dim>::Complete;
  using CompleteArray = typename PerfectGas<Dim>::CompleteArray;

  explicit ExactRiemannSolver(const PerfectGas<Dim>& equation)
      : equation_{equation} {}

  /// Returns either left or right, depending on the upwind velocity.
  void SolveRiemannProblem(Complete& state, const Complete& left,
                           const Complete& right, Direction dir);

  void SolveRiemannProblem(CompleteArray& state, const CompleteArray& left,
                           const CompleteArray& right, Direction dir);

  void SolveRiemannProblem(CompleteArray& state, const CompleteArray& left,
                           const CompleteArray& right, MaskArray mask,
                           Direction dir);

  /// Returns the upwind velocity in the specified direction.
  std::array<double, 2> ComputeSignals(const Complete&, const Complete&,
                                       Direction dir);

  std::array<Array1d, 2> ComputeSignals(const CompleteArray&,
                                        const CompleteArray&, Direction dir);

  std::array<double, 2> ComputeMiddleState(const Complete& left,
                                           const Complete& right,
                                           Direction dir);

private:
  PerfectGas<Dim> equation_;
};

extern template class ExactRiemannSolver<PerfectGas<1>>;
extern template class ExactRiemannSolver<PerfectGas<2>>;
extern template class ExactRiemannSolver<PerfectGas<3>>;

}

#endif