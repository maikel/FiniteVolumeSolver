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

#ifndef FUB_SOLVER_DIMENSIONAL_SPLIT_HYPERBOLIC_TIME_INTEGRATOR_HPP
#define FUB_SOLVER_DIMENSIONAL_SPLIT_HYPERBOLIC_TIME_INTEGRATOR_HPP

#include "fub/Direction.hpp"
#include "fub/Equation.hpp"

namespace fub {

template <typename Eq> class DimensionalSplitHyperbolicTimeIntegrator {
public:
  using Equation = Eq;
  using State = SimdifyT<typename Equation::State>;
  using Cons = SimdifyT<typename Equation::Cons>;

  template <typename T> using ConsSpan = ConsView<T, Equation>;
  template <typename T> using StateSpan = StateView<T, Equation>;
  template <typename T>
  using StateStridedSpan = StateView<T, Equation, layout_stride>;

  DimensionalSplitHyperbolicTimeIntegrator(const Equation& eq) : equation{eq} {}

  void AdvanceTimeOnSpans(StateSpan<double> next,
                          StateStridedSpan<const double> prev,
                          ConsSpan<const double> fluxes, double dt, double dx,
                          Direction dir = Direction::X);

  void AdvanceInTime(State& next, const State& old, const Cons& fL,
                     const Cons& fR, const Scalar& lambda);

private:
  Equation equation;

  State next_state;
  State prev_state;
  Cons flux_left;
  Cons flux_right;
};

template <typename Equation>
void DimensionalSplitHyperbolicTimeIntegrator<Equation>::AdvanceInTime(
    State& next, const State& prev, const Cons& flux_left,
    const Cons& flux_right, const Array1d& lambda) {
  ForEachVariable(
      [lambda](auto&& next, auto prev, auto flux_left, auto flux_right) {
        next = prev + lambda * (flux_left - flux_right);
      },
      AsCons(next), AsCons(prev), flux_left, flux_right);

  equation.Reconstruct(next, AsCons(next));
}

template <typename Equation>
void DimensionalSplitHyperbolicTimeIntegrator<Equation>::AdvanceTimeOnSpans(
    StateSpan<double> next, StateStridedSpan<const double> prev,
    ConsSpan<const double> fluxes, double dt, double dx, Direction dir) {
  Array1d lambda;
  lambda.fill(dt / dx);
  ForEachRow(
      Extents(next),
      [&](auto next, auto fluxes, auto prev) {
        // Now compute the numeric flux for one side only reusing the previous
        // computation.
        for (int i = 0; i + 1 < Extents(prev).extent(0); i += kChunkSize) {
          Load(flux_left, fluxes, i);
          Load(flux_right, fluxes, i + 1);
          Load(prev_state, prev, i);
          AdvanceInTime(next_state, prev_state, flux_left, flux_right, lambda);
          Store(next, next_state, i);
        }
      },
      next, fluxes, prev);
}

} // namespace fub

#endif