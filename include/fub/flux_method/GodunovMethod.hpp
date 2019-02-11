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

#ifndef FUB_IDEAL_GAS_FLUX_GODUNOV_METHOD_HPP
#define FUB_IDEAL_GAS_FLUX_GODUNOV_METHOD_HPP

#include "fub/Equation.hpp"

namespace fub {

template <typename Equation> class GodunovMethod {
public:
  using State = SimdifyT<typename Equation::State>;
  using Cons = SimdifyT<typename Equation::Cons>;
  template <typename T> using ConsView = ConsView<T, Equation>;
  template <typename T> using StateView = StateView<T, Equation>;

  explicit GodunovMethod(const Equation& eq)
      : equation_{eq} {}

  static constexpr int stencil_width() noexcept { return 1; }

  void ComputeNumericFlux(Cons& numeric_flux, span<const State, 2> stencil,
                          double /* dt */, double /* dx */,
                          Direction dir = Direction::X) {
    equation_.SolveRiemannProblem(riemann_solution_, stencil[0], stencil[1],
                                  dir);
    equation_.Flux(numeric_flux, riemann_solution_, dir);
  }

  double ComputeStableDt(span<const State, 2> stencil, double dx) {
    const auto signals = equation_.ComputeSignals(stencil[0], stencil[1]);
    double s = 0;
    for (int c = 0; c < signals.cols(); ++c) {
      s = std::max(s, signals.col(c).abs().maxCoeff());
    }
    return 0.5 * dx / s;
  }

  double ComputeStableDtOnSpans(StateView<const double> states, double dx,
                                Direction dir = Direction::X) {
    double dt = std::numeric_limits<double>::infinity();
    ForEachRow(
        Extents(states),
        [&](auto state) {
          for (int i = 0; i < Extents(state).extent(0) - 1; i += kChunkSize) {
            Load(span(stencil_), state, i);
            dt = std::min(dt, ComputeStableDt(span(stencil_), dx));
          }
        },
        states);
    return dt;
  }

  void ComputeNumericFluxesOnSpans(ConsView<double> fluxes,
                                   StateView<const double> states, double dt,
                                   double dx, Direction dir = Direction::X) {
    ForEachRow(
        Extents(fluxes),
        [&](auto flux, auto state) {
          for (int i = 0; i < Extents(flux).extent(0); i += kChunkSize) {
            Load(span(stencil_), state, i);
            ComputeNumericFlux(numeric_flux_, stencil_, dt, dx, dir);
            Store(flux, numeric_flux_, i);
          }
        },
        fluxes, states);
  }

private:
  Equation equation_;
  std::array<State, 2> stencil_{State(equation_), State(equation_)};
  State riemann_solution_{equation_};
  Cons numeric_flux_{equation_};
};

} // namespace fub

#endif