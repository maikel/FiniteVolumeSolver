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

#ifndef FUB_IDEAL_GAS_FLUX_HLL_METHOD_HPP
#define FUB_IDEAL_GAS_FLUX_HLL_METHOD_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/core/span.hpp"
#include "fub/flux_method/FluxMethod.hpp"

#include <numeric>

namespace fub {

template <typename EquationT, typename SignalSpeeds> class Hll {
public:
  using Equation = EquationT;
  using Complete = typename Equation::Complete;
  using Conservative = typename Equation::Conservative;

  static constexpr int ChunkSize = 1;

  Hll(const Equation& equation, const SignalSpeeds& signals)
      : equation_{equation}, signal_speeds_{signals} {}

  void SolveRiemannProblem(Complete& solution, const Complete& left,
                           const Complete& right, Direction dir) {
    const auto signals = signal_speeds_(equation_, left, right, dir);
    equation_.Flux(flux_left, left, dir);
    equation_.Flux(flux_right, right, dir);
    const double sL = signals[0];
    const double sR = signals[1];
    const double ds = sR - sL;
    if (0.0 < sL) {
      solution = left;
    } else if (sR < 0.0) {
      solution = right;
    } else {
      ForEachComponent<Conservative>(
          [&](double& sol, double fL, double fR, double qL, double qR) {
            sol = (sR * qR - sL * qL + fL - fR) / ds;
          },
          (Conservative&)(solution), flux_left, flux_right, left, right);
      CompleteFromCons(equation_, solution, solution);
    }
  }

  void ComputeNumericFlux(Conservative& numeric_flux,
                          span<const Complete, 2> states, Duration /* dt */,
                          double /* dx */, Direction dir) {
    const Complete& left = states[0];
    const Complete& right = states[1];

    const auto signals = signal_speeds_(equation_, left, right, dir);

    equation_.Flux(flux_left, left, dir);
    equation_.Flux(flux_right, right, dir);

    const double sL = std::min(0.0, signals[0]);
    const double sR = std::max(0.0, signals[1]);
    const double sLsR = sL * sR;
    const double ds = sR - sL;
    FUB_ASSERT(ds > 0);

    ForEachComponent<Conservative>(
        [&](double& nf, double fL, double fR, double qL, double qR) {
          nf = (sR * fL - sL * fR + sLsR * (qR - qL)) / ds;
        },
        numeric_flux, flux_left, flux_right, states[0], states[1]);
  }

  double ComputeStableDt(span<const Complete, 2> states, double dx,
                         Direction dir) {
    const auto signals = signal_speeds_(equation_, states[0], states[1], dir);
    const double max = std::accumulate(
        signals.begin(), signals.end(), 0.0,
        [](double x, double y) { return std::max(x, std::abs(y)); });
    FUB_ASSERT(max > 0.0);
    return 0.5 * dx / max;
  }

  static constexpr int GetStencilWidth() noexcept { return 1; }

  const Equation& GetEquation() const noexcept { return equation_; }

private:
  Equation equation_;
  SignalSpeeds signal_speeds_;
  Conservative flux_left{equation_};
  Conservative flux_right{equation_};
};

template <typename Equation, typename Signals>
struct HllMethod : public FluxMethod<Hll<Equation, Signals>> {
  using FluxMethod<Hll<Equation, Signals>>::FluxMethod;
};

template <typename Equation, typename Signals>
HllMethod(const Equation& eq, const Signals& signals)
    ->HllMethod<Equation, Signals>;

} // namespace fub

#endif
