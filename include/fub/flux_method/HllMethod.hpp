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
#include "fub/ForEach.hpp"
#include "fub/core/span.hpp"
#include "fub/flux_method/FluxMethod.hpp"

#include <numeric>

namespace fub {

template <typename Eq, typename SignalSpeeds> class HllBase {
public:
  using Equation = Eq;
  using Complete = typename Equation::Complete;
  using Conservative = typename Equation::Conservative;

  HllBase(const Equation& eq, const SignalSpeeds& signals)
      : equation_{eq}, signal_speeds_{signals} {}

  /// \brief Returns the stencil width of this method
  static constexpr int GetStencilWidth() noexcept { return 1; }

  /// @{
  /// \brief Returns the strategy object which computes the signal speeds for
  /// this method.
  const SignalSpeeds& GetSignalSpeeds() const noexcept {
    return signal_speeds_;
  }
  SignalSpeeds& GetSignalSpeeds() noexcept { return signal_speeds_; }
  /// @}

  /// @{
  /// \brief Returns the underlying equations object.
  const Equation& GetEquation() const noexcept { return equation_; }
  Equation& GetEquation() noexcept { return equation_; }
  /// @}

  /// \brief Computes an approximate solution to the rieman problem defined by
  /// left and right states.
  ///
  /// \param[out] solution  The solution to the riemann problem will be stored
  ///                       here.
  /// \param[in] left  The left state of the riemann problem
  /// \param[in] right  The right state of the riemann problem.
  /// \param[in] dir  The direction in which the riemann problem will be solved.
  void SolveRiemannProblem(Complete& solution, const Complete& left,
                           const Complete& right, Direction dir);

  void ComputeNumericFlux(Conservative& numeric_flux,
                          span<const Complete, 2> states, Duration dt,
                          double dx, Direction dir);

  double ComputeStableDt(span<const Complete, 2> states, double dx,
                         Direction dir);

private:
  Equation equation_;
  SignalSpeeds signal_speeds_;
  Conservative flux_left_{equation_};
  Conservative flux_right_{equation_};
};

template <typename Eq, typename SignalSpeeds, bool Simdify> class HllArrayBase;

template <typename Equation, typename SignalSpeeds>
class HllArrayBase<Equation, SignalSpeeds, true>
    : public HllBase<Equation, SignalSpeeds> {
public:
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using Base = HllBase<Equation, SignalSpeeds>;

  using Base::HllBase;

  using Base::GetEquation;
  using Base::GetSignalSpeeds;

  using Base::SolveRiemannProblem;
  void SolveRiemannProblem(CompleteArray& solution, const CompleteArray& left,
                           const CompleteArray& right, Direction dir);

  using Base::ComputeNumericFlux;
  void ComputeNumericFlux(ConservativeArray& numeric_flux,
                          span<const CompleteArray, 2> states, Duration dt,
                          double dx, Direction dir);

  void ComputeNumericFlux(ConservativeArray& numeric_flux,
                          Array1d face_fraction,
                          span<const CompleteArray, 2> states,
                          span<const Array1d, 2> volume_fraction, Duration dt,
                          double dx, Direction dir);

  using Base::ComputeStableDt;
  Array1d ComputeStableDt(span<const CompleteArray, 2> states, double dx,
                          Direction dir);

private:
  ConservativeArray flux_left_array_{Base::GetEquation()};
  ConservativeArray flux_right_array_{Base::GetEquation()};
};

template <typename Equation, typename SignalSpeeds>
class HllArrayBase<Equation, SignalSpeeds, false>
    : public HllBase<Equation, SignalSpeeds> {
  using HllBase<Equation, SignalSpeeds>::HllBase;
};

template <typename Equation, typename SignalSpeeds>
class Hll : public HllArrayBase<Equation, SignalSpeeds,
                                HasVectorizedFlux<Equation>::value> {
  using Base =
      HllArrayBase<Equation, SignalSpeeds, HasVectorizedFlux<Equation>::value>;
  using Base::Base;
};

template <typename Equation, typename Signals>
Hll(const Equation& eq, const Signals& signals)->Hll<Equation, Signals>;

template <typename Equation, typename Signals>
struct HllMethod : public FluxMethod<Hll<Equation, Signals>> {
  using FluxMethod<Hll<Equation, Signals>>::FluxMethod;
};

template <typename Equation, typename Signals>
HllMethod(const Equation& eq, const Signals& signals)
    ->HllMethod<Equation, Signals>;

///////////////////////////////////////////////////////////////////////////////
//                            Implementation of HllBase<Equation, SignalSpeeds>

template <typename Equation, typename SignalSpeeds>
void HllBase<Equation, SignalSpeeds>::SolveRiemannProblem(Complete& solution,
                                                          const Complete& left,
                                                          const Complete& right,
                                                          Direction dir) {
  const auto signals = signal_speeds_(equation_, left, right, dir);
  Flux(GetEquation(), flux_left_, left, dir);
  Flux(GetEquation(), flux_right_, right, dir);
  const double sL = signals[0];
  const double sR = signals[1];
  const double ds = sR - sL;
  if (0.0 < sL) {
    solution = left;
  } else if (sR < 0.0) {
    solution = right;
  } else {
    ForEachComponent(
        [&](double& sol, double fL, double fR, double qL, double qR) {
          sol = (sR * qR - sL * qL + fL - fR) / ds;
        },
        AsCons(solution), flux_left_, flux_right_, AsCons(left), AsCons(right));
    CompleteFromCons(GetEquation(), solution, solution);
  }
}

template <typename Equation, typename SignalSpeeds>
void HllBase<Equation, SignalSpeeds>::ComputeNumericFlux(
    Conservative& numeric_flux, span<const Complete, 2> states,
    Duration /* dt */, double /* dx */, Direction dir) {
  const Complete& left = states[0];
  const Complete& right = states[1];

  const auto signals = signal_speeds_(GetEquation(), left, right, dir);

  Flux(GetEquation(), flux_left_, left, dir);
  Flux(GetEquation(), flux_right_, right, dir);

  const double sL = std::min(0.0, signals[0]);
  const double sR = std::max(0.0, signals[1]);
  const double sLsR = sL * sR;
  const double ds = sR - sL;
  FUB_ASSERT(ds > 0);

  ForEachComponent(
      [&](double& nf, double fL, double fR, double qL, double qR) {
        nf = (sR * fL - sL * fR + sLsR * (qR - qL)) / ds;
      },
      numeric_flux, flux_left_, flux_right_, states[0], states[1]);
}

template <typename Equation, typename SignalSpeeds>
double
HllBase<Equation, SignalSpeeds>::ComputeStableDt(span<const Complete, 2> states,
                                                 double dx, Direction dir) {
  const auto signals = signal_speeds_(GetEquation(), states[0], states[1], dir);
  const double max = std::accumulate(
      signals.begin(), signals.end(), 0.0,
      [](double x, double y) { return std::max(x, std::abs(y)); });
  FUB_ASSERT(max > 0.0);
  return dx / max;
}

///////////////////////////////////////////////////////////////////////////////
//                       Implementation of HllArrayBase<Equation, SignalSpeeds>

template <typename Equation, typename SignalSpeeds>
void HllArrayBase<Equation, SignalSpeeds, true>::SolveRiemannProblem(
    CompleteArray& solution, const CompleteArray& left,
    const CompleteArray& right, Direction dir) {
  const auto signals = GetSignalSpeeds()(GetEquation(), left, right, dir);
  Flux(GetEquation(), flux_left_array_, left, dir);
  Flux(GetEquation(), flux_right_array_, right, dir);
  const Array1d sL = signals[0];
  const Array1d sR = signals[1];
  const Array1d ds = sR - sL;
  ForEachComponent(
      [&](auto&& sol, Array1d fL, Array1d fR, Array1d qL, Array1d qR) {
        sol = (sR * qR - sL * qL + fL - fR) / ds;
        sol = (0.0 < sL).select(qL, (sR < 0.0).select(qR, sol));
      },
      AsCons(solution), flux_left_array_, flux_right_array_, AsCons(left),
      AsCons(right));
  //  CompleteFromCons(GetEquation(), solution, solution);
}

template <typename Equation, typename SignalSpeeds>
void HllArrayBase<Equation, SignalSpeeds, true>::ComputeNumericFlux(
    ConservativeArray& numeric_flux, span<const CompleteArray, 2> states,
    Duration /* dt */, double /* dx */, Direction dir) {
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  const auto signals = GetSignalSpeeds()(GetEquation(), left, right, dir);

  Flux(GetEquation(), flux_left_array_, left, dir);
  Flux(GetEquation(), flux_right_array_, right, dir);

  const Array1d zero = Array1d::Zero();
  const Array1d sL = signals[0].min(zero);
  const Array1d sR = signals[1].max(zero);
  const Array1d sLsR = sL * sR;
  const Array1d ds = sR - sL;

  ForEachComponent(
      [&](auto&& nf, Array1d fL, Array1d fR, Array1d qL, Array1d qR) {
        nf = (sR * fL - sL * fR + sLsR * (qR - qL)) / ds;
      },
      numeric_flux, flux_left_array_, flux_right_array_, AsCons(left),
      AsCons(right));
}

template <typename Equation, typename SignalSpeeds>
void HllArrayBase<Equation, SignalSpeeds, true>::ComputeNumericFlux(
    ConservativeArray& numeric_flux, Array1d face_fraction,
    span<const CompleteArray, 2> states, span<const Array1d, 2>, Duration /* dt */,
    double /* dx */, Direction dir) {
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  const auto signals = GetSignalSpeeds()(GetEquation(), left, right, dir);

  Flux(GetEquation(), flux_left_array_, left, dir);
  Flux(GetEquation(), flux_right_array_, right, dir);

  const Array1d zero = Array1d::Zero();
  const Array1d sL = signals[0].min(zero);
  const Array1d sR = signals[1].max(zero);
  const Array1d sLsR = sL * sR;
  const Array1d ds = sR - sL;

  MaskArray mask = face_fraction > 0.0;
  ForEachComponent(
      [&](auto&& nf, Array1d fL, Array1d fR, Array1d qL, Array1d qR) {
        nf = mask.select((sR * fL - sL * fR + sLsR * (qR - qL)) / ds, Array1d::Zero());
      },
      numeric_flux, flux_left_array_, flux_right_array_, AsCons(left),
      AsCons(right));
}

template <typename Equation, typename SignalSpeeds>
Array1d HllArrayBase<Equation, SignalSpeeds, true>::ComputeStableDt(
    span<const CompleteArray, 2> states, double dx, Direction dir) {
  const auto signals =
      GetSignalSpeeds()(GetEquation(), states[0], states[1], dir);
  Array1d zero = Array1d::Zero();
  const Array1d max =
      std::accumulate(signals.begin(), signals.end(), zero,
                      [](Array1d x, Array1d y) { return x.max(y.abs()); });
  return Array1d(dx) / max;
}

} // namespace fub

#endif
