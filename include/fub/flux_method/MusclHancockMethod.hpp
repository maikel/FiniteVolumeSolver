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

#ifndef FUB_IDEAL_GAS_FLUX_METHOD_MUSCL_HANCOCK_METHOD
#define FUB_IDEAL_GAS_FLUX_METHOD_MUSCL_HANCOCK_METHOD

#include "fub/CompleteFromCons.hpp"
#include "fub/Equation.hpp"
#include "fub/core/span.hpp"
#include "fub/flux_method/FluxMethod.hpp"
#include "fub/flux_method/GodunovMethod.hpp"

namespace fub {
struct MinMod {
  template <typename Equation>
  void ComputeLimitedSlope(Conservative<Equation>& cons,
                           span<const Complete<Equation>, 3> stencil) {
    ForEachComponent(
        [](double& cons, double qL, double qM, double qR) {
          const double sL = qM - qL;
          const double sR = qR - qM;
          if (sR > 0) {
            cons = std::max(0.0, std::min(sL, sR));
          } else {
            cons = std::min(0.0, std::max(sL, sR));
          }
        },
        cons, AsCons(stencil[0]), AsCons(stencil[1]), AsCons(stencil[2]));
  }

  template <typename Equation, int N>
  void ComputeLimitedSlope(ConservativeArray<Equation, N>& cons,
                           span<const CompleteArray<Equation, N>, 3> stencil) {
    ForEachComponent(
        [](auto&& cons, auto qL, auto qM, auto qR) {
          const Array<double, 1, N> sL = qM - qL;
          const Array<double, 1, N> sR = qR - qM;
          const Array<double, 1, N> zero = Array<double, 1, N>::Constant(0.0);
          cons = (sR > 0).select(zero.max(sL.min(sR)), zero.min(sL.max(sR)));
        },
        cons, AsCons(stencil[0]), AsCons(stencil[1]), AsCons(stencil[2]));
  }
};

struct VanLeer {
  template <typename Equation>
  void ComputeLimitedSlope(Conservative<Equation>& cons,
                           span<const Complete<Equation>, 3> stencil) {
    ForEachComponent(
        [](double& cons, double qL, double qM, double qR) {
          const double sL = qM - qL;
          const double sR = qR - qM;
          double r = 0.0;
          if (sL * sR > 0.0) {
            r = sL / sR;
          }
          if (r < 0.0) {
            cons = 0.0;
          } else {
            cons = 0.5 * std::min(2 * r / (1 + r), 2 / (1 + r)) * (sL + sR);
          }
        },
        cons, AsCons(stencil[0]), AsCons(stencil[1]), AsCons(stencil[2]));
  }
};

template <typename EquationT,
          typename BaseMethod =
              GodunovMethod<EquationT, ExactRiemannSolver<EquationT>>,
          typename SlopeLimiter = MinMod>
struct MusclHancock {
  using Equation = EquationT;
  using Complete = typename Equation::Complete;
  using Conservative = typename Equation::Conservative;

  using CompleteArray = ::fub::CompleteArray<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;

  explicit MusclHancock(const Equation& eq) : equation_{eq}, flux_method_{eq} {}

  MusclHancock(const Equation& eq, const BaseMethod& method)
      : equation_{eq}, flux_method_{method} {}

  static constexpr int GetStencilWidth() noexcept { return 2; }

  double ComputeStableDt(span<const Complete, 4> states, double dx,
                         Direction dir) noexcept {
    return flux_method_.ComputeStableDt(states.template subspan<1, 2>(), dx,
                                        dir);
  }

  Array1d ComputeStableDt(span<const CompleteArray, 4> states, double dx,
                          Direction dir) noexcept {
    return flux_method_.ComputeStableDt(states.template subspan<1, 2>(), dx,
                                        dir);
  }

  void ComputeNumericFlux(Conservative& flux, span<const Complete, 4> stencil,
                          Duration dt, double dx, Direction dir) {
    const double lambda_half = 0.5 * dt.count() / dx;

    ////////////////////////////////////////////////////////////////////////////
    // Compute Left Reconstructed Complete State

    slope_limiter_.ComputeLimitedSlope(slope_, stencil.template first<3>());

    ForEachComponent(
        [](double& qL, double& qR, double state, double slope) {
          qL = state - 0.5 * slope;
          qR = state + 0.5 * slope;
        },
        AsCons(q_left_), AsCons(q_right_), AsCons(stencil[1]), slope_);

    CompleteFromCons(equation_, q_left_, q_left_);
    CompleteFromCons(equation_, q_right_, q_right_);

    Flux(equation_, flux_left_, q_left_, dir);
    Flux(equation_, flux_right_, q_right_, dir);

    ForEachComponent(
        [&lambda_half](double& rec, double qR, double fL, double fR) {
          rec = qR + lambda_half * (fL - fR);
        },
        AsCons(rec_[0]), AsCons(q_right_), flux_left_, flux_right_);

    CompleteFromCons(equation_, rec_[0], rec_[0]);

    ///////////////////////////////////////////////////////////////////////////
    // Compute Right Reconstructed Complete State

    slope_limiter_.ComputeLimitedSlope(slope_, stencil.template last<3>());

    ForEachComponent(
        [](double& qL, double& qR, double state, double slope) {
          qL = state - 0.5 * slope;
          qR = state + 0.5 * slope;
        },
        AsCons(q_left_), AsCons(q_right_), AsCons(stencil[2]), slope_);

    CompleteFromCons(equation_, q_left_, q_left_);
    CompleteFromCons(equation_, q_right_, q_right_);

    Flux(equation_, flux_left_, q_left_, dir);
    Flux(equation_, flux_right_, q_right_, dir);

    ForEachComponent(
        [&lambda_half](double& rec, double qL, double fL, double fR) {
          rec = qL + lambda_half * (fL - fR);
        },
        AsCons(rec_[1]), AsCons(q_left_), flux_left_, flux_right_);

    CompleteFromCons(equation_, rec_[1], rec_[1]);

    ///////////////////////////////////////////////////////////////////////////
    // Invoke Lower Order Flux Method

    flux_method_.ComputeNumericFlux(flux, span{rec_}, dt, dx, dir);
  }

  void ComputeNumericFlux(ConservativeArray& flux,
                          span<const CompleteArray, 4> stencil, Duration dt,
                          double dx, Direction dir) {
    const Array1d lambda_half = Array1d::Constant(0.5 * dt.count() / dx);

    ////////////////////////////////////////////////////////////////////////////
    // Compute Left Reconstructed Complete State

    slope_limiter_.ComputeLimitedSlope(slope_arr_, stencil.template first<3>());

    ForEachComponent(
        [](auto&& qL, auto&& qR, const auto& state, const auto& slope) {
          qL = state - 0.5 * slope;
          qR = state + 0.5 * slope;
        },
        AsCons(q_left_arr_), AsCons(q_right_arr_), AsCons(stencil[1]),
        slope_arr_);

    CompleteFromCons(equation_, q_left_arr_, q_left_arr_);
    CompleteFromCons(equation_, q_right_arr_, q_right_arr_);

    Flux(equation_, flux_left_arr_, q_left_arr_, dir);
    Flux(equation_, flux_right_arr_, q_right_arr_, dir);

    ForEachComponent(
        [&lambda_half](auto&& rec, auto qR, auto fL, auto fR) {
          rec = qR + lambda_half * (fL - fR);
        },
        AsCons(rec_arr_[0]), AsCons(q_right_arr_), flux_left_arr_,
        flux_right_arr_);

    CompleteFromCons(equation_, rec_arr_[0], rec_arr_[0]);

    ///////////////////////////////////////////////////////////////////////////
    // Compute Right Reconstructed Complete State

    slope_limiter_.ComputeLimitedSlope(slope_arr_, stencil.template last<3>());

    ForEachComponent(
        [](auto&& qL, auto&& qR, const auto& state, const auto& slope) {
          qL = state - 0.5 * slope;
          qR = state + 0.5 * slope;
        },
        AsCons(q_left_arr_), AsCons(q_right_arr_), AsCons(stencil[2]),
        slope_arr_);

    CompleteFromCons(equation_, q_left_arr_, q_left_arr_);
    CompleteFromCons(equation_, q_right_arr_, q_right_arr_);

    Flux(equation_, flux_left_arr_, q_left_arr_, dir);
    Flux(equation_, flux_right_arr_, q_right_arr_, dir);

    ForEachComponent(
        [&lambda_half](auto&& rec, auto qL, auto fL, auto fR) {
          rec = qL + lambda_half * (fL - fR);
        },
        AsCons(rec_arr_[1]), AsCons(q_left_arr_), flux_left_arr_,
        flux_right_arr_);

    CompleteFromCons(equation_, rec_arr_[1], rec_arr_[1]);

    ///////////////////////////////////////////////////////////////////////////
    // Invoke Lower Order Flux Method

    flux_method_.ComputeNumericFlux(flux, span{rec_arr_}, dt, dx, dir);
  }

  const Equation& GetEquation() const noexcept { return equation_; }
  Equation& GetEquation() noexcept { return equation_; }

private:
  // These member variables control the behaviour of this method
  Equation equation_;
  BaseMethod flux_method_;
  SlopeLimiter slope_limiter_;

  // These member variables are used as an intermediate computation storage
  // and need to be allocated at construction time.
  std::array<Complete, 2> rec_{Complete{equation_}, Complete{equation_}};
  Conservative slope_{equation_};     //< Storage for the limited slopes
  Conservative flux_left_{equation_}; //< Storage for an intermediate left flux
  Conservative flux_right_{
      equation_};               //< Storage for an intermediate right flux
  Complete q_left_{equation_};  //< Storage for an intermediate left state
  Complete q_right_{equation_}; //< Storage for an intermediate right state

  std::array<CompleteArray, 2> rec_arr_{CompleteArray{equation_},
                                        CompleteArray{equation_}};
  ConservativeArray slope_arr_{equation_}; //< Storage for the limited slopes
  ConservativeArray flux_left_arr_{
      equation_}; //< Storage for an intermediate left flux
  ConservativeArray flux_right_arr_{
      equation_}; //< Storage for an intermediate right flux
  CompleteArray q_left_arr_{
      equation_}; //< Storage for an intermediate left state
  CompleteArray q_right_arr_{
      equation_}; //< Storage for an intermediate right state
};

template <typename Equation,
          typename BaseMethod =
              fub::GodunovMethod<Equation, ExactRiemannSolver<Equation>>,
          typename Slope = MinMod>
struct MusclHancockMethod
    : FluxMethod<MusclHancock<Equation, BaseMethod, Slope>> {
  using FluxMethod<MusclHancock<Equation, BaseMethod, Slope>>::FluxMethod;
};

template <typename Equation>
MusclHancockMethod(const Equation&)->MusclHancockMethod<Equation>;

template <typename Equation, typename Method>
MusclHancockMethod(const Equation&, const Method&)
    ->MusclHancockMethod<Equation, Method>;

} // namespace fub

#endif
