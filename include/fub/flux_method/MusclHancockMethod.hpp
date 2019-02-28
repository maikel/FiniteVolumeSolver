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

#include "fub/Equation.hpp"
#include "fub/core/span.hpp"

namespace fub {

template <typename Equation, typename FluxMethod, typename Limiter> 
struct MusclHancockMethod {
  using Complete = typename Equation::Complete;
  using Cons = typename Equation::Cons;

  MusclHancockMethod(const Equation& eq, const FluxMethod& method)
      : equation_{eq}, flux_method_{method} {}

  static void ComputeSlope(Cons& slope, span<const Complete, 3> stencil,
                           double cell_width) {
    ForEachComponent(
        [](double& slope, double qL, double qR) { slope = 0.5 * (qR - qL); }, 
        slope, stencil[0], stencil[2]);
  }

  void ComputeLimitedSlope(Cons& slope, span<const Complete, 3> stencil,
                           double cell_width) {
    ComputeSlope(slope, stencil, cell_width);
    limiter_.LimitSlope(slope, stencil, cell_width);
  }

  void ComputeNumericFlux(Cons& flux, span<const Complete, 4> stencil,
                          Duration dt, double dx) noexcept {
    const double dt_2 = 0.5 * dt;

    ////////////////////////////////////////////////////////////////////////////
    // Compute Left Reconstructed Complete State

    ComputeLimitedSlope(slope_, stencil.template first<3>(), dx);

    ForEachComponent<Conservative>(
        [](double& qL, double& qR, double state, double slope) {
          qL = state - slope;
          qR = state + slope;
        },
        q_left_, q_right_, stencil[1], slope_);
    Reconstruct(equation, q_left_, q_left_);
    Reconstruct(equation, q_right_, q_right_);

    equation_.Flux(flux_left_, q_left_);
    equation_.Flux(flux_right_, q_right_);

    ForEachComponent<Conservative>(
        [&dt_2](double& rec, double qR, double fL, double fR) {
          rec = qR + dt_2 * (fL - fR);
        },
        rec_[0], q_right_, stencil[1], slope_);
    Reconstruct(equation_, rec_[0], rec_[0]);

    ///////////////////////////////////////////////////////////////////////////
    // Compute Right Reconstructed Complete State

    ComputeLimitedSlope(slope_, stencil.template last<3>(), dx);

    ForEachComponent<Conservative>(
        [](auto qL, auto qR, auto state, auto slope) {
          qL = state - slope;
          qR = state + slope;
        },
        AsCons(q_left_), AsCons(q_right_), AsCons(stencil[2]), slope_);

    Reconstruct(equation_, q_left_, AsCons(q_left_));
    Reconstruct(equation_, q_right_, AsCons(q_right_));

    Flux(flux_left_, q_left_);
    Flux(flux_right_, q_right_);

    ForEachComponent<Conservative>(
      [&dt_2](double& rec, auto qL, auto fL, auto fR) { 
        rec = qL + dt_2 * (fL - fR); 
      }, rec_[1], q_left_, stencil[2], slope_);

    Reconstruct(equation_, rec_[1], AsCons(rec_[1]));

    ///////////////////////////////////////////////////////////////////////////
    // Invoke Lower Order Flux Method

    flux_method_.ComputeNumericFlux(flux, span{rec_}, dt, dx);
  }

private:
  // These member variables control the behaviour of this method
  Equation equation_;
  FluxMethod flux_method_;
  Limiter limiter_;

  // These member variables are used as an intermediate computation storage
  // and need to be allocated at construction time.
  std::array<Complete, 4> stencil_;
  std::array<Complete, 2> rec_;
  Cons numeric_flux_;
  Cons slope_;      //< Storage for the limited slopes
  Cons flux_left_;  //< Storage for an intermediate left flux
  Cons flux_right_; //< Storage for an intermediate right flux
  Complete q_left_;    //< Storage for an intermediate left state
  Complete q_right_;   //< Storage for an intermediate right state
};

} // namespace fub

#endif