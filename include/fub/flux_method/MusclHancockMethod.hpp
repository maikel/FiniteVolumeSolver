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

template <typename Equation, typename FluxMethod> struct MusclHancockMethod {
  template <typename T> using StateView = StateView<T, Equation>;
  template <typename T> using FluxView = ConsView<T, Equation>;

  using State = SimdifyT<typename Equation::State>;
  using Cons = SimdifyT<typename Equation::Cons>;

  MusclHancockMethod(const Equation& eq, const FluxMethod& method)
      : equation_{eq}, flux_method_{method} {}

  static void ComputeSlope(Cons& slope, span<const State, 3> stencil,
                           double cell_width) {
    ForEachVariable(
        [](auto slope, auto qL, auto qR) { slope = 0.5 * (qR - qL); }, slope,
        stencil[0], stencil[2]);
  }

  void ComputeNumericFlux(Cons& flux, span<const State, 4> stencil, double dt,
                          double dx) noexcept {
    const Array1d dt_2 = Array1d::Constant(0.5 * dt);

    ////////////////////////////////////////////////////////////////////////////
    // Compute Left Reconstructed State

    ComputeSlope(slope_, stencil.template first<3>(), dx);

    ForEachVariable(
        [](auto qL, auto qR, auto state, auto slope) {
          qL = state - slope;
          qR = state + slope;
        },
        AsCons(q_left_), AsCons(q_right_), AsCons(stencil[1]), slope_);

    equation_.Reconstruct(q_left_, AsCons(q_left_));
    equation_.Reconstruct(q_right_, AsCons(q_right_));

    equation_.Flux(flux_left_, q_left_);
    equation_.Flux(flux_right_, q_right_);

    ForEachVariable([&dt_2](auto rec, auto qR, auto fL,
                            auto fR) { rec = qR + dt_2 * (fL - fR); },
                    AsCons(rec_[0]), AsCons(q_right_), AsCons(stencil[1]),
                    slope_);

    equation_.Reconstruct(rec_[0], AsCons(rec_[0]));

    ///////////////////////////////////////////////////////////////////////////
    // Compute Right Reconstructed State

    ComputeSlope(slope_, stencil.template last<3>(), dx);

    ForEachVariable(
        [](auto qL, auto qR, auto state, auto slope) {
          qL = state - slope;
          qR = state + slope;
        },
        AsCons(q_left_), AsCons(q_right_), AsCons(stencil[2]), slope_);

    equation_.Reconstruct(q_left_, AsCons(q_left_));
    equation_.Reconstruct(q_right_, AsCons(q_right_));

    equation_.Flux(flux_left_, q_left_);
    equation_.Flux(flux_right_, q_right_);

    ForEachVariable([&dt_2](auto rec, auto qL, auto fL,
                            auto fR) { rec = qL + dt_2 * (fL - fR); },
                    AsCons(rec_[1]), AsCons(q_left_), AsCons(stencil[2]),
                    slope_);

    equation_.Reconstruct(rec_[1], AsCons(rec_[1]));

    ///////////////////////////////////////////////////////////////////////////
    // Invoke Lower Order Flux Method

    flux_method_.ComputeNumericFlux(flux, span{rec_}, dt, dx);
  }

  void ComputeNumericFluxesOnSpans(FluxView<double> fluxes,
                                   StateView<const double> states, double dt,
                                   double dx) noexcept {
    ForEachRow(
        Extents(fluxes),
        [&](auto flux, auto state) {
          for (int i = 0; i < Extents(flux).extent(0); i += kChunkSize) {
            Load(span(stencil_), state, i);
            ComputeNumericFlux(numeric_flux_, stencil_, dt, dx);
            Store(flux, numeric_flux_, i);
          }
        },
        fluxes, states);
  }

private:
  // These member variables control the behaviour of this method
  Equation equation_;
  FluxMethod flux_method_;

  // These member variables are used as an intermediate computation storage
  // and need to be allocated at construction time.
  std::array<State, 4> stencil_{State(equation_), State(equation_),
                                State(equation_), State(equation_)};
  std::array<State, 2> rec_{State(equation_), State(equation_)};
  Cons numeric_flux_{equation_};
  Cons slope_{equation_};      //< Storage for the limited slopes
  Cons flux_left_{equation_};  //< Storage for an intermediate left flux
  Cons flux_right_{equation_}; //< Storage for an intermediate right flux
  State q_left_{equation_};    //< Storage for an intermediate left state
  State q_right_{equation_};   //< Storage for an intermediate right state
};

} // namespace fub

#endif