// Copyright (c) 2020 Maikel Nadolski
// Copyright (c) 2020 Abbas Gholami Poshtehani
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

#ifndef FUB_EQUATIONS_PERFECT_GAS_THIRD_ORDER_RUNGE_KUTTA_HPP
#define FUB_EQUATIONS_PERFECT_GAS_THIRD_ORDER_RUNGE_KUTTA_HPP

#include "fub/equations/PerfectGas.hpp"
#include "fub/flux_method/FluxMethod.hpp"

namespace fub::perfect_gas {


/// \ingroup FluxMethod
template <int Dim> struct ThirdOrderRungeKutta {
  using Conservative = ::fub::Conservative<PerfectGas<Dim>>;
  using Complete = ::fub::Complete<PerfectGas<Dim>>;
  using ConservativeArray = ::fub::ConservativeArray<PerfectGas<Dim>>;
  using CompleteArray = ::fub::CompleteArray<PerfectGas<Dim>>;

  ThirdOrderRungeKutta(const PerfectGas<Dim>& equation) : equation_{equation} {}

  const PerfectGas<Dim>& GetEquation() const noexcept { return equation_; }

  static constexpr int GetStencilWidth() noexcept { return 6; }

  void ComputeNumericFlux(Conservative& flux, span<const Complete, 12> stencil,
                          Duration dt, double dx, Direction dir);

  void ComputeNumericFlux(ConservativeArray& flux, Array1d face_fractions,
                          span<const CompleteArray, 12> stencil,
                          span<const Array1d, 12> volume_fractions, Duration dt,
                          double dx, Direction dir);

  void ComputeNumericFlux(ConservativeArray& flux,
                          span<const CompleteArray, 12> stencil, Duration dt,
                          double dx, Direction dir);

  double ComputeStableDt(span<const Complete, 12> states, double dx,
                         Direction dir);

  Array1d ComputeStableDt(span<const CompleteArray, 12> states, double dx,
                          Direction dir);

  Array1d ComputeStableDt(span<const CompleteArray, 12> states,
                          Array1d face_fraction, span<const Array1d, 2>,
                          double dx, Direction dir);

  PerfectGas<Dim> equation_;
  double k_0{0.0};
  double eta_0{0.0};
};

extern template struct ThirdOrderRungeKutta<1>;

/// \ingroup FluxMethod
template <int Dim> using ThirdOrderRungeKuttaMethod = FluxMethod<ThirdOrderRungeKutta<Dim>>;

} // namespace fub::perfect_gas

namespace fub {

extern template class FluxMethod<perfect_gas::ThirdOrderRungeKutta<1>>;

} // namespace fub

#endif