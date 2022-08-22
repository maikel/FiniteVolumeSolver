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

#ifndef FUB_EQUATIONS_PERFECT_GAS_HLLEM_HPP
#define FUB_EQUATIONS_PERFECT_GAS_HLLEM_HPP

#include "fub/equations/PerfectGas.hpp"
#include "fub/flux_method/FluxMethod.hpp"

namespace fub::perfect_gas {

/// \ingroup FluxMethod
template <typename EulerEquation, bool Larrouturou = true> struct Hllem {
  using Conservative = ::fub::Conservative<EulerEquation>;
  using Complete = ::fub::Complete<EulerEquation>;
  using ConservativeArray = ::fub::ConservativeArray<EulerEquation>;
  using CompleteArray = ::fub::CompleteArray<EulerEquation>;

  Hllem(const EulerEquation& equation) : equation_{equation} {}

  const EulerEquation& GetEquation() const noexcept { return equation_; }

  static constexpr int GetStencilWidth() noexcept { return 1; }

  double ComputeNumericFlux(Conservative& flux, span<const Complete, 2> stencil,
                          Duration dt, double dx, Direction dir);

  void ComputeNumericFlux(ConservativeArray& flux, Array1d face_fractions,
                          span<const CompleteArray, 2> stencil,
                          span<const Array1d, 2> volume_fractions, Duration dt,
                          double dx, Direction dir);

  Array1d ComputeNumericFlux(ConservativeArray& flux,
                          span<const CompleteArray, 2> stencil, Duration dt,
                          double dx, Direction dir);

  double ComputeStableDt(span<const Complete, 2> states, double dx,
                         Direction dir);

  Array1d ComputeStableDt(span<const CompleteArray, 2> states, double dx,
                          Direction dir);

  Array1d ComputeStableDt(span<const CompleteArray, 2> states,
                          Array1d face_fraction, span<const Array1d, 2>,
                          double dx, Direction dir);

  void SolveRiemannProblem(Complete& solution, const Complete& left,
                           const Complete& right, Direction dir);

  EulerEquation equation_;
  Complete left_{equation_};
  Complete right_{equation_};
  Conservative fluxL_{equation_};
  Conservative fluxR_{equation_};
  Conservative flux_hlle_{equation_};
  Conservative w_hlle_{equation_};
  Conservative w_hllem_{equation_};

  CompleteArray left_array_{equation_};
  CompleteArray right_array_{equation_};
  ConservativeArray fluxL_array_{equation_};
  ConservativeArray fluxR_array_{equation_};
  ConservativeArray flux_hlle_array_{equation_};
};

extern template struct Hllem<PerfectGas<1>, false>;
extern template struct Hllem<PerfectGas<2>, false>;
extern template struct Hllem<PerfectGas<3>, false>;

extern template struct Hllem<PerfectGas<1>, true>;
extern template struct Hllem<PerfectGas<2>, true>;
extern template struct Hllem<PerfectGas<3>, true>;

/// \ingroup FluxMethod
template <typename EulerEquation, bool Larrouturou = true>
using HllemMethod = FluxMethod<Hllem<EulerEquation, Larrouturou>>;

} // namespace fub::perfect_gas

namespace fub {

extern template class FluxMethod<perfect_gas::Hllem<PerfectGas<1>, false>>;
extern template class FluxMethod<perfect_gas::Hllem<PerfectGas<2>, false>>;
extern template class FluxMethod<perfect_gas::Hllem<PerfectGas<3>, false>>;

extern template class FluxMethod<perfect_gas::Hllem<PerfectGas<1>, true>>;
extern template class FluxMethod<perfect_gas::Hllem<PerfectGas<2>, true>>;
extern template class FluxMethod<perfect_gas::Hllem<PerfectGas<3>, true>>;

} // namespace fub

#endif