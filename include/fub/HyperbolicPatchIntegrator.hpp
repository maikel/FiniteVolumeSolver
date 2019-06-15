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

#include "fub/Duration.hpp"
#include "fub/Execution.hpp"
#include "fub/PatchDataView.hpp"

namespace fub {

template <typename Tag> struct HyperbolicPatchIntegrator {
  template <typename T, int Rank>
  using PatchDataView = ::fub::PatchDataView<T, Rank, layout_stride>;

  HyperbolicPatchIntegrator(Tag);

  static void UpdateConservatively(const PatchDataView<double, 2>& next,
                                   const PatchDataView<const double, 2>& prev,
                                   const PatchDataView<const double, 2>& fluxes,
                                   Duration dt, double dx, Direction dir);

  static void UpdateConservatively(const PatchDataView<double, 3>& next,
                                   const PatchDataView<const double, 3>& prev,
                                   const PatchDataView<const double, 3>& fluxes,
                                   Duration dt, double dx, Direction dir);

  static void UpdateConservatively(const PatchDataView<double, 4>& next,
                                   const PatchDataView<const double, 4>& prev,
                                   const PatchDataView<const double, 4>& fluxes,
                                   Duration dt, double dx, Direction dir);
};

extern template struct HyperbolicPatchIntegrator<execution::SequentialTag>;
extern template struct HyperbolicPatchIntegrator<execution::SimdTag>;
extern template struct HyperbolicPatchIntegrator<execution::OpenMpTag>;
extern template struct HyperbolicPatchIntegrator<execution::OpenMpSimdTag>;

} // namespace fub

#endif