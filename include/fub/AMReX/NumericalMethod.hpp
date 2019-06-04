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

#ifndef FUB_AMREX_NUMERICAL_METHOD_HPP
#define FUB_AMREX_NUMERICAL_METHOD_HPP

#include "fub/AMReX/FluxMethod.hpp"
#include "fub/AMReX/HyperbolicSplitTimeIntegrator.hpp"
#include "fub/AMReX/Reconstruction.hpp"

namespace fub::amrex {
template <typename FM> struct FluxMethodWrapper;
/// \brief This struct summarizes the numerical method which is being deployed
/// by this integrator context.
struct NumericalMethod {
  NumericalMethod(FluxMethod, HyperbolicSplitTimeIntegrator, Reconstruction);

  template <typename FM, typename = std::enable_if_t<!std::is_same_v<
                             std::decay_t<FM>, NumericalMethod>>>
  explicit NumericalMethod(const FM& method)
      : flux_method(FluxMethodWrapper<std::decay_t<FM>>(method)),
        time_integrator{ForwardIntegrator{}},
        reconstruction{ReconstructEquationStates{method.GetEquation()}} {}

  FluxMethod flux_method;
  HyperbolicSplitTimeIntegrator time_integrator;
  Reconstruction reconstruction;
};

} // namespace fub::amrex

#endif
