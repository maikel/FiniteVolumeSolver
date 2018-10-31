// Copyright (c) 2018 Maikel Nadolski
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

#ifndef FUB_SOLVER_EULER_GODUNOV_METHOD_HPP
#define FUB_SOLVER_EULER_GODUNOV_METHOD_HPP

#include "fub/solver/DimensionalSplitFluxMethod.hpp"
#include "fub/solver/euler/IdealGas.hpp"

namespace fub {
namespace euler {

struct GodunovMethod
    : fub::DimensionalSplitFluxMethod<IdealGas::FluxStates,
                                      IdealGas::CompleteStates> {
  using FluxStates = IdealGas::FluxStates;
  using CompleteStates = IdealGas::CompleteStates;

  double
  estimateStableDtOnPatch(const StateData& states,
                          const SAMRAI::hier::Patch& patch) const override;

  void computeFluxesOnPatch(const FluxData& fluxes, const StateData& states,
                            const SAMRAI::hier::Patch& patch, double dt,
                            int dir) const override;

  SAMRAI::hier::IntVector
  getStencilWidth(const SAMRAI::tbox::Dimension& dim) const override;
};

}
} // namespace fub::euler

#endif