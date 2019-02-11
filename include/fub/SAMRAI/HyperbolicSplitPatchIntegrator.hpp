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

#ifndef FUB_SAMRAI_HYPERBOLIC_SPLIT_PATCH_INTEGRATOR_HPP
#define FUB_SAMRAI_HYPERBOLIC_SPLIT_PATCH_INTEGRATOR_HPP

#include "fub/core/span.hpp"

#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/OutersideData.h>
#include <SAMRAI/pdat/SideData.h>

namespace fub {
namespace samrai {

/// \ingroup Abstract
/// This is an abstract
struct HyperbolicSplitPatchIntegrator {
  virtual ~HyperbolicSplitPatchIntegrator() = default;

  virtual void
  ConservativeUpdateOnPatch(span<SAMRAI::pdat::CellData<double>*> next,
                            span<const SAMRAI::pdat::SideData<double>*> fluxes,
                            span<const SAMRAI::pdat::CellData<double>*> prev,
                            const SAMRAI::hier::Patch& patch, double dt,
                            int dir) = 0;

  virtual void
  ReconstructStatesFromCons(span<SAMRAI::pdat::CellData<double>*> states,
                            span<const SAMRAI::pdat::CellData<double>*> cons,
                            const SAMRAI::hier::Patch& patch, int dir) = 0;
};

} // namespace samrai
} // namespace fub

#endif