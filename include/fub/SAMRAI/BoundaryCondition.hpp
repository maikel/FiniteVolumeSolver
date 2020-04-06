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

#ifndef FUB_SAMRAI_BOUNDARY_CONDITION_HPP
#define FUB_SAMRAI_BOUNDARY_CONDITION_HPP

#include "fub/Duration.hpp"

#include <SAMRAI/xfer/RefinePatchStrategy.h>

#include <functional>
#include <type_traits>

namespace fub::samrai {

class BoundaryCondition : public SAMRAI::xfer::RefinePatchStrategy {
public:
  BoundaryCondition() = default;

  template <typename BC, typename = std::enable_if_t<!std::is_same_v<
                             std::decay_t<BC>, BoundaryCondition>>>
  BoundaryCondition(BC&& some_condition)
      : fn_{[cond = std::forward<BC>(some_condition)](
                SAMRAI::hier::Patch& patch, Duration t,
                const SAMRAI::hier::IntVector& ghosts) {
          return cond.SetPhysicalBoundaryConditions(patch, t, ghosts);
        }} {}

  void setPhysicalBoundaryConditions(
      SAMRAI::hier::Patch& patch, double fill_time,
      const SAMRAI::hier::IntVector& ghost_width_to_fill) override {
    if (fn_) {
      fn_(patch, Duration(fill_time), ghost_width_to_fill);
    }
  }

  void preprocessRefine(SAMRAI::hier::Patch&, const SAMRAI::hier::Patch&,
                        const SAMRAI::hier::Box&,
                        const SAMRAI::hier::IntVector&) override {}

  void postprocessRefine(SAMRAI::hier::Patch&, const SAMRAI::hier::Patch&,
                         const SAMRAI::hier::Box&,
                         const SAMRAI::hier::IntVector&) override {}

  SAMRAI::hier::IntVector
  getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const override {
    return SAMRAI::hier::IntVector::getZero(dim);
  }

private:
  using PhysicalBoundaryFunction = std::function<void(
      SAMRAI::hier::Patch&, Duration, const SAMRAI::hier::IntVector&)>;
  PhysicalBoundaryFunction fn_{};
};

} // namespace fub::samrai

#endif
