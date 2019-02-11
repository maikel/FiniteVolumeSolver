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

#ifndef FUB_SAMRAI_SPLIT_BOUNDARY_CONDITION_HPP
#define FUB_SAMRAI_SPLIT_BOUNDARY_CONDITION_HPP

#include "fub/SAMRAI/BoundaryCondition.hpp"
#include "fub/Direction.hpp"

#include <array>

namespace fub {

struct DimensionalSplitBoundaryCondition {
  virtual ~DimensionalSplitBoundaryCondition() = default;

  virtual void SetPhysicalBoundaryCondition(const SAMRAI::hier::Patch& patch,
                                            const SAMRAI::hier::Box& fill_box,
                                            double fill_time, Direction dir,
                                            int side) const = 0;

  virtual int GetStencilWidth() const = 0;

  virtual void
  PreFillGhostLayer(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&,
                    double) {}
};

class DimensionalSplitBoundaryConditions : public BoundaryCondition {
public:
  void SetPhysicalBoundaryConditions(
      const SAMRAI::hier::Patch& patch, double fill_time,
      const SAMRAI::hier::IntVector& ghost_width_to_fill) const override;

  SAMRAI::hier::IntVector
  GetStencilWidth(const SAMRAI::tbox::Dimension& dim) const override;

  void
  SetBoundaryCondition(std::shared_ptr<DimensionalSplitBoundaryCondition> condition,
                       Direction dir, int side);

  void
  PreFillGhostLayer(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hier,
                    double timepoint, Direction dir) override {
    const int d = static_cast<int>(dir);
    conds_[2 * d]->PreFillGhostLayer(hier, timepoint);
    conds_[2 * d + 1]->PreFillGhostLayer(hier, timepoint);
  }

private:
  static constexpr int max_n_conds = 2 * SAMRAI_MAXIMUM_DIMENSION;
  std::array<std::shared_ptr<DimensionalSplitBoundaryCondition>, max_n_conds> conds_{};
};

} // namespace fub

#endif
