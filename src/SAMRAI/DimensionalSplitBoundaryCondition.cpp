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

#include "fub/SAMRAI/DimensionalSplitBoundaryCondition.hpp"

namespace fub {

void DimensionalSplitBoundaryConditions::SetBoundaryCondition(
    std::shared_ptr<DimensionalSplitBoundaryCondition> condition, Direction dir,
    int side) {
  int i = 2 * static_cast<int>(dir) + side;
  conds_[i] = std::move(condition);
}

void DimensionalSplitBoundaryConditions::SetPhysicalBoundaryConditions(
    const SAMRAI::hier::Patch& patch, double fill_time,
    const SAMRAI::hier::IntVector& ghost_width_to_fill) const {
  const SAMRAI::hier::PatchGeometry& geometry = *patch.getPatchGeometry();
  const std::vector<SAMRAI::hier::BoundaryBox>& boundaries =
      geometry.getCodimensionBoundaries(1);
  for (const SAMRAI::hier::BoundaryBox& bbox : boundaries) {
    SAMRAI::hier::Box fill_box =
        geometry.getBoundaryFillBox(bbox, patch.getBox(), ghost_width_to_fill);
    if (fill_box.size() > 0) {
      const int direction = bbox.getLocationIndex() / 2;
      const int side = bbox.getLocationIndex() % 2;
      const int index = 2 * direction + side;
      if (conds_[index]) {
        conds_[index]->SetPhysicalBoundaryCondition(patch, fill_box, fill_time,
                                                    Direction(direction), side);
      }
    }
  }
}

SAMRAI::hier::IntVector DimensionalSplitBoundaryConditions::GetStencilWidth(
    const SAMRAI::tbox::Dimension& dim) const {
  SAMRAI::hier::IntVector gcw(dim, 0);
  for (int i = 0; i < dim.getValue(); ++i) {
    if (conds_[2 * i]) {
      gcw[i] = std::max(gcw[i], conds_[2 * i]->GetStencilWidth());
    }
    if (conds_[2 * i + 1]) {
      gcw[i] = std::max(gcw[i], conds_[2 * i + 1]->GetStencilWidth());
    }
  }
  return gcw;
}

} // namespace fub