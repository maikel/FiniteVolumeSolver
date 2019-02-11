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

#include "fub/SAMRAI/tagging/GradientDetector.hpp"

#include "fub/SAMRAI/utility.hpp"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

namespace fub {
namespace {
double computePartialGradient_(const SAMRAI::pdat::CellData<double>& q,
                               const SAMRAI::pdat::CellIndex& mid,
                               const SAMRAI::geom::CartesianPatchGeometry& geom,
                               Direction dir) {
  const SAMRAI::pdat::CellIndex left = shift(mid, dir, -1);
  const SAMRAI::pdat::CellIndex right = shift(mid, dir, +1);
  const double dx = geom.getDx()[static_cast<int>(dir)];
  return (q(left) - q(right)) / (2 * dx);
}
} // namespace

void GradientDetector::fillTagsForRefinement(const SAMRAI::hier::Patch& patch,
                                             int tag_id,
                                             double /* time_point */) const {
  SAMRAI::pdat::CellData<int>& tags =
      *std::static_pointer_cast<SAMRAI::pdat::CellData<int>>(
          patch.getPatchData(tag_id));

  auto threshhold_iterator = threshholds_.begin();
  for (int id : patch_data_ids_) {
    const double threshhold = *threshhold_iterator++;
    const double threshhold2 = threshhold * threshhold;
    const SAMRAI::pdat::CellData<double>& q =
        *std::static_pointer_cast<SAMRAI::pdat::CellData<double>>(
            patch.getPatchData(id));
    const SAMRAI::geom::CartesianPatchGeometry& geom =
        *GetCartesianPatchGeometry(patch);
    // Iterator through all indices in the box and compute the local gradient.
    const SAMRAI::hier::Box& box = patch.getBox();
    for (const SAMRAI::hier::Index& index : box) {
      SAMRAI::pdat::CellIndex mid(index);
      const double dq_x = computePartialGradient_(q, mid, geom, Direction::X);
      const double dq_y = computePartialGradient_(q, mid, geom, Direction::Y);
      const double dq_z = computePartialGradient_(q, mid, geom, Direction::Z);
      const double dq_norm2 = (dq_x * dq_x + dq_y * dq_y + dq_z * dq_z);
      const int tag_this_cell = (dq_norm2 < threshhold2);
      tags(mid) = tag_this_cell;
    }
  }
}
} // namespace fub