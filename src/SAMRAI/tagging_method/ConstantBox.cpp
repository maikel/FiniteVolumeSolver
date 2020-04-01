// Copyright (c) 2019 Maikel Nadolski
// Copyright (c) 2019 Patrick Denzler
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

#include "fub/SAMRAI/tagging/ConstantBox.hpp"

namespace fub::samrai {

ConstantBox::ConstantBox(const SAMRAI::hier::Box& box) : box_(box) {}

void ConstantBox::TagCellsForRefinement(GriddingAlgorithm& gridding, int level,
                                        int tag_id, Duration /* time_point */) {
  std::shared_ptr<SAMRAI::hier::PatchLevel> patch_level =
      gridding.GetPatchHierarchy().GetNative()->getPatchLevel(level);

  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : *patch_level) {
    SAMRAI::pdat::CellData<int>& tags =
        *static_cast<SAMRAI::pdat::CellData<int>*>( // NOLINT
            patch->getPatchData(tag_id).get());

    SAMRAI::hier::Box intersection(tags.getDim());
    box_.intersect(tags.getBox(), intersection);

    for (auto& i : intersection) {
      SAMRAI::pdat::CellIndex cell(i);
      tags(cell) = 1;
    }
  }
}

} // namespace fub::samrai