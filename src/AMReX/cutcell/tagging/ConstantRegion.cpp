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

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/cutcell/tagging/ConstantRegion.hpp"

namespace fub::amrex::cutcell {

ConstantBox::ConstantBox(const ::amrex::Box& coarse_box)
    : coarse_region_{coarse_box} {}

namespace {
::amrex::Box RefineToLevel(const ::amrex::Box& box, int level,
                           const GriddingAlgorithm& gridding) {
  ::amrex::Box result = box;
  for (int i = 1; i <= level; ++i) {
    result.refine(gridding.GetPatchHierarchy().GetRatioToCoarserLevel(level));
  }
  return result;
}
} // namespace

void ConstantBox::TagCellsForRefinement(::amrex::TagBoxArray& tags, Duration,
                                        int level,
                                        GriddingAlgorithm& gridding) const
    noexcept {
  const ::amrex::Box region = RefineToLevel(coarse_region_, level, gridding);
  ForEachFab(tags, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::Box where = region & mfi.tilebox();
    tags[mfi].setVal('\1', where);
  });
}

} // namespace fub::amrex::cutcell
