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

#include "fub/AMReX/cutcell/tagging/TagCutCells.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/Execution.hpp"
#include "fub/tagging/TagCutCells.hpp"

namespace fub::amrex::cutcell {

void TagCutCells::TagCellsForRefinement(::amrex::TagBoxArray& tags_array,
                                        Duration, int level,
                                        const GriddingAlgorithm& gridding) const
    noexcept {
  const ::amrex::MultiFab& volumes =
      gridding.GetPatchHierarchy().GetPatchLevel(level).factory->getVolFrac();
  ForEachFab(execution::openmp, tags_array, [&](const ::amrex::MFIter& mfi) {
    const IndexBox<AMREX_SPACEDIM> box =
        AsIndexBox<AMREX_SPACEDIM>(mfi.tilebox());
    PatchDataView<char, AMREX_SPACEDIM, layout_stride> tags =
        MakePatchDataView(tags_array[mfi], 0).Subview(box);
    PatchDataView<const double, AMREX_SPACEDIM, layout_stride> vols =
        MakePatchDataView(volumes[mfi], 0).Subview(box);
    ::fub::TagCutCells().TagCellsForRefinement(tags, vols);
  });
}

} // namespace fub::amrex::cutcell
