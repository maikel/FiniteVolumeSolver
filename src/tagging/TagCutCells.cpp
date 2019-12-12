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

#include "fub/tagging/TagCutCells.hpp"
#include "fub/ForEach.hpp"

namespace fub {
namespace {
template <int Rank>
void TagCellsForRefinement_(
    const PatchDataView<char, Rank, layout_stride>& tags,
    const PatchDataView<const double, Rank, layout_stride>& volume_fractions) {
  ForEachIndex(Intersect(tags.Box(), volume_fractions.Box()), [&](auto... is) {
    if (0.0 < volume_fractions(is...) && volume_fractions(is...) < 1.0) {
      tags(is...) = '\1';
    }
  });
}
} // namespace

void TagCutCells::TagCellsForRefinement(
    const PatchDataView<char, 1, layout_stride>& tags,
    const PatchDataView<const double, 1, layout_stride>& volume_fractions) {
  TagCellsForRefinement_(tags, volume_fractions);
}

void TagCutCells::TagCellsForRefinement(
    const PatchDataView<char, 2, layout_stride>& tags,
    const PatchDataView<const double, 2, layout_stride>& volume_fractions) {
  TagCellsForRefinement_(tags, volume_fractions);
}

void TagCutCells::TagCellsForRefinement(
    const PatchDataView<char, 3, layout_stride>& tags,
    const PatchDataView<const double, 3, layout_stride>& volume_fractions) {
  TagCellsForRefinement_(tags, volume_fractions);
}

} // namespace fub
