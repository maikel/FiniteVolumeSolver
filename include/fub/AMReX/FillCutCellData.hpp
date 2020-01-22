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

#ifndef FUB_AMREX_FILL_CUTCELL_DATA_HPP
#define FUB_AMREX_FILL_CUTCELL_DATA_HPP

#include "fub/CutCellData.hpp"
#include "fub/Direction.hpp"
#include "fub/Equation.hpp"
#include "fub/ForEach.hpp"
#include "fub/core/mdspan.hpp"

namespace fub::amrex::cutcell {

void FillCutCellData(
    std::array<StridedDataView<double, 2>, 2> unshielded,
    std::array<StridedDataView<double, 2>, 2> shielded_left,
    std::array<StridedDataView<double, 2>, 2> shielded_right,
    std::array<StridedDataView<double, 2>, 2> doubly_shielded,
    std::array<StridedDataView<double, 2>, 2> unshielded_rel,
    std::array<StridedDataView<double, 2>, 2> shielded_left_rel,
    std::array<StridedDataView<double, 2>, 2> shielded_right_rel,
    std::array<StridedDataView<double, 2>, 2> doubly_shielded_rel,
    const CutCellData<2>& data);

void FillCutCellData(
    std::array<StridedDataView<double, 3>, 3> unshielded,
    std::array<StridedDataView<double, 3>, 3> shielded_left,
    std::array<StridedDataView<double, 3>, 3> shielded_right,
    std::array<StridedDataView<double, 3>, 3> doubly_shielded,
    std::array<StridedDataView<double, 3>, 3> unshielded_rel,
    std::array<StridedDataView<double, 3>, 3> shielded_left_rel,
    std::array<StridedDataView<double, 3>, 3> shielded_right_rel,
    std::array<StridedDataView<double, 3>, 3> doubly_shielded_rel,
    const CutCellData<3>& data);

void FillCutCellData(
    std::array<PatchDataView<double, 2>, 2> unshielded,
    std::array<PatchDataView<double, 2>, 2> shielded_left,
    std::array<PatchDataView<double, 2>, 2> shielded_right,
    std::array<PatchDataView<double, 2>, 2> doubly_shielded,
    std::array<PatchDataView<double, 2>, 2> unshielded_rel,
    std::array<PatchDataView<double, 2>, 2> shielded_left_rel,
    std::array<PatchDataView<double, 2>, 2> shielded_right_rel,
    std::array<PatchDataView<double, 2>, 2> doubly_shielded_rel,
    const CutCellData<2>& data);

void FillCutCellData(
    std::array<PatchDataView<double, 3>, 3> unshielded,
    std::array<PatchDataView<double, 3>, 3> shielded_left,
    std::array<PatchDataView<double, 3>, 3> shielded_right,
    std::array<PatchDataView<double, 3>, 3> doubly_shielded,
    std::array<PatchDataView<double, 3>, 3> unshielded_rel,
    std::array<PatchDataView<double, 3>, 3> shielded_left_rel,
    std::array<PatchDataView<double, 3>, 3> shielded_right_rel,
    std::array<PatchDataView<double, 3>, 3> doubly_shielded_rel,
    const CutCellData<3>& data);

} // namespace fub::amrex::cutcell

#endif
