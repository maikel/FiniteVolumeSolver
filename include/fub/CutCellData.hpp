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

#ifndef FUB_CUTCELL_DATA_HPP
#define FUB_CUTCELL_DATA_HPP

#include "fub/PatchDataView.hpp"
#include "fub/ext/Eigen.hpp"

namespace fub {

template <int Rank> struct CutCellData {
  // The next member variables are given by AMReX
  PatchDataView<const double, Rank> volume_fractions;
  std::array<PatchDataView<const double, Rank>, Rank> face_fractions;
  PatchDataView<const double, Rank + 1> boundary_normals;
  PatchDataView<const double, Rank + 1> boundary_centeroids;
  // The following members need to be computed from the AMReX EB database
  PatchDataView<const double, Rank + 1> distance_to_boundary_left;
  PatchDataView<const double, Rank + 1> distance_to_boundary_right;
  std::array<PatchDataView<const double, Rank>, Rank> unshielded_fractions;
  std::array<PatchDataView<const double, Rank>, Rank> shielded_left_fractions;
  std::array<PatchDataView<const double, Rank>, Rank> shielded_right_fractions;
  std::array<PatchDataView<const double, Rank>, Rank> doubly_shielded_fractions;
  std::array<PatchDataView<const double, Rank>, Rank> unshielded_fractions_rel;
  std::array<PatchDataView<const double, Rank>, Rank>
      shielded_left_fractions_rel;
  std::array<PatchDataView<const double, Rank>, Rank>
      shielded_right_fractions_rel;
  std::array<PatchDataView<const double, Rank>, Rank>
      doubly_shielded_fractions_rel;
};

[[nodiscard]] bool IsCutCell(const CutCellData<2>& geom,
                             const std::array<std::ptrdiff_t, 2>& index);

[[nodiscard]] bool IsCutCell(const CutCellData<3>& geom,
                             const std::array<std::ptrdiff_t, 3>& index);

[[nodiscard]] Eigen::Vector2d
GetBoundaryNormal(const CutCellData<2>& ccdata,
                  const std::array<std::ptrdiff_t, 2>& index);

[[nodiscard]] Eigen::Vector3d
GetBoundaryNormal(const CutCellData<3>& ccdata,
                  const std::array<std::ptrdiff_t, 3>& index);

} // namespace fub

#endif
