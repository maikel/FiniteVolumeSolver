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

#include "fub/AMReX/FillCutCellData.hpp"

namespace fub::amrex::cutcell {
namespace {
template <int Rank>
void FillCutCellData_(
    std::array<PatchDataView<double, Rank>, Rank> unshielded,
    std::array<PatchDataView<double, Rank>, Rank> shielded_left,
    std::array<PatchDataView<double, Rank>, Rank> shielded_right,
    std::array<PatchDataView<double, Rank>, Rank> doubly_shielded,
    std::array<PatchDataView<double, Rank>, Rank> unshielded_rel,
    std::array<PatchDataView<double, Rank>, Rank> shielded_left_rel,
    std::array<PatchDataView<double, Rank>, Rank> shielded_right_rel,
    std::array<PatchDataView<double, Rank>, Rank> doubly_shielded_rel,
    const CutCellData<Rank>& geom) {
  for (int d = 0; d < Rank; ++d) {
    const Direction dir = static_cast<Direction>(d);
    const IndexBox<Rank> indices =
        Shrink(geom.face_fractions[d].Box(), dir, {1, 1});
    ForEachIndex(indices, [&](auto... is) {
      const std::array<std::ptrdiff_t, Rank> mid{is...};
      const std::array<std::ptrdiff_t, Rank> left = Shift(mid, dir, -1);
      const std::array<std::ptrdiff_t, Rank> right = Shift(mid, dir, +1);
      const double beta_left = geom.face_fractions[d](left);
      const double beta_mid = geom.face_fractions[d](mid);
      const double beta_right = geom.face_fractions[d](right);
      FUB_ASSERT(beta_mid >= 0.0);
      const double dBetaL = std::max(0.0, beta_mid - beta_left);
      const double dBetaR = std::max(0.0, beta_mid - beta_right);
      unshielded[d](mid) = std::max(0.0, beta_mid - std::max(dBetaL, dBetaR));
      unshielded_rel[d](mid) =
          beta_mid > 0.0 ? std::clamp(unshielded[d](mid) / beta_mid, 0.0, 1.0)
                         : 0.0;
      doubly_shielded[d](mid) = std::min(dBetaL, dBetaR);
      doubly_shielded_rel[d](mid) =
          beta_mid > 0.0
              ? std::clamp(doubly_shielded[d](mid) / beta_mid, 0.0, 1.0)
              : 0.0;
      shielded_left[d](mid) = std::max(0.0, dBetaL - doubly_shielded[d](mid));
      shielded_left_rel[d](mid) =
          beta_mid > 0.0
              ? std::clamp(shielded_left[d](mid) / beta_mid, 0.0, 1.0)
              : 0.0;
      shielded_right[d](mid) = std::max(0.0, dBetaR - doubly_shielded[d](mid));
      shielded_right_rel[d](mid) =
          beta_mid > 0.0
              ? std::clamp(shielded_right[d](mid) / beta_mid, 0.0, 1.0)
              : 0.0;
    });
  }
}

} // namespace

void FillCutCellData(
    std::array<PatchDataView<double, 2>, 2> unshielded,
    std::array<PatchDataView<double, 2>, 2> shielded_left,
    std::array<PatchDataView<double, 2>, 2> shielded_right,
    std::array<PatchDataView<double, 2>, 2> doubly_shielded,
    std::array<PatchDataView<double, 2>, 2> unshielded_rel,
    std::array<PatchDataView<double, 2>, 2> shielded_left_rel,
    std::array<PatchDataView<double, 2>, 2> shielded_right_rel,
    std::array<PatchDataView<double, 2>, 2> doubly_shielded_rel,
    const CutCellData<2>& geom) {
  FillCutCellData_<2>(unshielded, shielded_left, shielded_right,
                      doubly_shielded, unshielded_rel, shielded_left_rel,
                      shielded_right_rel, doubly_shielded_rel, geom);
}

void FillCutCellData(
    std::array<PatchDataView<double, 3>, 3> unshielded,
    std::array<PatchDataView<double, 3>, 3> shielded_left,
    std::array<PatchDataView<double, 3>, 3> shielded_right,
    std::array<PatchDataView<double, 3>, 3> doubly_shielded,
    std::array<PatchDataView<double, 3>, 3> unshielded_rel,
    std::array<PatchDataView<double, 3>, 3> shielded_left_rel,
    std::array<PatchDataView<double, 3>, 3> shielded_right_rel,
    std::array<PatchDataView<double, 3>, 3> doubly_shielded_rel,
    const CutCellData<3>& geom) {
  FillCutCellData_<3>(unshielded, shielded_left, shielded_right,
                      doubly_shielded, unshielded_rel, shielded_left_rel,
                      shielded_right_rel, doubly_shielded_rel, geom);
}

} // namespace fub::amrex::cutcell
