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

namespace fub {
namespace amrex {

void FillCutCellData(PatchDataView<double, 2> unshielded,
                     PatchDataView<double, 2> shielded_left,
                     PatchDataView<double, 2> shielded_right,
                     PatchDataView<double, 2> doubly_shielded,
                     const CutCellData<2>& data, Direction dir) {
  const IndexBox<2> indices = Shrink(data.face_fractions.Box(), dir, {1, 1});
  fub::ForEachIndex(indices, [&](std::ptrdiff_t i, std::ptrdiff_t j) {
    const std::array<std::ptrdiff_t, 2> mid{i, j};
    const std::array<std::ptrdiff_t, 2> left = Shift(mid, dir, -1);
    const std::array<std::ptrdiff_t, 2> right = Shift(mid, dir, +1);
    const double beta_left = data.face_fractions(left);
    const double beta_mid = data.face_fractions(mid);
    const double beta_right = data.face_fractions(right);
    FUB_ASSERT(beta_mid >= 0.0);
    const double dBetaL = std::max(0.0, beta_mid - beta_left);
    const double dBetaR = std::max(0.0, beta_mid - beta_right);
    unshielded(mid) = beta_mid - std::max(dBetaL, dBetaR);
    doubly_shielded(mid) = std::min(dBetaL, dBetaR);
    shielded_left(mid) = std::max(0.0, dBetaL - doubly_shielded(mid));
    shielded_right(mid) = std::max(0.0, dBetaR - doubly_shielded(mid));
  });
}

void FillCutCellData(PatchDataView<double, 3> unshielded,
                     PatchDataView<double, 3> shielded_left,
                     PatchDataView<double, 3> shielded_right,
                     PatchDataView<double, 3> doubly_shielded,
                     const CutCellData<3>& data, Direction dir) {
  const IndexBox<3> indices = Shrink(data.face_fractions.Box(), dir, {1, 1});
  fub::ForEachIndex(
      indices, [&](std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k) {
        const std::array<std::ptrdiff_t, 3> mid{i, j, k};
        const std::array<std::ptrdiff_t, 3> left = Shift(mid, dir, -1);
        const std::array<std::ptrdiff_t, 3> right = Shift(mid, dir, +1);
        const double beta_left = data.face_fractions(left);
        const double beta_mid = data.face_fractions(mid);
        const double beta_right = data.face_fractions(right);
        FUB_ASSERT(beta_mid >= 0.0);
        const double dBetaL = std::max(0.0, beta_mid - beta_left);
        const double dBetaR = std::max(0.0, beta_mid - beta_right);
        unshielded(mid) = beta_mid - std::max(dBetaL, dBetaR);
        doubly_shielded(mid) = std::min(dBetaL, dBetaR);
        shielded_left(mid) = std::max(0.0, dBetaL - doubly_shielded(mid));
        shielded_right(mid) = std::max(0.0, dBetaR - doubly_shielded(mid));
      });
}

} // namespace amrex
} // namespace fub
