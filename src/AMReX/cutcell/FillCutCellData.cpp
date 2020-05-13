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

#include "fub/StateRow.hpp"
#include "fub/ext/Vc.hpp"

namespace fub::amrex::cutcell {
namespace {
void FillCutCellData__(const span<double>& us, const span<double>& ssL,
                       const span<double>& ssR, const span<double>& ds,
                       const span<double>& us_rel, const span<double>& ssL_rel,
                       const span<double>& ssR_rel, const span<double>& ds_rel,
                       const span<const double>& beta_left,
                       const span<const double>& beta_mid,
                       const span<const double>& beta_right) {
  const auto n = static_cast<int>(us.size());
  int i = 0;
  const auto pack_size = static_cast<int>(Vc::double_v::size());
  for (; i + pack_size <= n; i += pack_size) {
    const Vc::double_v betaL(&beta_left[i], Vc::Unaligned);
    Vc::double_v betaM(&beta_mid[i], Vc::Unaligned);
    const Vc::double_v betaR(&beta_right[i], Vc::Unaligned);
    const Vc::double_v dBetaL = max(Vc::double_v(0.0), betaM - betaL);
    const Vc::double_v dBetaR = max(Vc::double_v(0.0), betaM - betaR);
    const Vc::double_v unshielded =
        max(Vc::double_v(0.0), betaM - max(dBetaL, dBetaR));
    const Vc::double_v doubly_shielded = min(dBetaL, dBetaR);
    const Vc::double_v shielded_left = dBetaL - doubly_shielded;
    const Vc::double_v shielded_right = dBetaR - doubly_shielded;
    Vc::double_v unshielded_rel(0);
    Vc::double_v doubly_shielded_rel(0);
    Vc::double_v shielded_left_rel(0);
    Vc::double_v shielded_right_rel(0);
    auto non_positive_beta = betaM <= 0;
    where(non_positive_beta, betaM) = Vc::double_v(1.0);
    unshielded_rel = unshielded / betaM;
    doubly_shielded_rel = doubly_shielded / betaM;
    shielded_left_rel = shielded_left / betaM;
    shielded_right_rel = shielded_right / betaM;
    unshielded.store(&us[i], Vc::Unaligned);
    unshielded_rel.store(&us_rel[i], Vc::Unaligned);
    doubly_shielded.store(&ds[i], Vc::Unaligned);
    doubly_shielded_rel.store(&ds_rel[i], Vc::Unaligned);
    shielded_left.store(&ssL[i], Vc::Unaligned);
    shielded_left_rel.store(&ssL_rel[i], Vc::Unaligned);
    shielded_right.store(&ssR[i], Vc::Unaligned);
    shielded_right_rel.store(&ssR_rel[i], Vc::Unaligned);
  }
  for (; i < n; ++i) {
    const double betaL = beta_left[i];
    double betaM = beta_mid[i];
    const double betaR = beta_right[i];
    FUB_ASSERT(betaM >= 0.0);
    const double dBetaL = std::max(0.0, betaM - betaL);
    const double dBetaR = std::max(0.0, betaM - betaR);
    const double unshielded = std::max(0.0, betaM - std::max(dBetaL, dBetaR));
    const double doubly_shielded = std::min(dBetaL, dBetaR);
    const double shielded_left = dBetaL - doubly_shielded;
    const double shielded_right = dBetaR - doubly_shielded;

    betaM = betaM > 0 ? betaM : 1.0;
    const double unshielded_rel = unshielded / betaM;
    const double doubly_shielded_rel = doubly_shielded / betaM;
    const double shielded_left_rel = shielded_left / betaM;
    const double shielded_right_rel = shielded_right / betaM;
    us[i] = unshielded;
    us_rel[i] = unshielded_rel;
    ds[i] = doubly_shielded;
    ds_rel[i] = doubly_shielded_rel;
    ssL[i] = shielded_left;
    ssL_rel[i] = shielded_left_rel;
    ssR[i] = shielded_right;
    ssR_rel[i] = shielded_right_rel;
  }
}

template <int Rank, typename Layout>
void FillCutCellData_(
    const PatchDataView<double, Rank, Layout>& unshielded,
    const PatchDataView<double, Rank, Layout>& shielded_left,
    const PatchDataView<double, Rank, Layout>& shielded_right,
    const PatchDataView<double, Rank, Layout>& doubly_shielded,
    const PatchDataView<double, Rank, Layout>& unshielded_rel,
    const PatchDataView<double, Rank, Layout>& shielded_left_rel,
    const PatchDataView<double, Rank, Layout>& shielded_right_rel,
    const PatchDataView<double, Rank, Layout>& doubly_shielded_rel,
    const PatchDataView<const double, Rank>& beta, Direction dir) {
  IndexBox<Rank> box = unshielded.Box();
  IndexBox<Rank> betaL_box = Grow(box, dir, {1, -1});
  IndexBox<Rank> betaR_box = Grow(box, dir, {-1, 1});
  PatchDataView<const double, Rank, Layout> betaL = beta.Subview(betaL_box);
  PatchDataView<const double, Rank, Layout> betaM = beta.Subview(box);
  PatchDataView<const double, Rank, Layout> betaR = beta.Subview(betaR_box);
  ForEachRow(std::tuple{unshielded, shielded_left, shielded_right,
                        doubly_shielded, unshielded_rel, shielded_left_rel,
                        shielded_right_rel, doubly_shielded_rel, betaL, betaM,
                        betaR},
             &FillCutCellData__);
}

} // namespace

void FillCutCellData(const StridedDataView<double, 2>& unshielded,
                     const StridedDataView<double, 2>& shielded_left,
                     const StridedDataView<double, 2>& shielded_right,
                     const StridedDataView<double, 2>& doubly_shielded,
                     const StridedDataView<double, 2>& unshielded_rel,
                     const StridedDataView<double, 2>& shielded_left_rel,
                     const StridedDataView<double, 2>& shielded_right_rel,
                     const StridedDataView<double, 2>& doubly_shielded_rel,
                     const PatchDataView<const double, 2>& beta,
                     Direction dir) {
  FillCutCellData_(unshielded, shielded_left, shielded_right, doubly_shielded,
                   unshielded_rel, shielded_left_rel, shielded_right_rel,
                   doubly_shielded_rel, beta, dir);
}

void FillCutCellData(const StridedDataView<double, 3>& unshielded,
                     const StridedDataView<double, 3>& shielded_left,
                     const StridedDataView<double, 3>& shielded_right,
                     const StridedDataView<double, 3>& doubly_shielded,
                     const StridedDataView<double, 3>& unshielded_rel,
                     const StridedDataView<double, 3>& shielded_left_rel,
                     const StridedDataView<double, 3>& shielded_right_rel,
                     const StridedDataView<double, 3>& doubly_shielded_rel,
                     const PatchDataView<const double, 3>& beta,
                     Direction dir) {
  FillCutCellData_(unshielded, shielded_left, shielded_right, doubly_shielded,
                   unshielded_rel, shielded_left_rel, shielded_right_rel,
                   doubly_shielded_rel, beta, dir);
}

} // namespace fub::amrex::cutcell
