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

#include "fub/cutcell_method/MyStabilisation.hpp"
#include "fub/StateRow.hpp"
#include "fub/ext/Vc.hpp"

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>

namespace fub {
namespace {
template <typename T, typename S> struct Fluxes {
  T stable;
  T shielded_left;
  T shielded_right;
  S regular;
};

template <typename T> struct CutCellGeometry {
  T betaUS;
  T betaL;
  T betaR;
};

void ComputeStableFluxes_Row(const Fluxes<double*, const double*>& fluxes,
                             const CutCellGeometry<const double*>& geom,
                             fub::span<const double>::index_type n,
                             Duration /* dt */, double /* dx */) {
  int face = 0;
  const int simd_width = Vc::double_v::size();
  for (face = 0; face + simd_width <= n; face += simd_width) {
    const Vc::double_v f(fluxes.regular + face, Vc::Unaligned);
    const Vc::double_v fsL(fluxes.shielded_left + face, Vc::Unaligned);
    const Vc::double_v fsR(fluxes.shielded_right + face, Vc::Unaligned);
    const Vc::double_v betaL(geom.betaL + face, Vc::Unaligned);
    const Vc::double_v betaR(geom.betaR + face, Vc::Unaligned);
    const Vc::double_v betaUS(geom.betaUS + face, Vc::Unaligned);
    const Vc::double_v f_stable = betaUS * f + betaL * fsL + betaR * fsR;
    FUB_ASSERT(none_of(isnan(f_stable)));
    f_stable.store(fluxes.stable + face, Vc::Unaligned);
  }
  for (; face < n; ++face) {
    FUB_ASSERT(!std::isnan(fluxes.regular[face]));
    FUB_ASSERT(!std::isnan(fluxes.shielded_left[face]));
    FUB_ASSERT(!std::isnan(fluxes.shielded_right[face]));
    const double f = fluxes.regular[face];
    const double fsL = fluxes.shielded_left[face];
    const double fsR = fluxes.shielded_right[face];
    const double betaL = geom.betaL[face];
    const double betaR = geom.betaR[face];
    const double betaUS = geom.betaUS[face];
    fluxes.stable[face] = betaUS * f + betaL * fsL + betaR * fsR;
    FUB_ASSERT(!std::isnan(fluxes.stable[face]));
  }
}

template <int Rank>
void ComputeStableFluxComponents_View(
    const PatchDataView<double, Rank, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, Rank, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, Rank, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, Rank, layout_stride>& regular_fluxes,
    const CutCellData<Rank>& geom, Duration dt, double dx, Direction dir) {
  IndexBox<Rank> faces = regular_fluxes.Box();
  // const int d = static_cast<int>(dir);
  const std::size_t r = static_cast<std::size_t>(dir);
  PatchDataView<const double, Rank, layout_stride> betaUs =
      geom.unshielded_fractions_rel[r].Subview(faces);
  PatchDataView<const double, Rank, layout_stride> betaL =
      geom.shielded_left_fractions_rel[r].Subview(faces);
  PatchDataView<const double, Rank, layout_stride> betaR =
      geom.shielded_right_fractions_rel[r].Subview(faces);
  ForEachRow(std::tuple{stabilised_fluxes, shielded_left_fluxes,
                        shielded_right_fluxes, regular_fluxes,
                        betaUs, betaL, betaR},
             [dt, dx](span<double> fs, span<double> fsL, span<double> fsR,
                      span<const double> f, span<const double> betaUs,
                      span<const double> betaL, span<const double> betaR) {
               Fluxes<double*, const double*> fluxes;
               fluxes.stable = fs.data();
               fluxes.shielded_left = fsL.data();
               fluxes.shielded_right = fsR.data();
               fluxes.regular = f.data();
               CutCellGeometry<const double*> geom;
               geom.betaL = betaL.data();
               geom.betaR = betaR.data();
               geom.betaUS = betaUs.data();
               ComputeStableFluxes_Row(fluxes, geom, f.size(), dt, dx);
             });
}
} // namespace

void MyStab_ComputeStableFluxComponents(
    const PatchDataView<double, 3, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, 3, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, 3, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, 3, layout_stride>& regular_fluxes,
    const PatchDataView<const double, 3, layout_stride>& /* boundary_fluxes */,
    const CutCellData<3>& geom, Duration dt, double dx, Direction dir) {
  ComputeStableFluxComponents_View<3>(
      stabilised_fluxes, shielded_left_fluxes, shielded_right_fluxes,
      regular_fluxes, geom, dt, dx, dir);
}

void MyStab_ComputeStableFluxComponents(
    const PatchDataView<double, 2, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, 2, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, 2, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, 2, layout_stride>& regular_fluxes,
    const PatchDataView<const double, 2, layout_stride>& /* boundary_fluxes */,
    const CutCellData<2>& geom, Duration dt, double dx, Direction dir) {
  ComputeStableFluxComponents_View<2>(
      stabilised_fluxes, shielded_left_fluxes, shielded_right_fluxes,
      regular_fluxes, geom, dt, dx, dir);
}

} // namespace fub
