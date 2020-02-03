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

#include "fub/cutcell_method/KbnStabilisation.hpp"
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
  S boundaryL;
  S boundaryR;
};

template <typename T> struct CutCellGeometry {
  T betaUS;
  T betaL;
  T betaR;
  T centerL;
  T centerR;
};

void ComputeStableFluxes_Row(const Fluxes<double*, const double*>& fluxes,
                             const CutCellGeometry<const double*>& geom, int n,
                             Duration /* dt */, double /* dx */) {
  int face = 0;
  const int simd_width = Vc::double_v::size();
  for (face = 0; face + simd_width <= n; face += simd_width) {
    const Vc::double_v centerL(geom.centerL + face, Vc::Unaligned);
    const Vc::double_v centerR(geom.centerR + face, Vc::Unaligned);
    const Vc::double_v dL = Vc::double_v(0.5) - centerL;
    const Vc::double_v dR = Vc::double_v(0.5) + centerR;

    const Vc::double_v f(fluxes.regular + face, Vc::Unaligned);
    const Vc::double_v fbL(fluxes.boundaryL + face, Vc::Unaligned);
    const Vc::double_v fbR(fluxes.boundaryR + face, Vc::Unaligned);
    const Vc::double_v fsL = dL * f + (Vc::double_v(1.0) - dL) * fbL;
    const Vc::double_v fsR = dR * f + (Vc::double_v(1.0) - dR) * fbR;
    
    const Vc::double_v betaL(geom.betaL + face, Vc::Unaligned);
    const Vc::double_v betaR(geom.betaR + face, Vc::Unaligned);
    const Vc::double_v betaUS(geom.betaUS + face, Vc::Unaligned);
    const Vc::double_v f_stable = betaUS * f + betaL * fsL + betaR * fsR;

    fsL.store(fluxes.shielded_left + face, Vc::Unaligned);
    fsR.store(fluxes.shielded_right + face, Vc::Unaligned);
    f_stable.store(fluxes.stable + face, Vc::Unaligned);
  }

  for (; face < n; ++face) {
    const double dL = 0.5 - geom.centerL[face];
    const double dR = 0.5 + geom.centerR[face];
    const double f = fluxes.regular[face];
    const double fsL = dL * f + (1.0 - dL) * fluxes.boundaryL[face];
    const double fsR = dR * f + (1.0 - dR) * fluxes.boundaryR[face];
    const double betaL = geom.betaL[face];
    const double betaR = geom.betaR[face];
    const double betaUS = geom.betaUS[face];
    fluxes.shielded_left[face] = fsL;  // betaL ? fsL : 0.0;
    fluxes.shielded_right[face] = fsR; // betaR ? fsR : 0.0;
    fluxes.stable[face] = betaUS * f + betaL * fsL + betaR * fsR;
    FUB_ASSERT(!std::isnan(fluxes.shielded_left[face]));
    FUB_ASSERT(!std::isnan(fluxes.shielded_right[face]));
    FUB_ASSERT(!std::isnan(fluxes.stable[face]));
  }
}

template <int Rank>
void ComputeStableFluxComponents_View(
    const PatchDataView<double, Rank, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, Rank, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, Rank, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, Rank, layout_stride>& regular_fluxes,
    const PatchDataView<const double, Rank, layout_stride>& boundary_fluxes,
    const CutCellData<Rank>& geom, Duration dt, double dx, Direction dir) {
  IndexBox<Rank> faces = regular_fluxes.Box();
  IndexBox<Rank> cells_to_right = faces;
  IndexBox<Rank> cells_to_left = Grow(cells_to_right, dir, {1, -1});
  const int d = static_cast<int>(dir);
  PatchDataView<const double, Rank, layout_stride> boundaryL =
      boundary_fluxes.Subview(cells_to_left);
  PatchDataView<const double, Rank, layout_stride> boundaryR =
      boundary_fluxes.Subview(cells_to_right);
  PatchDataView<const double, Rank, layout_stride> betaUs =
      geom.unshielded_fractions_rel[d].Subview(faces);
  PatchDataView<const double, Rank, layout_stride> betaL =
      geom.shielded_left_fractions_rel[d].Subview(faces);
  PatchDataView<const double, Rank, layout_stride> betaR =
      geom.shielded_right_fractions_rel[d].Subview(faces);
  PatchDataView<const double, Rank, layout_stride> centerL =
      SliceLast(geom.boundary_centeroids.Subview(
          Embed<Rank + 1>(cells_to_left, {d, d + 1})));
  PatchDataView<const double, Rank, layout_stride> centerR =
      SliceLast(geom.boundary_centeroids.Subview(
          Embed<Rank + 1>(cells_to_right, {d, d + 1})));
  ForEachRow(std::tuple{stabilised_fluxes, shielded_left_fluxes,
                        shielded_right_fluxes, regular_fluxes, boundaryL,
                        boundaryR, betaUs, betaL, betaR, centerL, centerR},
             [dt, dx](span<double> fs, span<double> fsL, span<double> fsR,
                      span<const double> f, span<const double> fBL,
                      span<const double> fBR, span<const double> betaUs,
                      span<const double> betaL, span<const double> betaR,
                      span<const double> centerL, span<const double> centerR) {
               Fluxes<double*, const double*> fluxes;
               fluxes.stable = fs.data();
               fluxes.shielded_left = fsL.data();
               fluxes.shielded_right = fsR.data();
               fluxes.regular = f.data();
               fluxes.boundaryL = fBL.data();
               fluxes.boundaryR = fBR.data();
               CutCellGeometry<const double*> geom;
               geom.betaL = betaL.data();
               geom.betaR = betaR.data();
               geom.betaUS = betaUs.data();
               geom.centerR = centerR.data();
               geom.centerL = centerL.data();
               ComputeStableFluxes_Row(fluxes, geom, f.size(), dt, dx);
             });
}
} // namespace

void ComputeStableFluxComponents(
    const PatchDataView<double, 3, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, 3, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, 3, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, 3, layout_stride>& regular_fluxes,
    const PatchDataView<const double, 3, layout_stride>& boundary_fluxes,
    const CutCellData<3>& geom, Duration dt, double dx, Direction dir) {
  return ComputeStableFluxComponents_View<3>(
      stabilised_fluxes, shielded_left_fluxes, shielded_right_fluxes,
      regular_fluxes, boundary_fluxes, geom, dt, dx, dir);
}

void ComputeStableFluxComponents(
    const PatchDataView<double, 2, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, 2, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, 2, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, 2, layout_stride>& regular_fluxes,
    const PatchDataView<const double, 2, layout_stride>& boundary_fluxes,
    const CutCellData<2>& geom, Duration dt, double dx, Direction dir) {
  return ComputeStableFluxComponents_View<2>(
      stabilised_fluxes, shielded_left_fluxes, shielded_right_fluxes,
      regular_fluxes, boundary_fluxes, geom, dt, dx, dir);
}

} // namespace fub
