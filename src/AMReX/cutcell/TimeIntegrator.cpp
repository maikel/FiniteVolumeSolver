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

#include "fub/AMReX/cutcell/TimeIntegrator.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/ForEach.hpp"
#include "fub/HyperbolicPatchIntegrator.hpp"
#include "fub/StateRow.hpp"
#include "fub/ext/Vc.hpp"
#include "fub/AMReX/output/DebugOutput.hpp"

#include "fub/equations/PerfectGas.hpp"

#include <algorithm>

namespace fub::amrex::cutcell {
template <typename T> struct FluxData {
  T stabilized;
  T unshielded;
  T shielded_left;
  T shielded_right;
};

template <typename T> struct CellData {
  T alpha;
  T center;
};

template <typename T> struct FaceData {
  T beta;
  T beta_us;
  T beta_sL;
  T beta_sR;
};

template <typename T> struct CutCellFractions {
  CellData<T> cell;
  FaceData<T> faceL;
  FaceData<T> faceR;
};

Vc::double_v clamp(Vc::double_v x, Vc::double_v lo, Vc::double_v hi) noexcept {
  where(x < lo, x) = lo;
  where(hi < x, x) = hi;
  return x;
}

void UpdateConservatively_Row(double* next, const double* prev,
                              const FluxData<const double*>& fluxL,
                              const FluxData<const double*>& fluxR,
                              const double* fluxB,
                              const CutCellFractions<const double*>& cc, int n,
                              double dt_over_dx) {

  std::ptrdiff_t i = 0;
  const auto pack_size = static_cast<std::ptrdiff_t>(Vc::double_v::size());
  for (i = 0; i + pack_size <= n; i += pack_size) {
    ////////////////////////////////////////////////////////////////////////////////////
    //                beta_L == beta_R  => regular update

    Vc::double_v alpha(&cc.cell.alpha[i], Vc::Unaligned);
    if (all_of(alpha == 0.0)) {
      Vc::double_v u_next(0.0);
      u_next.store(&next[i], Vc::Unaligned);
      continue;
    }
    where(alpha == 0.0, alpha) = 1.0;
    const Vc::double_v betaL(&cc.faceL.beta[i], Vc::Unaligned);
    const Vc::double_v betaR(&cc.faceR.beta[i], Vc::Unaligned);
    const Vc::double_v cell_center(&cc.cell.center[i], Vc::Unaligned);

    const Vc::double_v fL(&fluxL.stabilized[i], Vc::Unaligned);
    const Vc::double_v fR(&fluxR.stabilized[i], Vc::Unaligned);
    const Vc::double_v fB(&fluxB[i], Vc::Unaligned);

    const Vc::double_v u_prev(&prev[i], Vc::Unaligned);

    FUB_ASSERT(none_of(isnan(fL)));
    FUB_ASSERT(none_of(isnan(fR)));
    FUB_ASSERT(none_of(isnan(fB)));
    FUB_ASSERT(none_of(isnan(u_prev)));
    const Vc::double_v regular = u_prev + dt_over_dx * (fL - fR);

    ////////////////////////////////////////////////////////////////////////////////////
    //                 beta_L > beta_R

    const Vc::double_v left_greater = [&] {
      const Vc::double_v betaL_us(&cc.faceL.beta_us[i], Vc::Unaligned);
      const Vc::double_v betaL_us_over_alpha = betaL_us / alpha;
      const Vc::double_v betaL_sL(&cc.faceL.beta_sL[i], Vc::Unaligned);
      const Vc::double_v betaL_sL_over_alpha = betaL_sL / alpha;
      const Vc::double_v distance_to_boundary = Vc::double_v(0.5) + cell_center;
      const Vc::double_v betaL_sR(&cc.faceL.beta_sR[i], Vc::Unaligned);
      const Vc::double_v dBeta_over_alpha =
          clamp((betaL_sR * distance_to_boundary) / alpha, Vc::double_v(0.0),
                Vc::double_v(1.0));
      const Vc::double_v fLus(&fluxL.unshielded[i], Vc::Unaligned);
      const Vc::double_v fLssL(&fluxL.shielded_left[i], Vc::Unaligned);
      const Vc::double_v dF_us = fLus - fR;
      const Vc::double_v dF_sL = fLssL - fR;
      const Vc::double_v dB = fLus - fB;
      const Vc::double_v next = u_prev +
                                dt_over_dx * betaL_us_over_alpha * dF_us +
                                dt_over_dx * betaL_sL_over_alpha * dF_sL +
                                dt_over_dx * dBeta_over_alpha * dB;
      return next;
    }();

    ////////////////////////////////////////////////////////////////////////////////////
    //                 beta_L < beta_R

    const Vc::double_v right_greater = [&] {
      const Vc::double_v betaR_us(&cc.faceR.beta_us[i], Vc::Unaligned);
      const Vc::double_v betaR_us_over_alpha = betaR_us / alpha;
      const Vc::double_v betaR_sR(&cc.faceR.beta_sR[i], Vc::Unaligned);
      const Vc::double_v betaR_sR_over_alpha = betaR_sR / alpha;
      const Vc::double_v betaR_sL(&cc.faceR.beta_sL[i], Vc::Unaligned);
      const Vc::double_v distance_to_boundary = Vc::double_v(0.5) - cell_center;
      const Vc::double_v dBeta_over_alpha =
          clamp((betaR_sL * distance_to_boundary) / alpha, Vc::double_v(0.0),
                Vc::double_v(1.0));
      const Vc::double_v fRus(&fluxR.unshielded[i], Vc::Unaligned);
      const Vc::double_v fRssR(&fluxR.shielded_right[i], Vc::Unaligned);
      const Vc::double_v dF_us = fL - fRus;
      const Vc::double_v dF_sR = fL - fRssR;
      const Vc::double_v dB = fB - fRus;
      const Vc::double_v next = u_prev +
                                dt_over_dx * betaR_us_over_alpha * dF_us +
                                dt_over_dx * betaR_sR_over_alpha * dF_sR +
                                dt_over_dx * dBeta_over_alpha * dB;
      return next;
    }();

    Vc::double_v u_next(0.0);
    where(alpha > 0.0, u_next) = regular;
    where(alpha > 0.0 && betaL > betaR, u_next) = left_greater;
    where(alpha > 0.0 && betaL < betaR, u_next) = right_greater;
    u_next.store(&next[i], Vc::Unaligned);
  }
  for (; i < n; ++i) {
    if (cc.cell.alpha[i] == 0.0) {
      next[i] = 0.0;
      continue;
    }
    ////////////////////////////////////////////////////////////////////////////////////
    //                beta_L == beta_R  => regular update

    const double alpha = cc.cell.alpha[i];
    const double betaL = cc.faceL.beta[i];
    const double betaR = cc.faceR.beta[i];

    const double regular =
        prev[i] + dt_over_dx * (fluxL.stabilized[i] - fluxR.stabilized[i]);

    ////////////////////////////////////////////////////////////////////////////////////
    //                 beta_L > beta_R

    const double left_greater = [&] {
      const double betaL_us = cc.faceL.beta_us[i];
      const double betaL_us_over_alpha = betaL_us / alpha;
      const double betaL_sL = cc.faceL.beta_sL[i];
      const double betaL_sL_over_alpha = betaL_sL / alpha;
      const double distance_to_boundary = 0.5 + cc.cell.center[i];
      const double betaL_sR = cc.faceL.beta_sR[i];
      const double dBeta_over_alpha =
          std::clamp((betaL_sR * distance_to_boundary) / alpha, 0.0, 1.0);
      const double dF_us = fluxL.unshielded[i] - fluxR.stabilized[i];
      const double dF_sL = fluxL.shielded_left[i] - fluxR.stabilized[i];
      const double dB = fluxL.unshielded[i] - fluxB[i];
      const double next = prev[i] + dt_over_dx * betaL_us_over_alpha * dF_us +
                          dt_over_dx * betaL_sL_over_alpha * dF_sL +
                          dt_over_dx * dBeta_over_alpha * dB;
      return next;
    }();

    ////////////////////////////////////////////////////////////////////////////////////
    //                 beta_L < beta_R

    const double right_greater = [&] {
      const double betaR_us = cc.faceR.beta_us[i];
      const double betaR_us_over_alpha = betaR_us / alpha;
      const double betaR_sR = cc.faceR.beta_sR[i];
      const double betaR_sR_over_alpha = betaR_sR / alpha;
      const double betaR_sL = cc.faceR.beta_sL[i];
      const double distance_to_boundary = 0.5 - cc.cell.center[i];
      const double dBeta_over_alpha =
          std::clamp((betaR_sL * distance_to_boundary) / alpha, 0.0, 1.0);
      const double dF_us = fluxL.stabilized[i] - fluxR.unshielded[i];
      const double dF_sR = fluxL.stabilized[i] - fluxR.shielded_right[i];
      const double dB = fluxB[i] - fluxR.unshielded[i];
      const double next = prev[i] + dt_over_dx * betaR_us_over_alpha * dF_us +
                          dt_over_dx * betaR_sR_over_alpha * dF_sR +
                          dt_over_dx * dBeta_over_alpha * dB;
      return next;
    }();

    next[i] =
        betaL == betaR ? regular : betaL > betaR ? left_greater : right_greater;
  }

  for (i = 0; i < n; ++i) {
    FUB_ASSERT(cc.cell.alpha[i] == 0.0 || !std::isnan(next[i]));
  }
}

template <typename T>
using View = PatchDataView<T, AMREX_SPACEDIM, layout_stride>;

template <typename T>
View<T> MakeView(const PatchDataView<T, AMREX_SPACEDIM + 1>& pdview,
                 const IndexBox<AMREX_SPACEDIM> box, int c) {
  const IndexBox<AMREX_SPACEDIM + 1> pdbox = pdview.Box();
  const std::array<std::ptrdiff_t, 2> limits{pdbox.lower[AMREX_SPACEDIM],
                                             pdbox.upper[AMREX_SPACEDIM]};
  return SliceLast(pdview.Subview(Embed<AMREX_SPACEDIM + 1>(box, limits)), c);
}

View<const double> MakeView(const ::amrex::FArrayBox& fab,
                            const IndexBox<AMREX_SPACEDIM> box, int c) {
  return MakeView(MakePatchDataView(fab), box, c);
}
View<double> MakeView(::amrex::FArrayBox& fab,
                      const IndexBox<AMREX_SPACEDIM> box, int c) {
  return MakeView(MakePatchDataView(fab), box, c);
}

void UpdateConservatively_View(
    const View<double>& scratch, const FluxData<View<const double>>& fluxesL,
    const FluxData<View<const double>>& fluxesR,
    const View<const double>& fluxes_boundary,
    const CutCellFractions<View<const double>>& fractions, double dt_over_dx) {
  std::tuple views{scratch,
                   fluxesL.unshielded,
                   fluxesL.stabilized,
                   fluxesL.shielded_left,
                   fluxesL.shielded_right,
                   fluxesR.unshielded,
                   fluxesR.stabilized,
                   fluxesR.shielded_left,
                   fluxesR.shielded_right,
                   fluxes_boundary,
                   fractions.cell.alpha,
                   fractions.cell.center,
                   fractions.faceL.beta,
                   fractions.faceL.beta_us,
                   fractions.faceL.beta_sL,
                   fractions.faceL.beta_sR,
                   fractions.faceR.beta,
                   fractions.faceR.beta_us,
                   fractions.faceR.beta_sL,
                   fractions.faceR.beta_sR};
  ForEachRow(views,
             [&](span<double> scratch_, span<const double> fL,
                 span<const double> fLs, span<const double> fLsL,
                 span<const double> fLsR, span<const double> fR,
                 span<const double> fRs, span<const double> fRsL,
                 span<const double> fRsR, span<const double> fB,
                 span<const double> alpha, span<const double> center,
                 span<const double> betaL, span<const double> betaLus,
                 span<const double> betaLsL, span<const double> betaLsR,
                 span<const double> betaR, span<const double> betaRus,
                 span<const double> betaRsL, span<const double> betaRsR) {
               double* next = scratch_.data();
               int size = static_cast<int>(scratch_.size());
               const double* prev = next;

               FluxData<const double*> fluxesL;
               fluxesL.unshielded = fL.data();
               fluxesL.stabilized = fLs.data();
               fluxesL.shielded_left = fLsL.data();
               fluxesL.shielded_right = fLsR.data();

               FluxData<const double*> fluxesR;
               fluxesR.unshielded = fR.data();
               fluxesR.stabilized = fRs.data();
               fluxesR.shielded_left = fRsL.data();
               fluxesR.shielded_right = fRsR.data();

               const double* fluxB = fB.data();

               CutCellFractions<const double*> fractions;
               fractions.cell.alpha = alpha.data();
               fractions.cell.center = center.data();

               fractions.faceL.beta = betaL.data();
               fractions.faceL.beta_us = betaLus.data();
               fractions.faceL.beta_sL = betaLsL.data();
               fractions.faceL.beta_sR = betaLsR.data();

               fractions.faceR.beta = betaR.data();
               fractions.faceR.beta_us = betaRus.data();
               fractions.faceR.beta_sL = betaRsL.data();
               fractions.faceR.beta_sR = betaRsR.data();

               UpdateConservatively_Row(next, prev, fluxesL, fluxesR, fluxB,
                                        fractions, size, dt_over_dx);
             });
}

void TimeIntegrator::UpdateConservatively(IntegratorContext& context, int level,
                                          Duration dt, Direction dir) {
  ::amrex::MultiFab& scratch = context.GetScratch(level);
  const ::amrex::MultiFab& unshielded = context.GetFluxes(level, dir);
  const ::amrex::MultiFab& stabilized = context.GetStabilizedFluxes(level, dir);
  const ::amrex::MultiFab& shielded_left =
      context.GetShieldedFromLeftFluxes(level, dir);
  const ::amrex::MultiFab& shielded_right =
      context.GetShieldedFromRightFluxes(level, dir);
  const ::amrex::MultiCutFab& boundary = context.GetBoundaryFluxes(level);
  const PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const int n_cons = hierarchy.GetDataDescription().n_cons_components;

  const double dx = context.GetDx(level, dir);
  const double dt_over_dx = dt.count() / dx;
  const auto d = static_cast<std::size_t>(dir);
  const auto dir_v = static_cast<int>(dir);

  //  const int gcw =
  //  context.GetHyperbolicMethod().flux_method.GetStencilWidth();

  ForEachFab(execution::openmp, scratch, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FabType type = context.GetFabType(level, mfi);
    const ::amrex::Box tilebox = mfi.growntilebox();
    const ::amrex::Box all_faces_tilebox =
        ::amrex::surroundingNodes(tilebox, dir_v);
    const ::amrex::Box all_fluxes_box = unshielded[mfi].box();
    const ::amrex::Box flux_box = all_faces_tilebox & all_fluxes_box;
    const ::amrex::Box cell_box = enclosedCells(flux_box);
    const IndexBox<AMREX_SPACEDIM> cells = AsIndexBox<AMREX_SPACEDIM>(cell_box);
    const IndexBox<AMREX_SPACEDIM> faces = AsIndexBox<AMREX_SPACEDIM>(flux_box);
    if (type == ::amrex::FabType::singlevalued) {
      const IndexBox<AMREX_SPACEDIM> facesL = Grow(faces, dir, {0, -1});
      const IndexBox<AMREX_SPACEDIM> facesR = Grow(faces, dir, {-1, 0});
      for (int comp = 0; comp < n_cons; ++comp) {
        View<double> scratch_view = MakeView(scratch[mfi], cells, comp);

        FluxData<View<const double>> fluxesL;
        fluxesL.stabilized = MakeView(stabilized[mfi], facesL, comp);
        fluxesL.unshielded = MakeView(unshielded[mfi], facesL, comp);
        fluxesL.shielded_left = MakeView(shielded_left[mfi], facesL, comp);
        fluxesL.shielded_right = MakeView(shielded_right[mfi], facesL, comp);

        FluxData<View<const double>> fluxesR;
        fluxesR.stabilized = MakeView(stabilized[mfi], facesR, comp);
        fluxesR.unshielded = MakeView(unshielded[mfi], facesR, comp);
        fluxesR.shielded_left = MakeView(shielded_left[mfi], facesR, comp);
        fluxesR.shielded_right = MakeView(shielded_right[mfi], facesR, comp);

        CutCellData<AMREX_SPACEDIM> geom = hierarchy.GetCutCellData(level, mfi);
        CutCellFractions<View<const double>> fractions;
        fractions.cell.alpha = geom.volume_fractions.Subview(cells);
        fractions.cell.center =
            MakeView(geom.boundary_centeroids, cells, static_cast<int>(dir));
        fractions.faceL.beta = geom.face_fractions[d].Subview(facesL);
        fractions.faceL.beta_us = geom.unshielded_fractions[d].Subview(facesL);
        fractions.faceL.beta_sL =
            geom.shielded_left_fractions[d].Subview(facesL);
        fractions.faceL.beta_sR =
            geom.shielded_right_fractions[d].Subview(facesL);
        fractions.faceR.beta = geom.face_fractions[d].Subview(facesR);
        fractions.faceR.beta_us = geom.unshielded_fractions[d].Subview(facesR);
        fractions.faceR.beta_sL =
            geom.shielded_left_fractions[d].Subview(facesR);
        fractions.faceR.beta_sR =
            geom.shielded_right_fractions[d].Subview(facesR);

        View<const double> fluxesB = MakeView(boundary[mfi], cells, comp);

        UpdateConservatively_View(scratch_view, fluxesL, fluxesR, fluxesB,
                                  fractions, dt_over_dx);
      }
    } else if (type == ::amrex::FabType::regular) {
      ::amrex::FArrayBox& next = scratch[mfi];
      const ::amrex::FArrayBox& prev = scratch[mfi];
      const ::amrex::FArrayBox& flux = unshielded[mfi];
      const ::amrex::Box tilebox = mfi.growntilebox();
      const ::amrex::Box all_faces_tilebox =
          ::amrex::surroundingNodes(tilebox, dir_v);
      const ::amrex::Box all_fluxes_box = flux.box();
      const ::amrex::Box flux_box = all_faces_tilebox & all_fluxes_box;
      const ::amrex::Box cell_box = enclosedCells(flux_box);

      const IndexBox<AMREX_SPACEDIM + 1> cells = Embed<AMREX_SPACEDIM + 1>(
          AsIndexBox<AMREX_SPACEDIM>(cell_box), {0, n_cons});
      auto nv = MakePatchDataView(next).Subview(cells);
      auto pv = MakePatchDataView(prev).Subview(cells);
      const auto faces = Embed<AMREX_SPACEDIM + 1>(
          AsIndexBox<AMREX_SPACEDIM>(flux_box), {0, n_cons});
      auto fv = MakePatchDataView(flux).Subview(faces);
      HyperbolicPatchIntegrator<execution::OpenMpSimdTag>::UpdateConservatively(
          nv, pv, fv, dt, dx, dir);
    }
  });
}

void TimeIntegrator2::UpdateConservatively(IntegratorContext& context,
                                           int level, Duration dt,
                                           Direction dir) {
  ::amrex::MultiFab& scratch = context.GetScratch(level);
  const ::amrex::MultiFab& unshielded = context.GetFluxes(level, dir);
  const ::amrex::MultiFab& stabilized = context.GetStabilizedFluxes(level, dir);
  const ::amrex::MultiCutFab& boundary = context.GetBoundaryFluxes(level);
  const PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const int n_cons = hierarchy.GetDataDescription().n_cons_components;

  const double dx = context.GetDx(level, dir);
  const double dt_over_dx = dt.count() / dx;
  // const auto d = static_cast<std::size_t>(dir);
  const auto dir_v = static_cast<int>(dir);

  ForEachFab(execution::openmp, scratch, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FabType type = context.GetFabType(level, mfi);
    const ::amrex::Box tilebox = mfi.growntilebox();
    const ::amrex::Box all_faces_tilebox =
        ::amrex::surroundingNodes(tilebox, dir_v);
    const ::amrex::Box all_fluxes_box = unshielded[mfi].box();
    const ::amrex::Box flux_box = all_faces_tilebox & all_fluxes_box;
    const ::amrex::Box cell_box = enclosedCells(flux_box);
    const IndexBox<AMREX_SPACEDIM> cells = AsIndexBox<AMREX_SPACEDIM>(cell_box);
    const IndexBox<AMREX_SPACEDIM> faces = AsIndexBox<AMREX_SPACEDIM>(flux_box);
    if (type == ::amrex::FabType::singlevalued) {
      CutCellData<AMREX_SPACEDIM> geom = hierarchy.GetCutCellData(level, mfi);
      const IndexBox<AMREX_SPACEDIM> facesL = Grow(faces, dir, {0, -1});
      const IndexBox<AMREX_SPACEDIM> facesR = Grow(faces, dir, {-1, 0});
      for (int comp = 0; comp < n_cons; ++comp) {
        View<double> scratch_view = MakeView(scratch[mfi], cells, comp);
        View<const double> alpha = geom.volume_fractions.Subview(cells);
        View<const double> betaL = geom.face_fractions[dir_v].Subview(facesL);
        View<const double> betaR = geom.face_fractions[dir_v].Subview(facesR);
        View<const double> fluxL = MakeView(stabilized[mfi], facesL, comp);
        View<const double> fluxR = MakeView(stabilized[mfi], facesR, comp);
        View<const double> fluxB = MakeView(boundary[mfi], cells, comp);
        ForEachIndex(cells, [&](auto... is) {
          const Index<AMREX_SPACEDIM> cell{is...};
          const double volfrac = alpha(cell);
          if (volfrac > 0.0) {
            const Index<AMREX_SPACEDIM> faceL = cell;
            const Index<AMREX_SPACEDIM> faceR = Shift(cell, dir, 1);
            const double bL = betaL(faceL);
            const double bR = betaR(faceR);
            if (bL == bR) {
              FUB_ASSERT(volfrac == bL);
              const double fL = fluxL(faceL);
              const double fR = fluxR(faceR);
              const double prev = scratch_view(cell);
              const double next = prev + dt_over_dx * (fL - fR);
              scratch_view(cell) = next;
            } else {
              const double dB = bL - bR;
              const double fL = fluxL(faceL);
              const double fR = fluxR(faceR);
              const double fB = fluxB(cell);
              const double prev = scratch_view(cell);
              const double next =
                  prev + dt_over_dx / volfrac * (bL * fL - bR * fR - dB * fB);
              scratch_view(cell) = next;
            }
          }
        });
      }
    } else if (type == ::amrex::FabType::regular) {
      ::amrex::FArrayBox& next = scratch[mfi];
      const ::amrex::FArrayBox& prev = scratch[mfi];
      const ::amrex::FArrayBox& flux = unshielded[mfi];
      const ::amrex::Box tilebox = mfi.growntilebox();
      const ::amrex::Box all_faces_tilebox =
          ::amrex::surroundingNodes(tilebox, dir_v);
      const ::amrex::Box all_fluxes_box = flux.box();
      const ::amrex::Box flux_box = all_faces_tilebox & all_fluxes_box;
      const ::amrex::Box cell_box = enclosedCells(flux_box);

      const IndexBox<AMREX_SPACEDIM + 1> cells = Embed<AMREX_SPACEDIM + 1>(
          AsIndexBox<AMREX_SPACEDIM>(cell_box), {0, n_cons});
      auto nv = MakePatchDataView(next).Subview(cells);
      auto pv = MakePatchDataView(prev).Subview(cells);
      const auto faces = Embed<AMREX_SPACEDIM + 1>(
          AsIndexBox<AMREX_SPACEDIM>(flux_box), {0, n_cons});
      auto fv = MakePatchDataView(flux).Subview(faces);
      HyperbolicPatchIntegrator<execution::OpenMpSimdTag>::UpdateConservatively(
          nv, pv, fv, dt, dx, dir);
    }
  });
}

} // namespace fub::amrex::cutcell
