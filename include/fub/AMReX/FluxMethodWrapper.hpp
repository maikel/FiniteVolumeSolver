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

#ifndef FUB_AMREX_FLUX_METHOD_WRAPPER_HPP
#define FUB_AMREX_FLUX_METHOD_WRAPPER_HPP

#include "fub/AMReX/FluxMethod.hpp"
#include "fub/AMReX/HyperbolicSplitIntegratorContext.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"

#include "fub/State.hpp"

#include <AMReX_MultiFab.H>

#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#include <omp.h>
#endif

#include <memory>
#include <vector>

namespace fub::amrex {

template <typename FM> struct FluxMethodWrapper : FluxMethodBase {
  FluxMethodWrapper() = default;

  FluxMethodWrapper(const FM& fm);
  FluxMethodWrapper(FM&& fm) noexcept;

  std::unique_ptr<FluxMethodBase> Clone() const override;

  Duration ComputeStableDt(HyperbolicSplitIntegratorContext& context, int level,
                           Direction dir) override;

  void ComputeNumericFluxes(HyperbolicSplitIntegratorContext& context,
                            int level, Duration dt, Direction dir) override;

  int GetStencilWidth() const override;

  OmpLocal<FM> flux_method_{};
};

template <typename FM>
FluxMethodWrapper<FM>::FluxMethodWrapper(const FM& fm) : flux_method_{fm} {}

template <typename FM>
FluxMethodWrapper<FM>::FluxMethodWrapper(FM&& fm) noexcept
    : flux_method_{std::move(fm)} {}

template <typename FM>
std::unique_ptr<FluxMethodBase> FluxMethodWrapper<FM>::Clone() const {
  return std::make_unique<FluxMethodWrapper<FM>>(flux_method_.Get());
}

template <typename T, typename... Args>
using ComputeStableDt_t =
    decltype(std::declval<T>().ComputeStableDt(std::declval<Args>()...));

template <typename FM>
Duration FluxMethodWrapper<FM>::ComputeStableDt(
    HyperbolicSplitIntegratorContext& context, int level, Direction dir) {
  if constexpr (!is_detected<ComputeStableDt_t, FM&,
                             HyperbolicSplitIntegratorContext&, int,
                             Direction>::value) {
    const ::amrex::MultiFab& scratch = context.GetScratch(level, dir);
    const ::amrex::Geometry& geom = context.GetGeometry(level);
    double min_dt = std::numeric_limits<double>::infinity();
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel reduction(min : min_dt)
#endif
    for (::amrex::MFIter mfi(scratch, true); mfi.isValid(); ++mfi) {
      const ::amrex::FArrayBox& data = scratch[mfi];
      const ::amrex::Box& box = mfi.tilebox();
      using Eq = std::decay_t<decltype(flux_method_.Get().GetEquation())>;
      static constexpr int Rank = Eq::Rank();
      Eq& equation = flux_method_->GetEquation();
      const double dx = geom.CellSize(int(dir));
      const IndexBox<Rank> cells = AsIndexBox<Rank>(box);
      View<const Complete<Eq>> states =
          MakeView<const Complete<Eq>>(data, equation, cells);
      const Duration dt(
          flux_method_->ComputeStableDt(execution::simd, states, dx, dir));
      min_dt = std::min(min_dt, dt.count());
    }
    return Duration(min_dt);
  } else {
    return flux_method_->ComputeStableDt(context, level, dir);
  }
}

template <typename T, typename... Args>
using ComputeNumericFluxes_t =
    decltype(std::declval<T>().ComputeNumericFluxes(std::declval<Args>()...));

template <typename FM>
void FluxMethodWrapper<FM>::ComputeNumericFluxes(
    HyperbolicSplitIntegratorContext& context, int level, Duration dt,
    Direction dir) {
  if constexpr (!is_detected<ComputeNumericFluxes_t, FM&,
                             HyperbolicSplitIntegratorContext&, int, Duration,
                             Direction>::value) {
    const ::amrex::MultiFab& scratch = context.GetScratch(level, dir);
    ::amrex::MultiFab& fluxes = context.GetFluxes(level, dir);
    const ::amrex::Geometry& geom = context.GetGeometry(level);
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
    for (::amrex::MFIter mfi(scratch, true); mfi.isValid(); ++mfi) {
      ::amrex::FArrayBox& fdata = fluxes[mfi];
      const ::amrex::FArrayBox& data = scratch[mfi];
      const ::amrex::Box& box = mfi.tilebox();
      using Eq = std::decay_t<decltype(flux_method_.Get().GetEquation())>;
      static const int Rank = Eq::Rank();
      Eq& equation = flux_method_->GetEquation();
      const double dx = geom.CellSize(int(dir));

      const int gcw = flux_method_->GetStencilWidth() + 1;

      const IndexBox<Rank> cells = Grow(AsIndexBox<Rank>(box), dir, {gcw, gcw});
      View<const Complete<Eq>> states =
          MakeView<const Complete<Eq>>(data, equation, cells);

      const IndexBox<Rank> faces =
          Grow(AsIndexBox<Rank>(::amrex::surroundingNodes(box, int(dir))), dir,
               {1, 1});
      View<Conservative<Eq>> flux =
          MakeView<Conservative<Eq>>(fdata, equation, faces);

      flux_method_->ComputeNumericFluxes(execution::simd, flux, states, dt, dx,
                                         dir);
    }
  } else {
    flux_method_->ComputeNumericFluxes(context, level, dir, dir);
  }
}

template <typename FM> int FluxMethodWrapper<FM>::GetStencilWidth() const {
  return flux_method_->GetStencilWidth();
}

} // namespace fub::amrex

#endif

