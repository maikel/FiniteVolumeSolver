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

#ifndef FUB_AMREX_FLUX_METHOD_HPP
#define FUB_AMREX_FLUX_METHOD_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/ForEachFab.hpp"

#include "fub/Execution.hpp"
#include "fub/State.hpp"

#include <AMReX_MultiFab.H>

namespace fub::amrex {

template <typename Tag, typename FM> struct FluxMethod {
  using Equation = std::decay_t<decltype(std::declval<FM&>().GetEquation())>;
  static const int Rank = Equation::Rank();

  FluxMethod(Tag, const FM& fm);
  // FluxMethod(Tag, FM&& fm) noexcept;

  Duration ComputeStableDt(IntegratorContext& context, int level,
                           Direction dir);

  void ComputeNumericFluxes(IntegratorContext& context,
                            int level, Duration dt, Direction dir);

  int GetStencilWidth() const;

  Local<Tag, FM> flux_method_{};
};

template <typename Tag, typename FM>
FluxMethod<Tag, FM>::FluxMethod(Tag, const FM& fm) : flux_method_{fm} {}

// template <typename Tag, typename FM>
// FluxMethod<Tag, FM>::FluxMethod(Tag, FM&& fm) noexcept
//     : flux_method_{std::move(fm)} {}

template <typename Tag, typename FM>
Duration
FluxMethod<Tag, FM>::ComputeStableDt(IntegratorContext& context,
                                     int level, Direction dir) {
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  const double dx = geom.CellSize(int(dir));
  double min_dt = std::numeric_limits<double>::infinity();
  const ::amrex::MultiFab& scratch = context.GetScratch(level, dir);
  if constexpr (std::is_base_of<execution::OpenMpTag, Tag>::value) {
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel reduction(min : min_dt)
#endif
    for (::amrex::MFIter mfi(scratch,
                             ::amrex::IntVect(AMREX_D_DECL(1024000, 8, 8)));
         mfi.isValid(); ++mfi) {
      const ::amrex::Box& box = mfi.tilebox();
      Equation& equation = flux_method_->GetEquation();
      const IndexBox<Rank> cells = AsIndexBox<Rank>(box);
      const ::amrex::MultiFab& scratch = context.GetScratch(level, dir);
      const ::amrex::FArrayBox& data = scratch[mfi];
      View<const Complete<Equation>> states =
          MakeView<const Complete<Equation>>(data, equation, cells);
      const Duration dt(flux_method_->ComputeStableDt(Tag(), states, dx, dir));
      min_dt = std::min(min_dt, dt.count());
    }
    return Duration(min_dt);
  } else {
    for (::amrex::MFIter mfi(scratch); mfi.isValid(); ++mfi) {
      const ::amrex::Box& box = mfi.tilebox();
      Equation& equation = flux_method_->GetEquation();
      const IndexBox<Rank> cells = AsIndexBox<Rank>(box);
      const ::amrex::MultiFab& scratch = context.GetScratch(level, dir);
      const ::amrex::FArrayBox& data = scratch[mfi];
      View<const Complete<Equation>> states =
          MakeView<const Complete<Equation>>(data, equation, cells);
      const Duration dt(flux_method_->ComputeStableDt(Tag(), states, dx, dir));
      min_dt = std::min(min_dt, dt.count());
    }
    return Duration(min_dt);
  }
}

template <typename T, typename... Args>
using ComputeNumericFluxes_t =
    decltype(std::declval<T>().ComputeNumericFluxes(std::declval<Args>()...));

template <typename Tag, typename FM>
void FluxMethod<Tag, FM>::ComputeNumericFluxes(
    IntegratorContext& context, int level, Duration dt,
    Direction dir) {
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  ::amrex::MultiFab& fluxes = context.GetFluxes(level, dir);
  const ::amrex::MultiFab& scratch = context.GetScratch(level, dir);
  const double dx = geom.CellSize(int(dir));
  ForEachFab(Tag(), scratch, [&](::amrex::MFIter& mfi) {
    // Get a view of all complete state variables
    const ::amrex::Box& box = mfi.tilebox();
    Equation& equation = flux_method_->GetEquation();
    const int gcw = flux_method_->GetStencilWidth() + 1;
    const IndexBox<Rank> cells = Grow(AsIndexBox<Rank>(box), dir, {gcw, gcw});
    const ::amrex::FArrayBox& data = scratch[mfi];
    View<const Complete<Equation>> states =
        MakeView<const Complete<Equation>>(data, equation, cells);
    // Get a view of all flux variables
    const IndexBox<Rank> faces =
        Grow(AsIndexBox<Rank>(::amrex::surroundingNodes(box, int(dir))), dir,
             {1, 1});
    ::amrex::FArrayBox& fdata = fluxes[mfi];
    View<Conservative<Equation>> flux =
        MakeView<Conservative<Equation>>(fdata, equation, faces);
    // Pass views to implementation
    flux_method_->ComputeNumericFluxes(Tag(), flux, states, dt, dx, dir);
  });
}

template <typename Tag, typename FM>
int FluxMethod<Tag, FM>::GetStencilWidth() const {
  return flux_method_->GetStencilWidth();
}

} // namespace fub::amrex

#endif
