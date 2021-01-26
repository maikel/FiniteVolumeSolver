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

#ifndef FUB_AMREX_AXIAL_FLUX_METHOD_HPP
#define FUB_AMREX_AXIAL_FLUX_METHOD_HPP

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"

#include "fub/Execution.hpp"
#include "fub/State.hpp"

#include <AMReX_MultiFab.H>

namespace fub::amrex {

/// \ingroup FluxMethod
template <typename Tag, typename FM> struct AxialFluxMethodAdapter {
  using Equation = std::decay_t<decltype(std::declval<FM&>().GetEquation())>;
  static const int Rank = Equation::Rank();
  static const int StencilSize = FM::StencilSize;

  AxialFluxMethodAdapter(const FM& fm) : AxialFluxMethodAdapter(Tag(), fm) {}
  AxialFluxMethodAdapter(Tag, const FM& fm);

  Duration ComputeStableDt(IntegratorContext& context, int level,
                           Direction dir);

  void ComputeNumericFluxes(IntegratorContext& context, int level, Duration dt,
                            Direction dir);

  void
  ComputeNumericFluxes(FM& flux_method,
                       const View<Conservative<Equation>>& fluxes,
                       const StridedDataView<double, AMREX_SPACEDIM>& pressures,
                       const View<const Complete<Equation>>& states,
                       Duration dt, double dx, Direction dir);

  int GetStencilWidth() const;

  void AllocatePressureIfNeeded(const IntegratorContext& context, int level);

  const std::shared_ptr<std::vector<::amrex::MultiFab>>&
  SharedInterfacePressure() const noexcept;

  Local<Tag, FM> flux_method_{};
  std::shared_ptr<std::vector<::amrex::MultiFab>> pressure_;
};

template <typename FM>
AxialFluxMethodAdapter(const FM& fm)
    ->AxialFluxMethodAdapter<execution::OpenMpSimdTag, FM>;

template <typename Tag, typename FM>
AxialFluxMethodAdapter(Tag, const FM& fm)->AxialFluxMethodAdapter<Tag, FM>;

template <typename Tag, typename FM>
AxialFluxMethodAdapter<Tag, FM>::AxialFluxMethodAdapter(Tag, const FM& fm)
    : flux_method_{fm},
      pressure_{std::make_shared<std::vector<::amrex::MultiFab>>()} {}

template <typename T, typename... Args>
using ComputeStableDt_t =
    decltype(std::declval<T>().ComputeStableDt(std::declval<Args>()...));

template <typename Tag, typename FM>
Duration
AxialFluxMethodAdapter<Tag, FM>::ComputeStableDt(IntegratorContext& context,
                                                 int level, Direction dir) {
  if constexpr (is_detected<ComputeStableDt_t, FM&, IntegratorContext&, int,
                            Direction>::value) {
    return flux_method_->ComputeStableDt(context, level, dir);
  } else {
    const ::amrex::Geometry& geom = context.GetGeometry(level);
    const double dx = geom.CellSize(int(dir));
    double min_dt = std::numeric_limits<double>::max();
    const ::amrex::MultiFab& scratch = context.GetScratch(level);
    FUB_ASSERT(!scratch.contains_nan());
    if constexpr (std::is_base_of<execution::OpenMpTag, Tag>::value) {
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel reduction(min : min_dt)
#endif
      for (::amrex::MFIter mfi(scratch,
                               ::amrex::IntVect(AMREX_D_DECL(1024000, 8, 8)));
           mfi.isValid(); ++mfi) {
        const ::amrex::Box cell_box = mfi.growntilebox();
        auto&& equation = flux_method_->GetEquation();
        View<const Complete<Equation>> states =
            MakeView<const Complete<Equation>>(scratch[mfi], equation,
                                               cell_box);
        const Duration dt(
            flux_method_->ComputeStableDt(Tag(), states, dx, dir));
        min_dt = std::min(min_dt, dt.count());
      }
      double local_count = min_dt;
      return Duration(local_count);
    } else {
      for (::amrex::MFIter mfi(scratch); mfi.isValid(); ++mfi) {
        const ::amrex::Box cell_box = mfi.growntilebox();
        auto&& equation = flux_method_->GetEquation();
        View<const Complete<Equation>> states =
            MakeView<const Complete<Equation>>(scratch[mfi], equation,
                                               cell_box);
        const Duration dt(
            flux_method_->ComputeStableDt(Tag(), states, dx, dir));
        min_dt = std::min(min_dt, dt.count());
      }
      double local_count = min_dt;
      return Duration(local_count);
    }
  }
}

template <typename Tag, typename FM>
void AxialFluxMethodAdapter<Tag, FM>::AllocatePressureIfNeeded(
    const IntegratorContext& context, int level) {
  if (pressure_->size() <= static_cast<std::size_t>(level)) {
    int nlevels =
        context.GetGriddingAlgorithm()->GetPatchHierarchy().GetNumberOfLevels();
    pressure_->resize(static_cast<std::size_t>(nlevels));
  }
  const ::amrex::MultiFab& fluxes = context.GetFluxes(level, Direction::X);
  ::amrex::MultiFab& pressure = (*pressure_)[level];
  if (pressure.boxArray() != fluxes.boxArray() ||
      pressure.DistributionMap() != fluxes.DistributionMap() ||
      pressure.nGrowVect() != fluxes.nGrowVect()) {
    pressure.define(fluxes.boxArray(), fluxes.DistributionMap(), 1,
                    fluxes.nGrowVect());
  }
}

template <typename Tag, typename FM>
void AxialFluxMethodAdapter<Tag, FM>::ComputeNumericFluxes(
    IntegratorContext& context, int level, Duration dt, Direction dir) {
  AllocatePressureIfNeeded(context, level);
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  ::amrex::MultiFab& fluxes = context.GetFluxes(level, dir);
  const ::amrex::MultiFab& scratch = context.GetScratch(level);
  const int dir_v = int(dir);
  const double dx = geom.CellSize(dir_v);
  FUB_ASSERT(!scratch.contains_nan());
  const int stencil = GetStencilWidth();
  ForEachFab(Tag(), fluxes, [&](::amrex::MFIter& mfi) {
    // Get a view of all complete state variables
    const ::amrex::Box face_tilebox = mfi.growntilebox();
    const ::amrex::Box cell_validbox = scratch[mfi].box();
    const auto [cell_box, face_box] = GetCellsAndFacesInStencilRange(
        cell_validbox, face_tilebox, stencil, dir);
    auto&& equation = flux_method_->GetEquation();
    View<const Complete<Equation>> states =
        MakeView<const Complete<Equation>>(scratch[mfi], equation, cell_box);
    View<Conservative<Equation>> flux =
        MakeView<Conservative<Equation>>(fluxes[mfi], equation, face_box);
    StridedDataView<double, AMREX_SPACEDIM> pressures =
        MakePatchDataView((*pressure_)[level][mfi], 0, face_box);
    // Pass views to implementation
    ComputeNumericFluxes(*flux_method_, flux, pressures, states, dt, dx, dir);
  });
  FUB_ASSERT(!fluxes.contains_nan());
}

template <typename Tag, typename FM>
int AxialFluxMethodAdapter<Tag, FM>::GetStencilWidth() const {
  return flux_method_->GetStencilWidth();
}

template <typename Tag, typename FM>
void AxialFluxMethodAdapter<Tag, FM>::ComputeNumericFluxes(
    FM& flux_method, const View<Conservative<Equation>>& fluxes,
    const StridedDataView<double, AMREX_SPACEDIM>& pressures,
    const View<const Complete<Equation>>& states, Duration dt, double dx,
    Direction dir) {
  using Complete = Complete<Equation>;
  using Conservative = Conservative<Equation>;
  IndexBox<Rank> fluxbox = Box<0>(fluxes);
  static constexpr int kWidth = flux_method.GetStencilWidth();
  IndexBox<Rank> cellbox = Grow(fluxbox, dir, {kWidth, kWidth - 1});
  BasicView base = Subview(states, cellbox);
  std::array<View<const Complete>, StencilSize> stencil_views;
  for (std::size_t i = 0; i < StencilSize; ++i) {
    stencil_views[i] =
        Shrink(base, dir,
               {static_cast<std::ptrdiff_t>(i),
                static_cast<std::ptrdiff_t>(StencilSize - i) - 1});
  }
  std::tuple views = std::apply(
      [&fluxes, &pressures](const auto&... vs) {
        return std::tuple{fluxes, pressures, vs...};
      },
      stencil_views);
  ForEachRow(views, [this, dt, dx, dir,
                     &flux_method](const Row<Conservative>& fluxes,
                                   span<double> pressure, auto... rows) {
    ViewPointer fit = Begin(fluxes);
    ViewPointer fend = End(fluxes);
    double* p = pressure.begin();
    std::array<ViewPointer<const Complete>, StencilSize> states{Begin(rows)...};
    int n = static_cast<int>(get<0>(fend) - get<0>(fit));
    while (n >= kDefaultChunkSize) {
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Load(flux_method.stencil_array_[i], states[i]);
      }
      const Array1d pr = flux_method.ComputeNumericFlux(
          flux_method.numeric_flux_array_, flux_method.stencil_array_, dt, dx,
          dir);
      Store(fit, flux_method.numeric_flux_array_);
      Array1d::Map(p) = pr;
      p += kDefaultChunkSize;
      Advance(fit, kDefaultChunkSize);
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Advance(states[i], kDefaultChunkSize);
      }
      n = static_cast<int>(get<0>(fend) - get<0>(fit));
    }
    for (std::size_t i = 0; i < StencilSize; ++i) {
      LoadN(flux_method.stencil_array_[i], states[i], n);
    }
    const Array1d pr =
        flux_method.ComputeNumericFlux(flux_method.numeric_flux_array_,
                                       flux_method.stencil_array_, dt, dx, dir);
    StoreN(fit, flux_method.numeric_flux_array_, n);
    for (int i = 0; i < n; ++i) {
      p[i] = pr[i];
    }
  });
}

template <typename Tag, typename FM>
const std::shared_ptr<std::vector<::amrex::MultiFab>>&
AxialFluxMethodAdapter<Tag, FM>::SharedInterfacePressure() const noexcept {
  return pressure_;
}

} // namespace fub::amrex

#endif
