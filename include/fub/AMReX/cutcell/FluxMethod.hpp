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

#ifndef FUB_AMREX_CUTCELL_FLUX_METHOD_HPP
#define FUB_AMREX_CUTCELL_FLUX_METHOD_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/ext/omp.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/cutcell/IntegratorContext.hpp"
#include "fub/AMReX/output/DebugOutput.hpp"

#include "fub/AMReX/cutcell/tagging_method/TagBuffer.hpp"

#include <AMReX_EBAmrUtil.H>

#include <memory>

namespace fub::amrex::cutcell {

inline ::amrex::Vector<std::string>
AddPrefix(const ::amrex::Vector<std::string>& names,
          const std::string& prefix) {
  ::amrex::Vector<std::string> new_names{};
  new_names.reserve(names.size());
  for (const std::string& name : names) {
    new_names.push_back(fmt::format("{}_{}", prefix, name));
  }
  return new_names;
}

template <typename Tag, typename Base> class FluxMethod {
public:
  using Equation = typename Base::Equation;

  explicit FluxMethod(Base&& fm) : FluxMethod(Tag(), std::move(fm)) {}
  explicit FluxMethod(const Base& fm) : FluxMethod(Tag(), fm) {}
  FluxMethod(Tag, const Base& fm);
  FluxMethod(Tag, Base&& fm);

  static constexpr int GetStencilWidth() noexcept;

  const Base& GetBase() const noexcept;

  void PreAdvanceHierarchy(IntegratorContext& context);

  void ComputeNumericFluxes(IntegratorContext& context, int level, Duration dt,
                            Direction dir);

  Duration ComputeStableDt(IntegratorContext& context, int level,
                           Direction dir);

private:
  Local<Tag, Base> flux_method_;
  ::amrex::MultiFab gradient_x;
  ::amrex::MultiFab gradient_y;
  ::amrex::MultiFab gradient_z;
};

template <typename F>
FluxMethod(F&&) -> FluxMethod<execution::OpenMpSimdTag, std::decay_t<F>>;

template <typename Tag, typename FM>
FluxMethod(Tag, FM&&) -> FluxMethod<Tag, std::decay_t<FM>>;

template <typename Tag, typename FM>
FluxMethod<Tag, FM>::FluxMethod(Tag, FM&& flux_method)
    : flux_method_(std::move(flux_method)) {}

template <typename Tag, typename FM>
FluxMethod<Tag, FM>::FluxMethod(Tag, const FM& flux_method)
    : flux_method_(flux_method) {}

template <typename Tag, typename FM>
constexpr int FluxMethod<Tag, FM>::GetStencilWidth() noexcept {
  return FM::GetStencilWidth();
}

template <typename Tag, typename FM>
const FM& FluxMethod<Tag, FM>::GetBase() const noexcept {
  return *flux_method_;
}

template <typename Tag, typename FM>
void FluxMethod<Tag, FM>::PreAdvanceHierarchy(IntegratorContext& context) {
  if constexpr (is_detected<detail::PreAdvanceHierarchyT, FM&,
                            View<Complete<Equation>>,
                            View<const Complete<Equation>>,
                            CutCellData<AMREX_SPACEDIM>>::value) {
    const PatchHierarchy& hierarchy = context.GetPatchHierarchy();
    for (int level = 0; level < hierarchy.GetNumberOfLevels(); ++level) {
      ::amrex::MultiFab& references = context.GetReferenceStates(level);
      context.GetGriddingAlgorithm()->FillMultiFabFromLevel(references, level);
      ForEachFab(Tag(), references, [&](const ::amrex::MFIter& mfi) {
        if (context.GetFabType(level, mfi) == ::amrex::FabType::singlevalued) {
          const Equation& eq = flux_method_->GetEquation();
          const ::amrex::Box box = mfi.growntilebox();
          View<Complete<Equation>> refs =
              MakeView<Complete<Equation>>(references[mfi], eq, box);
          View<const Complete<Equation>> states =
              MakeView<const Complete<Equation>>(references[mfi], eq, box);
          CutCellData<AMREX_SPACEDIM> geom =
              hierarchy.GetCutCellData(level, mfi);
          flux_method_->PreAdvanceHierarchy(refs, states, geom);
        }
      });
    }
  }
}

namespace detail {
template <typename T, typename... Args>
using ComputeBoundaryFluxesT =
    decltype(std::declval<T>().ComputeBoundaryFluxes(std::declval<Args>()...));

template <typename T, typename... Args>
using ComputeCutCellFluxes =
    decltype(std::declval<T>().ComputeCutCellFluxes(std::declval<Args>()...));

template <typename T, typename... Args>
using ComputeGradients =
    decltype(std::declval<T>().ComputeGradients(std::declval<Args>()...));
} // namespace detail

template <typename Tag, typename FM>
void FluxMethod<Tag, FM>::ComputeNumericFluxes(IntegratorContext& context,
                                               int level, Duration dt,
                                               Direction dir) {
  const PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  ::amrex::MultiFab& references = context.GetReferenceStates(level);
  const ::amrex::MultiFab& scratch = context.GetScratch(level);
  ::amrex::MultiCutFab& boundary_fluxes = context.GetBoundaryFluxes(level);

  ::amrex::MultiFab& fluxes = context.GetFluxes(level, dir);
  ::amrex::MultiFab& fluxes_s = context.GetStabilizedFluxes(level, dir);
  ::amrex::MultiFab& fluxes_sL = context.GetShieldedFromLeftFluxes(level, dir);
  ::amrex::MultiFab& fluxes_sR = context.GetShieldedFromRightFluxes(level, dir);
  ::amrex::MultiFab& fluxes_ds = context.GetDoublyShieldedFluxes(level, dir);

  const double dx = context.GetDx(level, dir);
  const Eigen::Matrix<double, AMREX_SPACEDIM, 1> dx_vec{AMREX_D_DECL(
      context.GetDx(level, Direction::X), context.GetDx(level, Direction::Y),
      context.GetDx(level, Direction::Z))};

  boundary_fluxes.setVal(0.0);
  if constexpr (is_detected<detail::ComputeBoundaryFluxesT, FM&,
                            View<Conservative<Equation>>,
                            View<const Complete<Equation>>,
                            View<const Complete<Equation>>,
                            CutCellData<AMREX_SPACEDIM>, Duration, double,
                            Direction>::value) {
    ForEachFab(Tag(), scratch, [&](const ::amrex::MFIter& mfi) {
      const Equation& equation = flux_method_->GetEquation();
      const ::amrex::Box box = mfi.growntilebox();
      const ::amrex::FabType type = context.GetFabType(level, mfi);
      if (type == ::amrex::FabType::singlevalued) {
        CutCellData<AMREX_SPACEDIM> geom = hierarchy.GetCutCellData(level, mfi);
        auto boundary_flux = MakeView<Conservative<Equation>>(
            boundary_fluxes[mfi], equation, box);
        auto states =
            MakeView<const Complete<Equation>>(scratch[mfi], equation, box);
        auto refs =
            MakeView<const Complete<Equation>>(references[mfi], equation, box);
        flux_method_->ComputeBoundaryFluxes(boundary_flux, states, refs, geom,
                                            dt, dx, dir);
      }
    });
  }

  if constexpr (is_detected<
                    detail::ComputeGradients, FM&, View<Conservative<Equation>>,
                    View<Conservative<Equation>>, View<Conservative<Equation>>,
                    View<const Conservative<Equation>>,
                    StridedDataView<const char, AMREX_SPACEDIM>,
                    CutCellData<AMREX_SPACEDIM>,
                    Eigen::Matrix<double, AMREX_SPACEDIM, 1>>::value) {
    const ::amrex::BoxArray ba = scratch.boxArray();
    const ::amrex::DistributionMapping dm = scratch.DistributionMap();
    const int ncons = fluxes.nComp();
    const ::amrex::IntVect ngrow = scratch.nGrowVect();
    gradient_x.define(ba, dm, ncons, ngrow);
    gradient_y.define(ba, dm, ncons, ngrow);
    gradient_z.define(ba, dm, ncons, ngrow);
    ::amrex::TagBoxArray limiter_flags(ba, dm, ngrow);
    ::amrex::TagCutCells(limiter_flags, context.GetData(level));
    TagBuffer(2).TagCellsForRefinement(limiter_flags);

    ForEachFab(Tag(), scratch, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::Box box = mfi.growntilebox();
      const ::amrex::FabType type = context.GetFabType(level, mfi);
      if (type == ::amrex::FabType::singlevalued) {
        CutCellData<AMREX_SPACEDIM> geom = hierarchy.GetCutCellData(level, mfi);
        auto flags = MakePatchDataView(limiter_flags[mfi], 0, box);
        const Equation& equation = flux_method_->GetEquation();
        auto states =
            MakeView<const Conservative<Equation>>(scratch[mfi], equation, box);
        auto grad_x =
            MakeView<Conservative<Equation>>(gradient_x[mfi], equation, box);
        auto grad_y =
            MakeView<Conservative<Equation>>(gradient_y[mfi], equation, box);
        auto grad_z =
            MakeView<Conservative<Equation>>(gradient_z[mfi], equation, box);
        flux_method_->ComputeGradients(grad_x, grad_y, grad_z, states, flags,
                                       geom, dx_vec);
      }
    });

    DebugStorage& debug = *hierarchy.GetDebugStorage();
    const ::amrex::Geometry& geom = hierarchy.GetGeometry(level);
    const Equation& equation = flux_method_->GetEquation();
    const auto names =
        VarNames<Conservative<Equation>, ::amrex::Vector<std::string>>(
            equation);
    DebugSnapshotProxy snapshot =
        debug.AddSnapshot(fmt::format("Gradients_{}", int(dir)));
    if (snapshot) {
      snapshot.SaveData(gradient_x, AddPrefix(names, "GradX_"), geom);
      snapshot.SaveData(gradient_y, AddPrefix(names, "GradY_"), geom);
      snapshot.SaveData(gradient_z, AddPrefix(names, "GradZ_"), geom);
    }
  }

  static constexpr int gcw = GetStencilWidth();

  ForEachFab(Tag(), fluxes, [&](const ::amrex::MFIter& mfi) {
    const Equation& equation = flux_method_->GetEquation();
    const ::amrex::Box face_tilebox = mfi.growntilebox();
    const ::amrex::Box cell_validbox = scratch[mfi].box();
    auto [cell_box, face_box] =
        GetCellsAndFacesInStencilRange(cell_validbox, face_tilebox, gcw, dir);
    const ::amrex::FabType type = context.GetFabType(level, mfi);
    if (type == ::amrex::FabType::singlevalued) {
      CutCellData<AMREX_SPACEDIM> geom = hierarchy.GetCutCellData(level, mfi);
      {
        auto flux =
            MakeView<Conservative<Equation>>(fluxes[mfi], equation, face_box);
        auto states = MakeView<const Complete<Equation>>(scratch[mfi], equation,
                                                         cell_box);
        flux_method_->ComputeRegularFluxes(flux, states, geom, dt, dx, dir);
      }
      auto flux_s =
          MakeView<Conservative<Equation>>(fluxes_s[mfi], equation, face_box);
      auto flux_sL =
          MakeView<Conservative<Equation>>(fluxes_sL[mfi], equation, face_box);
      auto flux_sR =
          MakeView<Conservative<Equation>>(fluxes_sR[mfi], equation, face_box);
      auto flux_ds =
          MakeView<Conservative<Equation>>(fluxes_ds[mfi], equation, face_box);
      auto flux_B = MakeView<Conservative<Equation>>(boundary_fluxes[mfi],
                                                     equation, cell_box);
      auto states =
          MakeView<const Complete<Equation>>(scratch[mfi], equation, cell_box);
      if constexpr (is_detected<detail::ComputeCutCellFluxes, FM&,
                                View<Conservative<Equation>>,
                                View<Conservative<Equation>>,
                                View<Conservative<Equation>>,
                                View<Conservative<Equation>>,
                                View<Conservative<Equation>>,
                                View<Conservative<Equation>>,
                                View<const Conservative<Equation>>,
                                View<const Conservative<Equation>>,
                                View<const Conservative<Equation>>,
                                View<const Complete<Equation>>,
                                CutCellData<AMREX_SPACEDIM>, Duration,
                                Eigen::Matrix<double, AMREX_SPACEDIM, 1>,
                                Direction>::value) {
        auto flux =
            MakeView<Conservative<Equation>>(fluxes[mfi], equation, face_box);
        auto grad_x = MakeView<const Conservative<Equation>>(
            gradient_x[mfi], equation, cell_box);
        auto grad_y = MakeView<const Conservative<Equation>>(
            gradient_y[mfi], equation, cell_box);
        auto grad_z = MakeView<const Conservative<Equation>>(
            gradient_z[mfi], equation, cell_box);
        flux_method_->ComputeCutCellFluxes(flux_s, flux_sL, flux_sR, flux_ds,
                                           flux, flux_B, grad_x, grad_y, grad_z,
                                           states, geom, dt, dx_vec, dir);
      } else {
        auto flux = MakeView<const Conservative<Equation>>(fluxes[mfi],
                                                           equation, face_box);
        flux_method_->ComputeCutCellFluxes(flux_s, flux_sL, flux_sR, flux_ds,
                                           flux_B, flux, states, geom, dt, dx,
                                           dir);
      }
    } else if (type == ::amrex::FabType::regular) {
      auto flux =
          MakeView<Conservative<Equation>>(fluxes[mfi], equation, face_box);
      auto states =
          MakeView<const Complete<Equation>>(scratch[mfi], equation, cell_box);
      flux_method_->ComputeNumericFluxes(Tag(), flux, states, dt, dx, dir);
    }
  });

  DebugStorage& debug = *hierarchy.GetDebugStorage();
  const ::amrex::Geometry& geom = hierarchy.GetGeometry(level);
  const Equation& equation = flux_method_->GetEquation();
  DebugSnapshotProxy snapshot =
      debug.AddSnapshot(fmt::format("Fluxes_{}", int(dir)));
  if (snapshot) {
    const auto names =
      VarNames<Conservative<Equation>, ::amrex::Vector<std::string>>(equation);
    snapshot.SaveData(fluxes, AddPrefix(names, "RegularFlux_"), geom);
    snapshot.SaveData(fluxes_s, AddPrefix(names, "StableFlux_"), geom);
    snapshot.SaveData(fluxes_sL, AddPrefix(names, "ShieldedFromLeftFlux_"),
                      geom);
    snapshot.SaveData(fluxes_sR, AddPrefix(names, "ShieldedFromRightFlux_"),
                      geom);
    snapshot.SaveData(boundary_fluxes.ToMultiFab(0.0, 0.0),
                      AddPrefix(names, "BoundaryFlux_"), geom);
  }
}

template <typename Tag, typename FM>
Duration FluxMethod<Tag, FM>::ComputeStableDt(IntegratorContext& context,
                                              int level, Direction dir) {
  const ::amrex::MultiFab& scratch = context.GetScratch(level);
  const ::amrex::MultiFab& fluxes = context.GetFluxes(level, dir);
  const double dx = context.GetDx(level, dir);
  Local<Tag, Duration> min_dt{
      Duration(std::numeric_limits<double>::infinity())};

  static constexpr int gcw = GetStencilWidth();

  ForEachFab(Tag(), fluxes, [&](const ::amrex::MFIter& mfi) {
    const Equation& equation = flux_method_->GetEquation();
    const ::amrex::Box face_tilebox = mfi.growntilebox();
    const ::amrex::Box cell_validbox = scratch[mfi].box();
    [[maybe_unused]] auto [cell_box, face_box] =
        GetCellsAndFacesInStencilRange(cell_validbox, face_tilebox, gcw, dir);
    auto states =
        MakeView<const Complete<Equation>>(scratch[mfi], equation, cell_box);
    const ::amrex::FabType type = context.GetFabType(level, mfi);
    if (type == ::amrex::FabType::singlevalued) {
      CutCellData<AMREX_SPACEDIM> geom =
          context.GetPatchHierarchy().GetCutCellData(level, mfi);
      *min_dt = std::min(*min_dt, Duration(flux_method_->ComputeStableDt(
                                      states, geom, dx, dir)));
    } else if (type == ::amrex::FabType::regular) {
      *min_dt = std::min(*min_dt, Duration(flux_method_->ComputeStableDt(
                                      Tag(), states, dx, dir)));
    }
  });
  return Min(min_dt);
}

} // namespace fub::amrex::cutcell

#endif
