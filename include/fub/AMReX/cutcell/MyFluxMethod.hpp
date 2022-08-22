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

#ifndef FUB_AMREX_CUTCELL_MY_FLUX_METHOD_HPP
#define FUB_AMREX_CUTCELL_MY_FLUX_METHOD_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/ext/omp.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/cutcell/IntegratorContext.hpp"
#include "fub/AMReX/cutcell/TimeIntegrator.hpp"
#include "fub/AMReX/output/DebugOutput.hpp"

#include "fub/AMReX/cutcell/FluxMethod.hpp"
#include "fub/AMReX/cutcell/tagging_method/TagBuffer.hpp"

#include <AMReX_EBAmrUtil.H>

#include <memory>

namespace fub::amrex::cutcell {

template <typename Tag, typename Base> class MyFluxMethod {
public:
  using Equation = typename Base::Equation;

  explicit MyFluxMethod(Base&& fm) : MyFluxMethod(Tag(), std::move(fm)) {}
  explicit MyFluxMethod(const Base& fm) : MyFluxMethod(Tag(), fm) {}
  MyFluxMethod(Tag, const Base& fm);
  MyFluxMethod(Tag, Base&& fm);

  MyFluxMethod(const MyFluxMethod& other) : flux_method_(other.flux_method_), complete_from_cons_(flux_method_->GetEquation()) {}
  MyFluxMethod(MyFluxMethod&& other) = default;

  MyFluxMethod& operator=(const MyFluxMethod& other) {
    flux_method_ = other.flux_method_;
    complete_from_cons_ = other.complete_from_cons_;
  }
  MyFluxMethod& operator=(MyFluxMethod&& other) = default;

  static constexpr int GetStencilWidth() noexcept;

  const Base& GetBase() const noexcept;

  void PreSplitStep(IntegratorContext& context, int level, Duration dt,
                    Direction dir, std::pair<int, int> subcycle);

  void ComputeGradients(IntegratorContext& context,
                        const ::amrex::MultiFab& source, int level);

  void ComputeReferenceStates(IntegratorContext& context,
                              ::amrex::MultiFab& dest, const ::amrex::MultiFab& src, int level, Duration dt, Direction dir);

  void IntegrateInTime(IntegratorContext& context, ::amrex::MultiFab& dest,
                       const ::amrex::MultiFab& source, int level, Duration dt,
                       Direction dir);

  void UpdateStates(IntegratorContext& context, ::amrex::MultiFab& dest,
                    const ::amrex::MultiFab& source, int level, Duration dt,
                    Direction dir);

  void ComputeNumericFluxes(IntegratorContext& context,
                            const ::amrex::MultiFab& scratch, int level,
                            Duration dt, Direction dir);

  void ComputeNumericFluxes(IntegratorContext& context, int level, Duration dt,
                            Direction dir);

  Duration ComputeStableDt(IntegratorContext& context, int level,
                           Direction dir);

private:
  Local<Tag, Base> flux_method_;
  TimeIntegrator2 time_integrator_;
  Reconstruction<Tag, Equation> complete_from_cons_;
  std::array<::amrex::MultiFab, 3> gradients_{};
};

template <typename F>
MyFluxMethod(F&&) -> MyFluxMethod<execution::OpenMpSimdTag, std::decay_t<F>>;

template <typename Tag, typename FM>
MyFluxMethod(Tag, FM&&) -> MyFluxMethod<Tag, std::decay_t<FM>>;

template <typename Tag, typename FM>
MyFluxMethod<Tag, FM>::MyFluxMethod(Tag, FM&& flux_method)
    : flux_method_(std::move(flux_method)), complete_from_cons_(flux_method_->GetEquation()) {}

template <typename Tag, typename FM>
MyFluxMethod<Tag, FM>::MyFluxMethod(Tag, const FM& flux_method)
    : flux_method_(flux_method), complete_from_cons_(flux_method_->GetEquation()) {}

template <typename Tag, typename FM>
constexpr int MyFluxMethod<Tag, FM>::GetStencilWidth() noexcept {
  return FM::GetStencilWidth();
}

template <typename Tag, typename FM>
const FM& MyFluxMethod<Tag, FM>::GetBase() const noexcept {
  return *flux_method_;
}

template <typename Tag, typename FM>
void MyFluxMethod<Tag, FM>::PreSplitStep(IntegratorContext& context, int level,
                                         Duration dt, Direction dir,
                                         std::pair<int, int> subcycle) {
  const ::amrex::MultiFab& scratch = context.GetScratch(level);
  const ::amrex::IntVect ngrow = scratch.nGrowVect();
  ComputeGradients(context, scratch, level);
  if ((subcycle.first % AMREX_SPACEDIM) == 0) {
    ::amrex::MultiFab& references = context.GetReferenceStates(level);
    ComputeReferenceStates(context, references, scratch, level, dt, dir);
  }
}

template <typename Tag, typename FM>
void MyFluxMethod<Tag, FM>::ComputeReferenceStates(IntegratorContext& context,
                                                   ::amrex::MultiFab& dest,
                                                   const ::amrex::MultiFab& src,
                                                   int level, Duration dt, Direction dir) {
  PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const Eigen::Matrix<double, AMREX_SPACEDIM, 1> h{AMREX_D_DECL(
      context.GetDx(level, Direction::X), context.GetDx(level, Direction::Y),
      context.GetDx(level, Direction::Z))};
  const ::amrex::BoxArray ba = dest.boxArray();
  const ::amrex::DistributionMapping dm = dest.DistributionMap();
  ::amrex::MultiFab& ref_gradients_x = context.GetReferenceGradients(level, Direction::X);
  ::amrex::MultiFab& ref_gradients_y = context.GetReferenceGradients(level, Direction::Y);
  ::amrex::MultiFab& ref_gradients_z = context.GetReferenceGradients(level, Direction::Z);

  ForEachFab(Tag(), dest, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::Box box = mfi.growntilebox(src.nGrow() - 1);
    const ::amrex::FabType type = context.GetFabType(level, mfi);
    if (type == ::amrex::FabType::singlevalued) {
      CutCellData<AMREX_SPACEDIM> geom = hierarchy.GetCutCellData(level, mfi);
      const Equation& equation = flux_method_->GetEquation();
      auto refs = MakeView<Complete<Equation>>(dest[mfi], equation, box);
      auto ref_grad_x = MakeView<typename FM::Gradient>(ref_gradients_x[mfi], equation, box);
      auto ref_grad_y = MakeView<typename FM::Gradient>(ref_gradients_y[mfi], equation, box);
      auto ref_grad_z = MakeView<typename FM::Gradient>(ref_gradients_z[mfi], equation, box);
      auto states = MakeView<const Complete<Equation>>(src[mfi], equation, src[mfi].box());
      auto grad_x = MakeView<const typename FM::Gradient>(gradients_[0][mfi],
                                                          equation, gradients_[0][mfi].box());
      auto grad_y = MakeView<const typename FM::Gradient>(gradients_[1][mfi],
                                                          equation, gradients_[1][mfi].box());
      auto grad_z = MakeView<const typename FM::Gradient>(gradients_[2][mfi],
                                                          equation, gradients_[2][mfi].box());
      std::array<View<const typename FM::Gradient>, 3> grads{grad_x, grad_y,
                                                             grad_z};
      std::array<View<typename FM::Gradient>, 3> ref_grads{ref_grad_x, ref_grad_y,
                                                             ref_grad_z};
      flux_method_->ComputeBoundaryStates(refs, span{ref_grads}, states, span{grads}, geom, h, dt, dir);
    }
  });
}

template <typename Tag, typename FM>
void MyFluxMethod<Tag, FM>::ComputeGradients(IntegratorContext& context,
                                             const ::amrex::MultiFab& source,
                                             int level) {
  PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const Eigen::Matrix<double, AMREX_SPACEDIM, 1> h{AMREX_D_DECL(
      context.GetDx(level, Direction::X), context.GetDx(level, Direction::Y),
      context.GetDx(level, Direction::Z))};
  const ::amrex::BoxArray ba = source.boxArray();
  const ::amrex::DistributionMapping dm = source.DistributionMap();
  const int ncons = context.GetFluxes(level, Direction::X).nComp();
  const ::amrex::IntVect ngrow = source.nGrowVect();
  for (::amrex::MultiFab& gradient : gradients_) {
    gradient.define(ba, dm, ncons, ngrow);
    gradient.setVal(0.0);
  }

  ForEachFab(Tag(), source, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::Box box = mfi.growntilebox(source.nGrow() - 1);
    const ::amrex::FabType type = context.GetFabType(level, mfi);
    if (type == ::amrex::FabType::singlevalued) {
      CutCellData<AMREX_SPACEDIM> geom = hierarchy.GetCutCellData(level, mfi);
      const Equation& equation = flux_method_->GetEquation();
      auto states = MakeView<const Complete<Equation>>(source[mfi], equation,
                                                       source[mfi].box());
      auto grad_x =
          MakeView<typename FM::Gradient>(gradients_[0][mfi], equation, box);
      auto grad_y =
          MakeView<typename FM::Gradient>(gradients_[1][mfi], equation, box);
      auto grad_z =
          MakeView<typename FM::Gradient>(gradients_[2][mfi], equation, box);
      flux_method_->ComputeGradients(grad_x, grad_y, grad_z, states, geom, h);
    }
  });

  DebugStorage& debug = *hierarchy.GetDebugStorage();
  const ::amrex::Geometry& geom = hierarchy.GetGeometry(level);
  const Equation& equation = flux_method_->GetEquation();
  const auto names =
      VarNames<typename FM::Gradient, ::amrex::Vector<std::string>>(equation);
  DebugSnapshotProxy snapshot = debug.AddSnapshot("Gradients");
  if (snapshot) {
    snapshot.SaveData(gradients_[0], AddPrefix(names, "GradX_"), geom);
    snapshot.SaveData(gradients_[1], AddPrefix(names, "GradY_"), geom);
    snapshot.SaveData(gradients_[2], AddPrefix(names, "GradZ_"), geom);
  }
}

template <typename Tag, typename FM>
void MyFluxMethod<Tag, FM>::IntegrateInTime(IntegratorContext& context,
                                            ::amrex::MultiFab& dest,
                                            const ::amrex::MultiFab& source,
                                            int level, Duration dt,
                                            Direction dir) {
  ComputeNumericFluxes(context, source, level, dt, dir);
  time_integrator_.UpdateConservatively(dest, source, context, level, dt, dir);
  complete_from_cons_.CompleteFromCons(context, dest ,source, level, dt);
}

template <typename Tag, typename FM>
void MyFluxMethod<Tag, FM>::ComputeNumericFluxes(
    IntegratorContext& context, const ::amrex::MultiFab& scratch, int level,
    Duration dt, Direction dir) {
  const PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const Eigen::Matrix<double, AMREX_SPACEDIM, 1> h{AMREX_D_DECL(
      context.GetDx(level, Direction::X), context.GetDx(level, Direction::Y),
      context.GetDx(level, Direction::Z))};
  const ::amrex::MultiFab& references = context.GetReferenceStates(level);
  const ::amrex::MultiFab& ref_gradients_x = context.GetReferenceGradients(level, Direction::X);
  const ::amrex::MultiFab& ref_gradients_y = context.GetReferenceGradients(level, Direction::Y);
  const ::amrex::MultiFab& ref_gradients_z = context.GetReferenceGradients(level, Direction::Z);

  ::amrex::MultiCutFab& boundary_fluxes = context.GetBoundaryFluxes(level);
  ::amrex::MultiFab& fluxes = context.GetFluxes(level, dir);
  ::amrex::MultiFab& fluxes_s = context.GetStabilizedFluxes(level, dir);
  ::amrex::MultiFab& fluxes_sL = context.GetShieldedFromLeftFluxes(level, dir);
  ::amrex::MultiFab& fluxes_sR = context.GetShieldedFromRightFluxes(level, dir);
  ::amrex::MultiFab& fluxes_ds = context.GetDoublyShieldedFluxes(level, dir);

  const double dx = context.GetDx(level, dir);
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
      auto flux =
          MakeView<Conservative<Equation>>(fluxes[mfi], equation, face_box);
      auto states =
          MakeView<const Complete<Equation>>(scratch[mfi], equation, cell_box);
      auto grad_x = MakeView<const typename FM::Gradient>(gradients_[0][mfi],
                                                          equation, cell_box);
      auto grad_y = MakeView<const typename FM::Gradient>(gradients_[1][mfi],
                                                          equation, cell_box);
      auto grad_z = MakeView<const typename FM::Gradient>(gradients_[2][mfi],
                                                          equation, cell_box);
      flux_method_->ComputeRegularFluxes(flux, states, grad_x, grad_y, grad_z,
                                         geom, dt, dx, dir);

    } else if (type == ::amrex::FabType::regular) {
      auto flux =
          MakeView<Conservative<Equation>>(fluxes[mfi], equation, face_box);
      auto states =
          MakeView<const Complete<Equation>>(scratch[mfi], equation, cell_box);
      flux_method_->ComputeNumericFluxes(Tag(), flux, states, dt, dx, dir);
    }
  });

  ForEachFab(Tag(), scratch, [&](const ::amrex::MFIter& mfi) {
    const Equation& equation = flux_method_->GetEquation();
    const ::amrex::Box tilebox = mfi.growntilebox();
    const ::amrex::Box all_faces = ::amrex::surroundingNodes(tilebox, int(dir));
    const ::amrex::Box face_box = fluxes_s[mfi].box();
    const ::amrex::Box cell_box =
        ::amrex::enclosedCells(all_faces & fluxes_s[mfi].box());

    const ::amrex::FabType type = context.GetFabType(level, mfi);
    if (type == ::amrex::FabType::singlevalued) {
      CutCellData<AMREX_SPACEDIM> geom = hierarchy.GetCutCellData(level, mfi);
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
      auto flux =
          MakeView<Conservative<Equation>>(fluxes[mfi], equation, face_box);
      auto states = MakeView<const Complete<Equation>>(scratch[mfi], equation,
                                                       scratch[mfi].box());
      auto grad_x = MakeView<const typename FM::Gradient>(
          gradients_[0][mfi], equation, gradients_[0][mfi].box());
      auto grad_y = MakeView<const typename FM::Gradient>(
          gradients_[1][mfi], equation, gradients_[1][mfi].box());
      auto grad_z = MakeView<const typename FM::Gradient>(
          gradients_[2][mfi], equation, gradients_[2][mfi].box());
      auto eb_ref = MakeView<const Complete<Equation>>(
          references[mfi], equation, references[mfi].box());
      auto ref_grad_x = MakeView<const typename FM::Gradient>(
          ref_gradients_x[mfi], equation, ref_gradients_x[mfi].box());
      auto ref_grad_y = MakeView<const typename FM::Gradient>(
          ref_gradients_y[mfi], equation, ref_gradients_y[mfi].box());
      auto ref_grad_z = MakeView<const typename FM::Gradient>(
          ref_gradients_z[mfi], equation, ref_gradients_z[mfi].box());
      flux_method_->ComputeCutCellFluxes(flux_s, flux_sL, flux_sR, flux_ds,
                                         flux, flux_B, eb_ref, ref_grad_x, ref_grad_y, ref_grad_z, grad_x, grad_y,
                                         grad_z, states, geom, dt, h, dir);
    }
  });

  DebugStorage& debug = *hierarchy.GetDebugStorage();
  const ::amrex::Geometry& geom = hierarchy.GetGeometry(level);
  const Equation& equation = flux_method_->GetEquation();
  DebugSnapshotProxy snapshot =
      debug.AddSnapshot(fmt::format("Fluxes_{}", int(dir)));
  if (snapshot) {
    const auto names =
        VarNames<Conservative<Equation>, ::amrex::Vector<std::string>>(
            equation);
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
void MyFluxMethod<Tag, FM>::ComputeNumericFluxes(IntegratorContext& context,
                                                 int level, Duration dt,
                                                 Direction dir) {
  ComputeNumericFluxes(context, context.GetScratch(level), level, dt, dir);
}

template <typename Tag, typename FM>
Duration MyFluxMethod<Tag, FM>::ComputeStableDt(IntegratorContext& context,
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
    auto cell_box = GetCellsAndFacesInStencilRange(cell_validbox, face_tilebox,
                                                   gcw, dir)[0];
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
