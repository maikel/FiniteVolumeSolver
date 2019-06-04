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

#include "fub/AMReX/HyperbolicSplitIntegratorContext.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxReg_C.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>

#include <fmt/format.h>

namespace fub::amrex {

////////////////////////////////////////////////////////////////////////////////
//                                   HyperbolicSplitIntegratorContext::LevelData

////////////////////////////////////////////////////////////////////////////////
// Move Assignment Operator

HyperbolicSplitIntegratorContext::LevelData&
HyperbolicSplitIntegratorContext::LevelData::
operator=(LevelData&& other) noexcept {
  if (other.coarse_fine.fineLevel() > 0) {
    // If we do not invoke clear in beforehand it will throw an error in AMReX
    coarse_fine.clear();
    coarse_fine.define(scratch[0].boxArray(), scratch[0].DistributionMap(),
                       other.coarse_fine.refRatio(),
                       other.coarse_fine.fineLevel(),
                       other.coarse_fine.nComp());
  }
  scratch = std::move(other.scratch);
  fluxes = std::move(other.fluxes);
  time_point = other.time_point;
  regrid_time_point = other.regrid_time_point;
  cycles = other.cycles;
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
//                                              HyperbolicSplitIntegratorContext

////////////////////////////////////////////////////////////////////////////////
// Constructor and Assignment Operators

HyperbolicSplitIntegratorContext::HyperbolicSplitIntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, NumericalMethod nm)
    : ghost_cell_width_{nm.flux_method.GetStencilWidth() + 1},
      gridding_{std::move(gridding)}, data_{}, numerical_method_{
                                                   std::move(nm)} {
  data_.reserve(
      static_cast<std::size_t>(GetPatchHierarchy().GetMaxNumberOfLevels()));
  // Allocate auxiliary data arrays for each refinement level in the hierarchy
  ResetHierarchyConfiguration();
}

HyperbolicSplitIntegratorContext::HyperbolicSplitIntegratorContext(
    const HyperbolicSplitIntegratorContext& other)
    : ghost_cell_width_{other.ghost_cell_width_}, gridding_{other.gridding_},
      data_(static_cast<std::size_t>(GetPatchHierarchy().GetNumberOfLevels())),
      numerical_method_{other.numerical_method_} {
  // Allocate auxiliary data arrays
  ResetHierarchyConfiguration();
  // Copy time stamps and cycle counters
  std::size_t n_levels =
      static_cast<std::size_t>(GetPatchHierarchy().GetNumberOfLevels());
  for (std::size_t i = 0; i < n_levels; ++i) {
    data_[i].cycles = other.data_[i].cycles;
    data_[i].time_point = other.data_[i].time_point;
    data_[i].regrid_time_point = other.data_[i].regrid_time_point;
  }
}

HyperbolicSplitIntegratorContext&
HyperbolicSplitIntegratorContext::HyperbolicSplitIntegratorContext::
operator=(const HyperbolicSplitIntegratorContext& other) {
  // We use the copy and move idiom to provide the strong exception guarantee.
  // If an exception occurs we do not change the original object.
  HyperbolicSplitIntegratorContext tmp{other};
  return (*this = std::move(tmp));
}

///////////////////////////////////////////////////////////////////////////////
//                                                            Member Accessors

const BoundaryCondition&
HyperbolicSplitIntegratorContext::GetBoundaryCondition(int level) const {
  return gridding_->GetBoundaryCondition(level);
}

BoundaryCondition&
HyperbolicSplitIntegratorContext::GetBoundaryCondition(int level) {
  return gridding_->GetBoundaryCondition(level);
}

const FluxMethod& HyperbolicSplitIntegratorContext::GetFluxMethod() const
    noexcept {
  return numerical_method_.flux_method;
}

const HyperbolicSplitTimeIntegrator&
HyperbolicSplitIntegratorContext::GetHyperbolicSplitTimeIntegrator() const
    noexcept {
  return numerical_method_.time_integrator;
}

const Reconstruction&
HyperbolicSplitIntegratorContext::GetReconstruction() const noexcept {
  return numerical_method_.reconstruction;
}

const std::shared_ptr<GriddingAlgorithm>&
HyperbolicSplitIntegratorContext::GetGriddingAlgorithm() const noexcept {
  return gridding_;
}

PatchHierarchy& HyperbolicSplitIntegratorContext::GetPatchHierarchy() noexcept {
  return gridding_->GetPatchHierarchy();
}

const PatchHierarchy&
HyperbolicSplitIntegratorContext::GetPatchHierarchy() const noexcept {
  return gridding_->GetPatchHierarchy();
}

MPI_Comm HyperbolicSplitIntegratorContext::GetMpiCommunicator() const noexcept {
  return ::amrex::ParallelContext::CommunicatorAll();
}

::amrex::MultiFab& HyperbolicSplitIntegratorContext::GetData(int level) {
  return GetPatchHierarchy().GetPatchLevel(level).data;
}

::amrex::MultiFab& HyperbolicSplitIntegratorContext::GetScratch(int level,
                                                                Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].scratch[d];
}

::amrex::MultiFab& HyperbolicSplitIntegratorContext::GetFluxes(int level,
                                                               Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].fluxes[d];
}

const ::amrex::Geometry&
HyperbolicSplitIntegratorContext::GetGeometry(int level) const {
  return GetPatchHierarchy().GetGeometry(level);
}

Duration HyperbolicSplitIntegratorContext::GetTimePoint(int level,
                                                        Direction dir) const {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].time_point[d];
}

std::ptrdiff_t
HyperbolicSplitIntegratorContext::GetCycles(int level, Direction dir) const {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].cycles[d];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                    Observers

bool HyperbolicSplitIntegratorContext::LevelExists(int level) const noexcept {
  return 0 <= level && level < GetPatchHierarchy().GetNumberOfLevels();
}

namespace {
constexpr std::ptrdiff_t ipow(int base, int exponent) {
  std::ptrdiff_t prod{1};
  while (exponent > 0) {
    prod *= base;
    exponent -= 1;
  }
  return prod;
}
} // namespace

int HyperbolicSplitIntegratorContext::GetRatioToCoarserLevel(
    int level, Direction dir) const noexcept {
  return GetPatchHierarchy().GetRatioToCoarserLevel(level, dir);
}

::amrex::IntVect
HyperbolicSplitIntegratorContext::GetRatioToCoarserLevel(int level) const
    noexcept {
  return GetPatchHierarchy().GetRatioToCoarserLevel(level);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                    Modifiers

void HyperbolicSplitIntegratorContext::ResetHierarchyConfiguration(
    std::shared_ptr<GriddingAlgorithm> gridding) {
  gridding_ = std::move(gridding);
  ResetHierarchyConfiguration();
}

void HyperbolicSplitIntegratorContext::ResetHierarchyConfiguration(
    int first_level) {
  const int n_cons_components =
      GetPatchHierarchy().GetDataDescription().n_cons_components;
  const int new_n_levels = GetPatchHierarchy().GetNumberOfLevels();
  data_.resize(static_cast<std::size_t>(new_n_levels));
  for (int level = first_level; level < new_n_levels; ++level) {
    LevelData& data = data_[static_cast<std::size_t>(level)];
    const ::amrex::BoxArray& ba =
        GetPatchHierarchy().GetPatchLevel(level).box_array;
    const ::amrex::DistributionMapping& dm =
        GetPatchHierarchy().GetPatchLevel(level).distribution_mapping;
    const int n_comp =
        GetPatchHierarchy().GetDataDescription().n_state_components;
    // We need to allocate data for each dimensions because of
    // amrex::FluxBoundary::OverwriteFlux.
    for (std::size_t d = 0; d < static_cast<std::size_t>(AMREX_SPACEDIM); ++d) {
      const ::amrex::IntVect unit =
          ::amrex::IntVect::TheDimensionVector(int(d));
      data.scratch[d].define(ba, dm, n_comp, ghost_cell_width_ * unit);
      data.fluxes[d].define(::amrex::convert(ba, unit), dm, n_cons_components,
                            1 * unit);
    }
    if (level > 0) {
      const ::amrex::IntVect ref_ratio =
          GetPatchHierarchy().GetRatioToCoarserLevel(level);
      data.coarse_fine.clear();
      data.coarse_fine.define(ba, dm, ref_ratio, level, n_cons_components);
    }
  }
}

void HyperbolicSplitIntegratorContext::SetTimePoint(Duration dt, int level,
                                                    Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  data_[l].time_point[d] = dt;
}

void HyperbolicSplitIntegratorContext::SetCycles(std::ptrdiff_t cycles,
                                                 int level, Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  data_[l].cycles[d] = cycles;
}

void HyperbolicSplitIntegratorContext::FillGhostLayerTwoLevels(int fine,
                                                               int coarse,
                                                               Direction dir) {
  FUB_ASSERT(coarse >= 0 && fine > coarse);
  ::amrex::MultiFab& scratch = GetScratch(fine, dir);
  ::amrex::Vector<::amrex::BCRec> bcr(
      static_cast<std::size_t>(scratch.nComp()));
  const int nc = scratch.nComp();
  const ::amrex::Vector<::amrex::MultiFab*> cmf{&GetData(coarse)};
  const ::amrex::Vector<::amrex::MultiFab*> fmf{&GetData(fine)};
  const ::amrex::Vector<double> ct{GetTimePoint(coarse, dir).count()};
  const ::amrex::Vector<double> ft{GetTimePoint(fine, dir).count()};
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  const ::amrex::Geometry& fgeom = GetGeometry(fine);
  const ::amrex::IntVect ratio = 2 * ::amrex::IntVect::TheUnitVector();
  ::amrex::Interpolater* mapper = &::amrex::pc_interp;
  BoundaryCondition& fine_condition = GetBoundaryCondition(fine);
  BoundaryCondition& coarse_condition = GetBoundaryCondition(coarse);
  ::amrex::FillPatchTwoLevels(scratch, ft[0], cmf, ct, fmf, ft, 0, 0, nc, cgeom,
                              fgeom, coarse_condition, 0, fine_condition, 0,
                              ratio, mapper, bcr, 0);
}

void HyperbolicSplitIntegratorContext::FillGhostLayerSingleLevel(
    int level, Direction dir) {
  ::amrex::MultiFab& scratch = GetScratch(level, dir);
  ::amrex::Vector<::amrex::BCRec> bcr(
      static_cast<std::size_t>(scratch.nComp()));
  const int nc = scratch.nComp();
  const ::amrex::Vector<::amrex::MultiFab*> smf{&GetData(level)};
  const ::amrex::Vector<double> stime{GetTimePoint(level, dir).count()};
  const ::amrex::Geometry& geom = GetGeometry(level);
  BoundaryCondition& condition = GetBoundaryCondition(level);
  ::amrex::FillPatchSingleLevel(scratch, stime[0], smf, stime, 0, 0, nc, geom,
                                condition, 0);
}

void HyperbolicSplitIntegratorContext::CoarsenConservatively(int fine_level,
                                                             int coarse_level,
                                                             Direction dir) {
  const int first =
      GetPatchHierarchy().GetDataDescription().first_cons_component;
  const int size = GetPatchHierarchy().GetDataDescription().n_cons_components;
  const ::amrex::MultiFab& fine = GetScratch(fine_level, dir);
  ::amrex::MultiFab& coarse = GetScratch(coarse_level, dir);
  FUB_ASSERT(!fine.contains_nan());
  FUB_ASSERT(!coarse.contains_nan());
  ::amrex::average_down(fine, coarse, GetGeometry(fine_level),
                        GetGeometry(coarse_level), first, size,
                        GetRatioToCoarserLevel(fine_level));
}

void HyperbolicSplitIntegratorContext::AccumulateCoarseFineFluxes(
    int level, Duration, Direction dir) {
  if (level > 0) {
    const ::amrex::MultiFab& fluxes = GetFluxes(level, dir);
    const int dim = GetPatchHierarchy().GetDataDescription().dimension;
    const double scale =
        1.0 /
        static_cast<double>(ipow(GetRatioToCoarserLevel(level, dir), dim));
    data_[static_cast<std::size_t>(level)].coarse_fine.FineAdd(
        fluxes, int(dir), 0, 0, fluxes.nComp(), scale);
  }
}

void HyperbolicSplitIntegratorContext::ResetCoarseFineFluxes(int fine,
                                                             int coarse,
                                                             Direction dir) {
  std::size_t sfine = static_cast<std::size_t>(fine);
  data_[sfine].coarse_fine.ClearInternalBorders(GetGeometry(coarse));
  const ::amrex::MultiFab& flux = GetFluxes(coarse, dir);
  const ::amrex::BoxArray& boxes = flux.boxArray();
  const ::amrex::DistributionMapping& distribution = flux.DistributionMap();
  const int ncomp = flux.nComp();
  ::amrex::MultiFab zero(boxes, distribution, ncomp, 0);
  zero.setVal(0.0);
  data_[sfine].coarse_fine.CrseInit(zero, int(dir), 0, 0, ncomp);
}

void HyperbolicSplitIntegratorContext::ApplyFluxCorrection(int fine, int coarse,
                                                           Duration,
                                                           Direction dir) {
  const std::size_t sfine = static_cast<std::size_t>(fine);
  const int ncomp = GetPatchHierarchy().GetDataDescription().n_cons_components;
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  std::array<::amrex::MultiFab*, AMREX_SPACEDIM> crse_fluxes{AMREX_D_DECL(
      &GetFluxes(coarse, Direction::X), &GetFluxes(coarse, Direction::Y),
      &GetFluxes(coarse, Direction::Z))};
  const auto size = crse_fluxes[std::size_t(dir)]->size();
  if (size > 0) {
    data_[sfine].coarse_fine.OverwriteFlux(crse_fluxes, 1.0, 0, 0, ncomp,
                                           cgeom);
  }
}

Duration HyperbolicSplitIntegratorContext::ComputeStableDt(int level,
                                                         Direction dir) {
  return numerical_method_.flux_method.ComputeStableDt(*this, level, dir);
}

void HyperbolicSplitIntegratorContext::ComputeNumericFluxes(int level,
                                                            Duration dt,
                                                            Direction dir) {
  numerical_method_.flux_method.ComputeNumericFluxes(*this, level, dt, dir);
}

void HyperbolicSplitIntegratorContext::CompleteFromCons(int level, Duration,
                                                        Direction dir) {
  const std::size_t l = static_cast<std::size_t>(level);
  const std::size_t d = static_cast<std::size_t>(dir);
  ::amrex::MultiFab& data = GetPatchHierarchy().GetPatchLevel(level).data;
  const ::amrex::MultiFab& scratch = data_[l].scratch[d];
  numerical_method_.reconstruction.CompleteFromCons(data, scratch);
}

void HyperbolicSplitIntegratorContext::UpdateConservatively(int level,
                                                            Duration dt,
                                                            Direction dir) {
  const std::size_t l = static_cast<std::size_t>(level);
  const std::size_t d = static_cast<std::size_t>(dir);
  ::amrex::MultiFab& scratch = data_[l].scratch[d];
  const ::amrex::MultiFab& fluxes = data_[l].fluxes[d];
  const ::amrex::Geometry& geom = GetGeometry(level);
  numerical_method_.time_integrator.UpdateConservatively(scratch, scratch,
                                                         fluxes, geom, dt, dir);
}

void HyperbolicSplitIntegratorContext::PreAdvanceLevel(int level_num,
                                                       Direction dir, Duration,
                                                       int subcycle) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level_num);
  if (subcycle == 0 && level_num > 0 &&
      data_[l].regrid_time_point[d] != data_[l].time_point[d]) {
    gridding_->RegridAllFinerlevels(level_num - 1);
    for (std::size_t lvl = l; lvl < data_.size(); ++lvl) {
      data_[lvl].regrid_time_point[d] = data_[lvl].time_point[d];
    }
    ResetHierarchyConfiguration(level_num);
    ResetCoarseFineFluxes(level_num, level_num - 1, dir);
  }
}

Result<void, TimeStepTooLarge>
HyperbolicSplitIntegratorContext::PostAdvanceLevel(int level_num, Direction dir,
                                                   Duration dt, int) {
  SetCycles(GetCycles(level_num, dir) + 1, level_num, dir);
  SetTimePoint(GetTimePoint(level_num, dir) + dt, level_num, dir);
  GetPatchHierarchy().GetPatchLevel(level_num).time_point =
      GetTimePoint(level_num, dir);
  GetPatchHierarchy().GetPatchLevel(level_num).cycles =
      GetCycles(level_num, dir);
  return boost::outcome_v2::success();
}

} // namespace fub::amrex
