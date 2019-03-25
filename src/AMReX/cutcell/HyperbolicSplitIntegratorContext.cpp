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

#include "fub/grid/AMReX/cutcell/HyperbolicSplitIntegratorContext.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"
#include "fub/grid/AMReX/utility.hpp"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxReg_C.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>

#include <fmt/format.h>

namespace fub {
namespace amrex {
namespace cutcell {

namespace {
PatchDataView<double, AMREX_SPACEDIM + 1>
GetPatchDataView(::amrex::MultiFab& multi_fab, PatchHandle patch) {
  ::amrex::FArrayBox& fab = multi_fab[*patch.iterator];
  return MakePatchDataView(fab);
}
} // namespace

bool HyperbolicSplitIntegratorContext::LevelExists(int level) const noexcept {
  return level < GetPatchHierarchy()->GetNumberOfLevels();
}

int HyperbolicSplitIntegratorContext::GetRatioToCoarserLevel(int level) const
    noexcept {
  if (level) {
    return 2;
  }
  return 1;
}

const std::shared_ptr<PatchHierarchy>&
HyperbolicSplitIntegratorContext::GetPatchHierarchy() const noexcept {
  return gridding_->GetPatchHierarchy();
}

BoundaryCondition
HyperbolicSplitIntegratorContext::GetBoundaryCondition(int level) const {
  const GriddingAlgorithm::BoundaryCondition& fn =
      gridding_->GetBoundaryCondition();
  BoundaryCondition bc(fn, GetGeometry(level), level);
  return bc;
}

int HyperbolicSplitIntegratorContext::GetGhostCellWidth(PatchHandle,
                                                        Direction) {
  return ghost_cell_width_;
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetData(PatchHandle patch) {
  return GetPatchDataView(GetPatchHierarchy()->GetPatchLevel(patch.level).data,
                          patch);
}

::amrex::MultiFab& HyperbolicSplitIntegratorContext::GetData(int level) {
  return GetPatchHierarchy()->GetPatchLevel(level).data;
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetScratch(PatchHandle patch, Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(patch.level);
  return GetPatchDataView(data_[level_num].scratch[d], patch);
}

::amrex::MultiFab& HyperbolicSplitIntegratorContext::GetScratch(int level,
                                                                Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].scratch[d];
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetFluxes(PatchHandle patch, Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(patch.level);
  return GetPatchDataView(data_[level_num].fluxes[d], patch);
}

::amrex::MultiFab& HyperbolicSplitIntegratorContext::GetFluxes(int level,
                                                               Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].fluxes[d];
}

::amrex::MultiFab&
HyperbolicSplitIntegratorContext::GetBoundaryFluxes(int level, Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].boundary_fluxes[d];
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetBoundaryFluxes(PatchHandle patch,
                                                    Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(patch.level);
  return GetPatchDataView(data_[level_num].boundary_fluxes[d], patch);
}

::amrex::MultiFab&
HyperbolicSplitIntegratorContext::GetReferenceStates(int level, Direction) {
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].reference_states;
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetReferenceStates(PatchHandle patch,
                                                     Direction /* dir */) {
  const std::size_t level_num = static_cast<std::size_t>(patch.level);
  return GetPatchDataView(data_[level_num].reference_states, patch);
}

::amrex::MultiFab&
HyperbolicSplitIntegratorContext::GetStabilizedFluxes(int level,
                                                      Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].stabilized_fluxes[d];
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetStabilizedFluxes(PatchHandle patch,
                                                      Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(patch.level);
  return GetPatchDataView(data_[level_num].stabilized_fluxes[d], patch);
}

CutCellData<AMREX_SPACEDIM>
HyperbolicSplitIntegratorContext::GetCutCellData(PatchHandle patch,
                                                 Direction dir) {
  return GetPatchHierarchy()->GetCutCellData(patch, dir);
}

const ::amrex::Geometry&
HyperbolicSplitIntegratorContext::GetGeometry(int level) const {
  return GetPatchHierarchy()->GetGeometry(level);
}

Duration HyperbolicSplitIntegratorContext::GetTimePoint(int level,
                                                        Direction dir) const {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].time_point[d];
}

std::ptrdiff_t
HyperbolicSplitIntegratorContext::GetCycles(int level, Direction dir) const {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].cycles[d];
}

void HyperbolicSplitIntegratorContext::SetTimePoint(Duration dt, int level,
                                                    Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  data_[level_num].time_point[d] = dt;
}

void HyperbolicSplitIntegratorContext::SetCycles(std::ptrdiff_t cycles,
                                                 int level, Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  data_[level_num].cycles[d] = cycles;
}

double HyperbolicSplitIntegratorContext::GetDx(PatchHandle patch,
                                               Direction dir) const {
  const int d = static_cast<int>(dir);
  return GetPatchHierarchy()->GetGeometry(patch.level).CellSize(d);
}

HyperbolicSplitIntegratorContext::HyperbolicSplitIntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, int gcw)
    : ghost_cell_width_{gcw + 1}, gridding_{std::move(gridding)},
  data_(static_cast<std::size_t>(GetPatchHierarchy()->GetMaxNumberOfLevels())) {
  ResetHierarchyConfiguration();
}

void HyperbolicSplitIntegratorContext::ResetHierarchyConfiguration(
    int first_level) {
  const int n_cons_components =
      GetPatchHierarchy()->GetDataDescription().n_cons_components;
  const int n_levels = GetPatchHierarchy()->GetMaxNumberOfLevels();
  for (int level = first_level; level < n_levels; ++level) {
    LevelData& data = data_[static_cast<std::size_t>(level)];
    const ::amrex::MultiFab& base =
        GetPatchHierarchy()->GetPatchLevel(level).data;
    const ::amrex::BoxArray& ba = base.boxArray();
    const ::amrex::DistributionMapping& dm = base.DistributionMap();
    const int n_comp = base.nComp();
    data.reference_states.define(ba, dm, n_comp, 0);
    for (std::size_t d = 0; d < static_cast<std::size_t>(AMREX_SPACEDIM); ++d) {
      const ::amrex::IntVect unit = ::amrex::IntVect::TheDimensionVector(static_cast<int>(d));
      data.scratch[d].define(ba, dm, n_comp, ghost_cell_width_ * unit);
      data.fluxes[d].define(::amrex::convert(ba, unit), dm, n_cons_components,
                            unit);
      data.boundary_fluxes[d].define(ba, dm, n_cons_components,
                                     ghost_cell_width_ * unit);
      data.stabilized_fluxes[d].define(::amrex::convert(ba, unit), dm,
                                       n_cons_components, unit);
    }
    if (level > 0) {
      const ::amrex::IntVect ref_ratio = 2 * ::amrex::IntVect::TheUnitVector();
      data.coarse_fine.clear();
      data.coarse_fine.define(ba, dm, ref_ratio, level, n_cons_components);
    }
  }
}

void HyperbolicSplitIntegratorContext::FillGhostLayerTwoLevels(int fine,
                                                               int coarse,
                                                               Direction dir) {
  FUB_ASSERT(coarse >= 0 && fine > coarse);
  ::amrex::MultiFab& scratch = GetScratch(fine, dir);
  const int nc = scratch.nComp();
  FUB_ASSERT(nc > 0);
  const std::size_t ncs = static_cast<std::size_t>(nc);
  ::amrex::Vector<::amrex::BCRec> bcr(ncs);
  const ::amrex::Vector<::amrex::MultiFab*> cmf{&GetData(coarse)};
  const ::amrex::Vector<::amrex::MultiFab*> fmf{&GetData(fine)};
  const ::amrex::Vector<double> ct{GetTimePoint(coarse, dir).count()};
  const ::amrex::Vector<double> ft{GetTimePoint(fine, dir).count()};
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  const ::amrex::Geometry& fgeom = GetGeometry(fine);
  const ::amrex::IntVect ratio = 2 * ::amrex::IntVect::TheUnitVector();
  ::amrex::Interpolater* mapper = &::amrex::pc_interp;
  BoundaryCondition fine_condition = GetBoundaryCondition(fine);
  BoundaryCondition coarse_condition = GetBoundaryCondition(coarse);
  ::amrex::FillPatchTwoLevels(scratch, ft[0], cmf, ct, fmf, ft, 0, 0, nc, cgeom,
                              fgeom, coarse_condition, 0, fine_condition, 0,
                              ratio, mapper, bcr, 0);
}

void HyperbolicSplitIntegratorContext::FillGhostLayerSingleLevel(
    int level, Direction dir) {
  ::amrex::MultiFab& scratch = GetScratch(level, dir);
  const int nc = scratch.nComp();
  FUB_ASSERT(nc > 0);
  const std::size_t ncs = static_cast<std::size_t>(nc);
  ::amrex::Vector<::amrex::BCRec> bcr(ncs);
  const ::amrex::Vector<::amrex::MultiFab*> smf{&GetData(level)};
  const ::amrex::Vector<double> stime{GetTimePoint(level, dir).count()};
  const ::amrex::Geometry& geom = GetGeometry(level);
  BoundaryCondition condition = GetBoundaryCondition(level);
  ::amrex::FillPatchSingleLevel(scratch, stime[0], smf, stime, 0, 0, nc, geom,
                                condition, 0);
}

void HyperbolicSplitIntegratorContext::CoarsenConservatively(int fine_level,
                                                             int coarse_level,
                                                             Direction dir) {
  const int first =
      GetPatchHierarchy()->GetDataDescription().first_cons_component;
  const int size = GetPatchHierarchy()->GetDataDescription().n_cons_components;
  ::amrex::average_down(GetScratch(fine_level, dir),
                        GetScratch(coarse_level, dir), GetGeometry(fine_level),
                        GetGeometry(coarse_level), first, size, 2);
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

void HyperbolicSplitIntegratorContext::AccumulateCoarseFineFluxes(int level,
                                                                  Direction dir,
                                                                  Duration) {
  if (level > 0) {
    const std::size_t level_num = static_cast<std::size_t>(level);
    const ::amrex::MultiFab& fluxes = GetFluxes(level, dir);
    const double scale =
        1.0 / ipow(GetRatioToCoarserLevel(level), AMREX_SPACEDIM);
    data_[level_num].coarse_fine.FineAdd(fluxes, static_cast<int>(dir), 0, 0,
                                     fluxes.nComp(), scale);
  }
}

void HyperbolicSplitIntegratorContext::ResetCoarseFineFluxes(int fine,
                                                             int coarse,
                                                             Direction dir) {
  FUB_ASSERT(fine > 0);
  const std::size_t f = static_cast<std::size_t>(fine);
  data_[f].coarse_fine.ClearInternalBorders(GetGeometry(coarse));
  const ::amrex::MultiFab& flux = GetFluxes(coarse, dir);
  const ::amrex::BoxArray& boxes = flux.boxArray();
  const ::amrex::DistributionMapping& distribution = flux.DistributionMap();
  const int ncomp = flux.nComp();
  ::amrex::MultiFab zero(boxes, distribution, ncomp, 0);
  zero.setVal(0.0);
  data_[f].coarse_fine.CrseInit(zero, static_cast<int>(dir), 0, 0, ncomp);
}

void HyperbolicSplitIntegratorContext::ApplyFluxCorrection(int fine, int coarse,
                                                           Duration,
                                                           Direction) {
  FUB_ASSERT(fine > 0);
  const std::size_t f = static_cast<std::size_t>(fine);
  const int ncomp = GetPatchHierarchy()->GetDataDescription().n_cons_components;
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  std::array<::amrex::MultiFab*, AMREX_SPACEDIM> crse_fluxes{AMREX_D_DECL(
      &GetFluxes(coarse, Direction::X), &GetFluxes(coarse, Direction::Y),
      &GetFluxes(coarse, Direction::Z))};
  data_[f].coarse_fine.OverwriteFlux(crse_fluxes, 1.0, 0, 0, ncomp, cgeom);
}

CartesianCoordinates HyperbolicSplitIntegratorContext::GetCartesianCoordinates(
    PatchHandle patch) const {
  const ::amrex::Geometry& geom = GetGeometry(patch.level);
  const ::amrex::Box& box = patch.iterator->tilebox();
  return ::fub::amrex::GetCartesianCoordinates(geom, box);
}

MPI_Comm HyperbolicSplitIntegratorContext::GetMpiCommunicator() const noexcept {
  return ::amrex::ParallelContext::CommunicatorAll();
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

void HyperbolicSplitIntegratorContext::PostAdvanceLevel(int level_num,
                                                        Direction dir,
                                                        Duration dt, int) {
  SetCycles(GetCycles(level_num, dir) + 1, level_num, dir);
  SetTimePoint(GetTimePoint(level_num, dir) + dt, level_num, dir);
  GetPatchHierarchy()->GetPatchLevel(level_num).time_point =
      GetTimePoint(level_num, dir);
  GetPatchHierarchy()->GetPatchLevel(level_num).cycles =
      GetCycles(level_num, dir);
}

::amrex::FabType HyperbolicSplitIntegratorContext::GetCutCellPatchType(
    PatchHandle handle, int gcw) const {
  const ::amrex::EBFArrayBoxFactory& eb_factory =
      GetPatchHierarchy()->GetEmbeddedBoundary(handle.level);
  const ::amrex::FabArray<::amrex::EBCellFlagFab>& flags =
      eb_factory.getMultiEBCellFlagFab();
  const ::amrex::EBCellFlagFab& ffab = flags[*handle.iterator];
  ::amrex::Box box = handle.iterator->tilebox();
  box.grow(gcw);
  return ffab.getType(box);
}

} // namespace cutcell
} // namespace amrex
} // namespace fub
