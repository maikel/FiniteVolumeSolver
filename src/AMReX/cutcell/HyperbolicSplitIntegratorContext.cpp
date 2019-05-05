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
#include "fub/grid/AMReX/cutcell/BoundaryCondition.hpp"
#include "fub/grid/AMReX/cutcell/IndexSpace.hpp"
#include "fub/grid/AMReX/utility.hpp"

#include <AMReX_EBMultiFabUtil.H>
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

PatchDataView<double, AMREX_SPACEDIM + 1>
GetPatchDataView(::amrex::MultiCutFab& multi_fab, PatchHandle patch) {
  ::amrex::FArrayBox& fab = multi_fab[*patch.iterator];
  return MakePatchDataView(fab);
}
} // namespace

bool HyperbolicSplitIntegratorContext::LevelExists(int level) const noexcept {
  return level < GetPatchHierarchy().GetNumberOfLevels();
}

int HyperbolicSplitIntegratorContext::GetRatioToCoarserLevel(
    int level, Direction dir) const noexcept {
  return GetPatchHierarchy().GetRatioToCoarserLevel(level, dir);
}

::amrex::IntVect
HyperbolicSplitIntegratorContext::GetRatioToCoarserLevel(int level) const
    noexcept {
  return GetPatchHierarchy().GetRatioToCoarserLevel(level);
}

PatchHierarchy& HyperbolicSplitIntegratorContext::GetPatchHierarchy() noexcept {
  return gridding_->GetPatchHierarchy();
}

const PatchHierarchy&
HyperbolicSplitIntegratorContext::GetPatchHierarchy() const noexcept {
  return gridding_->GetPatchHierarchy();
}

BoundaryCondition
HyperbolicSplitIntegratorContext::GetBoundaryCondition(int level) const {
  const GriddingAlgorithm::BoundaryCondition& fn =
      gridding_->GetBoundaryCondition();
  BoundaryCondition bc(fn, GetGeometry(level), level, GetPatchHierarchy());
  return bc;
}

int HyperbolicSplitIntegratorContext::GetGhostCellWidth(PatchHandle,
                                                        Direction) {
  return ghost_cell_width_;
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetData(PatchHandle patch) {
  return GetPatchDataView(GetData(patch.level), patch);
}

::amrex::MultiFab& HyperbolicSplitIntegratorContext::GetData(int level) {
  return GetPatchHierarchy().GetPatchLevel(level).data;
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

::amrex::MultiCutFab&
HyperbolicSplitIntegratorContext::GetBoundaryFluxes(int level, Direction) {
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].boundary_fluxes;
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetBoundaryFluxes(PatchHandle patch,
                                                    Direction) {
  const std::size_t level_num = static_cast<std::size_t>(patch.level);
  return GetPatchDataView(data_[level_num].boundary_fluxes, patch);
}

::amrex::MultiFab&
HyperbolicSplitIntegratorContext::GetReferenceStates(int level, Direction) {
  const std::size_t level_num = static_cast<std::size_t>(level);
  return *data_[level_num].reference_states;
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetReferenceStates(PatchHandle patch) {
  const std::size_t level_num = static_cast<std::size_t>(patch.level);
  return GetPatchDataView(*data_[level_num].reference_states, patch);
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

::amrex::MultiFab&
HyperbolicSplitIntegratorContext::GetShieldedLeftFluxes(int level,
                                                        Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].shielded_left_fluxes[d];
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetShieldedLeftFluxes(PatchHandle patch,
                                                        Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(patch.level);
  return GetPatchDataView(data_[level_num].shielded_left_fluxes[d], patch);
}

::amrex::MultiFab&
HyperbolicSplitIntegratorContext::GetShieldedRightFluxes(int level,
                                                         Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].shielded_right_fluxes[d];
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetShieldedRightFluxes(PatchHandle patch,
                                                         Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(patch.level);
  return GetPatchDataView(data_[level_num].shielded_right_fluxes[d], patch);
}

::amrex::MultiFab&
HyperbolicSplitIntegratorContext::GetDoublyShieldedFluxes(int level,
                                                          Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(level);
  return data_[level_num].doubly_shielded_fluxes[d];
}

PatchDataView<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetDoublyShieldedFluxes(PatchHandle patch,
                                                          Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t level_num = static_cast<std::size_t>(patch.level);
  return GetPatchDataView(data_[level_num].doubly_shielded_fluxes[d], patch);
}

CutCellData<AMREX_SPACEDIM>
HyperbolicSplitIntegratorContext::GetCutCellData(PatchHandle patch,
                                                 Direction dir) {
  return GetPatchHierarchy().GetCutCellData(patch, dir);
}

const ::amrex::Geometry&
HyperbolicSplitIntegratorContext::GetGeometry(int level) const {
  return GetPatchHierarchy().GetGeometry(level);
}

const std::shared_ptr<GriddingAlgorithm>&
HyperbolicSplitIntegratorContext::GetGriddingAlgorithm() const noexcept {
  return gridding_;
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
  return GetPatchHierarchy().GetGeometry(patch.level).CellSize(d);
}

HyperbolicSplitIntegratorContext::HyperbolicSplitIntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, int gcw)
    : ghost_cell_width_{gcw + 1}, gridding_{std::move(gridding)},
      data_(static_cast<std::size_t>(
          GetPatchHierarchy().GetMaxNumberOfLevels())) {
  ResetHierarchyConfiguration();
  const int nlevels = GetPatchHierarchy().GetNumberOfLevels();
  for (int level = 0; level < nlevels; ++level) {
    const Duration time_point = GetPatchHierarchy().GetTimePoint(level);
    const std::ptrdiff_t cycles = GetPatchHierarchy().GetCycles(level);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      SetTimePoint(time_point, level, Direction(d));
      SetCycles(cycles, level, Direction(d));
    }
  }
}

HyperbolicSplitIntegratorContext::HyperbolicSplitIntegratorContext(
    const HyperbolicSplitIntegratorContext& other)
    : ghost_cell_width_{other.ghost_cell_width_}, gridding_{other.gridding_},
      data_(static_cast<std::size_t>(
          GetPatchHierarchy().GetMaxNumberOfLevels())) {
  // Allocate data arrays
  ResetHierarchyConfiguration();
  // Copy relevant data
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

void HyperbolicSplitIntegratorContext::ResetHierarchyConfiguration(
    std::shared_ptr<GriddingAlgorithm> gridding) {
  gridding_ = std::move(gridding);
  ResetHierarchyConfiguration();
}

void HyperbolicSplitIntegratorContext::ResetHierarchyConfiguration(
    int first_level) {
  const int n_components =
      GetPatchHierarchy().GetDataDescription().n_state_components;
  const int n_cons_components =
      GetPatchHierarchy().GetDataDescription().n_cons_components;
  const int n_levels = GetPatchHierarchy().GetMaxNumberOfLevels();
  for (int level = first_level; level < n_levels; ++level) {
    LevelData& data = data_[static_cast<std::size_t>(level)];

    const ::amrex::BoxArray& ba =
        GetPatchHierarchy().GetPatchLevel(level).box_array;
    const ::amrex::DistributionMapping& dm =
        GetPatchHierarchy().GetPatchLevel(level).distribution_mapping;
    const std::shared_ptr<::amrex::EBFArrayBoxFactory>& ebf =
        GetPatchHierarchy().GetEmbeddedBoundary(level);

    { // Redistribute reference states
      std::unique_ptr<::amrex::MultiFab> refs =
          std::make_unique<::amrex::MultiFab>(
              ba, dm, n_components, ghost_cell_width_, ::amrex::MFInfo(), *ebf);
      if (data.reference_states) {
        refs->ParallelCopy(*data.reference_states, 0, 0, n_components,
                           ghost_cell_width_, ghost_cell_width_);
      }
      data.reference_states = std::move(refs);
    }

    data.boundary_fluxes.define(ba, dm, n_cons_components, ghost_cell_width_,
                                ebf->getMultiEBCellFlagFab());

    for (std::size_t d = 0; d < static_cast<std::size_t>(AMREX_SPACEDIM); ++d) {
      const ::amrex::IntVect unit =
          ::amrex::IntVect::TheDimensionVector(static_cast<int>(d));
      data.scratch[d].define(ba, dm, n_components, ghost_cell_width_ * unit);
      data.fluxes[d].define(::amrex::convert(ba, unit), dm, n_cons_components,
                            unit);
      data.stabilized_fluxes[d].define(::amrex::convert(ba, unit), dm,
                                       n_cons_components, unit);
      data.shielded_left_fluxes[d].define(::amrex::convert(ba, unit), dm,
                                          n_cons_components, unit);
      data.shielded_right_fluxes[d].define(::amrex::convert(ba, unit), dm,
                                           n_cons_components, unit);
      data.doubly_shielded_fluxes[d].define(::amrex::convert(ba, unit), dm,
                                            n_cons_components, unit);
    }

    data.eb_factory = ebf;
    if (level > 0) {
      const ::amrex::IntVect ref_ratio = GetRatioToCoarserLevel(level);
      data.coarse_fine.clear();
      data.coarse_fine.define(ba, dm, ref_ratio, level, n_cons_components);
    }
  }
}

namespace {
void FillGhostLayerTwoLevels_(::amrex::MultiFab& dest,
                              HyperbolicSplitIntegratorContext& context,
                              int fine, int coarse, Direction dir) {
  FUB_ASSERT(coarse >= 0 && fine > coarse);
  const int nc = dest.nComp();
  FUB_ASSERT(nc > 0);
  const std::size_t ncs = static_cast<std::size_t>(nc);
  ::amrex::Vector<::amrex::BCRec> bcr(ncs);
  const ::amrex::Vector<::amrex::MultiFab*> cmf{&context.GetData(coarse)};
  const ::amrex::Vector<::amrex::MultiFab*> fmf{&context.GetData(fine)};
  const ::amrex::Vector<double> ct{context.GetTimePoint(coarse, dir).count()};
  const ::amrex::Vector<double> ft{context.GetTimePoint(fine, dir).count()};
  const ::amrex::Geometry& cgeom = context.GetGeometry(coarse);
  const ::amrex::Geometry& fgeom = context.GetGeometry(fine);
  const ::amrex::IntVect ratio = context.GetRatioToCoarserLevel(fine);
  ::amrex::Interpolater* mapper = &::amrex::pc_interp;
  BoundaryCondition fine_condition = context.GetBoundaryCondition(fine);
  BoundaryCondition coarse_condition = context.GetBoundaryCondition(coarse);
  const ::amrex::EB2::IndexSpace* index_space = context.GetPatchHierarchy()
                                                    .GetPatchLevel(fine)
                                                    .factory->getEBIndexSpace();
  FUB_ASSERT(index_space != nullptr);
  ::amrex::FillPatchTwoLevels(
      dest, ft[0], *index_space, cmf, ct, fmf, ft, 0, 0, nc, cgeom, fgeom,
      coarse_condition, 0, fine_condition, 0, ratio, mapper, bcr, 0,
      ::amrex::NullInterpHook(), ::amrex::NullInterpHook());
}

void FillGhostLayerSingleLevel_(::amrex::MultiFab& dest,
                                HyperbolicSplitIntegratorContext& context,
                                int level, Direction dir) {
  const int nc = dest.nComp();
  FUB_ASSERT(nc > 0);
  const std::size_t ncs = static_cast<std::size_t>(nc);
  ::amrex::Vector<::amrex::BCRec> bcr(ncs);
  const ::amrex::Vector<::amrex::MultiFab*> smf{&context.GetData(level)};
  const ::amrex::Vector<double> stime{context.GetTimePoint(level, dir).count()};
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  BoundaryCondition condition = context.GetBoundaryCondition(level);
  ::amrex::FillPatchSingleLevel(dest, stime[0], smf, stime, 0, 0, nc, geom,
                                condition, 0);
}
} // namespace

void HyperbolicSplitIntegratorContext::FillGhostLayer(::amrex::MultiFab& dest,
                                                      int level) {
  if (level > 0) {
    FillGhostLayerTwoLevels_(dest, *this, level, level - 1, Direction::X);
  } else {
    FillGhostLayerSingleLevel_(dest, *this, level, Direction::X);
  }
}

void HyperbolicSplitIntegratorContext::FillGhostLayerTwoLevels(int fine,
                                                               int coarse,
                                                               Direction dir) {
  ::amrex::MultiFab& scratch = GetScratch(fine, dir);
  FillGhostLayerTwoLevels_(scratch, *this, fine, coarse, dir);
}

void HyperbolicSplitIntegratorContext::FillGhostLayerSingleLevel(
    int level, Direction dir) {
  ::amrex::MultiFab& scratch = GetScratch(level, dir);
  FillGhostLayerSingleLevel_(scratch, *this, level, dir);
}

void HyperbolicSplitIntegratorContext::CoarsenConservatively(int fine_level,
                                                             int coarse_level,
                                                             Direction dir) {
  const int first =
      GetPatchHierarchy().GetDataDescription().first_cons_component;
  const int size = GetPatchHierarchy().GetDataDescription().n_cons_components;
  ::amrex::EB_average_down(GetScratch(fine_level, dir),
                           GetScratch(coarse_level, dir), first, size, 2);
}

void HyperbolicSplitIntegratorContext::AccumulateCoarseFineFluxes(int level,
                                                                  Direction dir,
                                                                  Duration) {
  if (level > 0) {
    const std::size_t level_num = static_cast<std::size_t>(level);
    const ::amrex::MultiFab& fluxes = GetFluxes(level, dir);
    const ::amrex::IntVect ratio = GetRatioToCoarserLevel(level);
    const double scale =
        1.0 / static_cast<double>(AMREX_D_TERM(ratio[0], *ratio[1], *ratio[2]));
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
  const int ncomp = GetPatchHierarchy().GetDataDescription().n_cons_components;
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
}

::amrex::FabType
HyperbolicSplitIntegratorContext::GetCutCellPatchType(PatchHandle handle,
                                                      int gcw) const {
  const std::shared_ptr<::amrex::EBFArrayBoxFactory>& eb_factory =
      GetPatchHierarchy().GetEmbeddedBoundary(handle.level);
  const ::amrex::FabArray<::amrex::EBCellFlagFab>& flags =
      eb_factory->getMultiEBCellFlagFab();
  const ::amrex::EBCellFlagFab& ffab = flags[*handle.iterator];
  ::amrex::Box box = handle.iterator->growntilebox(gcw);
  return ffab.getType(box);
}

void HyperbolicSplitIntegratorContext::PreAdvanceHierarchy() {}

void HyperbolicSplitIntegratorContext::PostAdvanceHierarchy() {
  PatchHierarchy& hierarchy = GetPatchHierarchy();
  const int nlevels = hierarchy.GetNumberOfLevels();
  for (int level = 0; level < nlevels; ++level) {
    hierarchy.GetPatchLevel(level).time_point =
        GetTimePoint(level, Direction::X);
    hierarchy.GetPatchLevel(level).cycles = GetCycles(level, Direction::X);
  }
}

} // namespace cutcell
} // namespace amrex
} // namespace fub
