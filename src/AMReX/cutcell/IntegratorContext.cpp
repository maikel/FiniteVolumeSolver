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

#include "fub/AMReX/cutcell/IntegratorContext.hpp"

#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/cutcell/IndexSpace.hpp"

#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxReg_C.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>

namespace fub::amrex::cutcell {

////////////////////////////////////////////////////////////////////////////////
//                                                  IntegratorContext::LevelData

////////////////////////////////////////////////////////////////////////////////
//                                                      Move Assignment Operator

IntegratorContext::LevelData& IntegratorContext::LevelData::
operator=(LevelData&& other) noexcept {
  if (other.coarse_fine.fineLevel() > 0) {
    // If we do not invoke clear in beforehand it will throw an error in AMReX
    coarse_fine.clear();
    coarse_fine.define(scratch[0].boxArray(), scratch[0].DistributionMap(),
                       other.coarse_fine.refRatio(),
                       other.coarse_fine.fineLevel(),
                       other.coarse_fine.nComp());
  }
  eb_factory = std::move(other.eb_factory);
  reference_states = std::move(other.reference_states);
  boundary_fluxes = std::move(other.boundary_fluxes);
  scratch = std::move(other.scratch);
  fluxes = std::move(other.fluxes);
  stabilized_fluxes = std::move(other.stabilized_fluxes);
  shielded_left_fluxes = std::move(other.shielded_left_fluxes);
  shielded_right_fluxes = std::move(other.shielded_right_fluxes);
  doubly_shielded_fluxes = std::move(other.doubly_shielded_fluxes);
  time_point = other.time_point;
  regrid_time_point = other.regrid_time_point;
  cycles = other.cycles;
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
//                                                             IntegratorContext

////////////////////////////////////////////////////////////////////////////////
//                                          Constructor and Assignment Operators

IntegratorContext::IntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, HyperbolicMethod nm)
    : ghost_cell_width_{nm.flux_method.GetStencilWidth() + 1},
      gridding_{std::move(gridding)}, data_{}, method_{std::move(nm)} {
  data_.reserve(
      static_cast<std::size_t>(GetPatchHierarchy().GetMaxNumberOfLevels()));
  // Allocate auxiliary data arrays for each refinement level in the hierarchy
  ResetHierarchyConfiguration();
}

IntegratorContext::IntegratorContext(const IntegratorContext& other)
    : ghost_cell_width_{other.ghost_cell_width_}, gridding_{other.gridding_},
      data_(static_cast<std::size_t>(GetPatchHierarchy().GetNumberOfLevels())),
      method_{other.method_} {
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

IntegratorContext& IntegratorContext::IntegratorContext::
operator=(const IntegratorContext& other) {
  // We use the copy and move idiom to provide the strong exception guarantee.
  // If an exception occurs we do not change the original object.
  IntegratorContext tmp{other};
  return (*this = std::move(tmp));
}

///////////////////////////////////////////////////////////////////////////////
//                                                             Member Accessors

const BoundaryCondition&
IntegratorContext::GetBoundaryCondition(int level) const {
  return gridding_->GetBoundaryCondition(level);
}

BoundaryCondition& IntegratorContext::GetBoundaryCondition(int level) {
  return gridding_->GetBoundaryCondition(level);
}

const std::shared_ptr<GriddingAlgorithm>&
IntegratorContext::GetGriddingAlgorithm() const noexcept {
  return gridding_;
}

PatchHierarchy& IntegratorContext::GetPatchHierarchy() noexcept {
  return gridding_->GetPatchHierarchy();
}

const PatchHierarchy& IntegratorContext::GetPatchHierarchy() const noexcept {
  return gridding_->GetPatchHierarchy();
}

MPI_Comm IntegratorContext::GetMpiCommunicator() const noexcept {
  return ::amrex::ParallelContext::CommunicatorAll();
}

::amrex::MultiFab& IntegratorContext::GetData(int level) {
  return GetPatchHierarchy().GetPatchLevel(level).data;
}

::amrex::MultiFab& IntegratorContext::GetScratch(int level, Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].scratch[d];
}

::amrex::MultiFab& IntegratorContext::GetFluxes(int level, Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].fluxes[d];
}

::amrex::MultiCutFab& IntegratorContext::GetBoundaryFluxes(int level) {
  const std::size_t l = static_cast<std::size_t>(level);
  return *data_[l].boundary_fluxes;
}

::amrex::MultiFab& IntegratorContext::GetReferenceStates(int level) {
  const std::size_t l = static_cast<std::size_t>(level);
  return *data_[l].reference_states;
}

::amrex::MultiFab& IntegratorContext::GetShieldedFromLeftFluxes(int level,
                                                                Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].shielded_left_fluxes[d];
}

::amrex::MultiFab&
IntegratorContext::GetShieldedFromRightFluxes(int level, Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].shielded_right_fluxes[d];
}

::amrex::MultiFab& IntegratorContext::GetDoublyShieldedFluxes(int level,
                                                              Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].doubly_shielded_fluxes[d];
}

::amrex::MultiFab& IntegratorContext::GetStabilizedFluxes(int level,
                                                          Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].stabilized_fluxes[d];
}

const ::amrex::Geometry& IntegratorContext::GetGeometry(int level) const {
  return GetPatchHierarchy().GetGeometry(level);
}

Duration IntegratorContext::GetTimePoint(int level, Direction dir) const {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].time_point[d];
}

std::ptrdiff_t IntegratorContext::GetCycles(int level, Direction dir) const {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].cycles[d];
}

::amrex::EBFArrayBoxFactory& IntegratorContext::GetEmbeddedBoundary(int level) {
  return *GetPatchHierarchy().GetEmbeddedBoundary(level);
}

const ::amrex::EBFArrayBoxFactory&
IntegratorContext::GetEmbeddedBoundary(int level) const {
  return *GetPatchHierarchy().GetEmbeddedBoundary(level);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                    Observers

bool IntegratorContext::LevelExists(int level) const noexcept {
  return 0 <= level && level < GetPatchHierarchy().GetNumberOfLevels();
}

double IntegratorContext::GetDx(int level, Direction dir) const noexcept {
  return GetPatchHierarchy().GetGeometry(level).CellSize(int(dir));
}

int IntegratorContext::GetRatioToCoarserLevel(int level, Direction dir) const
    noexcept {
  return GetPatchHierarchy().GetRatioToCoarserLevel(level, dir);
}

::amrex::IntVect IntegratorContext::GetRatioToCoarserLevel(int level) const
    noexcept {
  return GetPatchHierarchy().GetRatioToCoarserLevel(level);
}

::amrex::FabType IntegratorContext::GetFabType(int level,
                                               const ::amrex::MFIter& mfi) const
    noexcept {
  return GetPatchHierarchy()
      .GetEmbeddedBoundary(level)
      ->getMultiEBCellFlagFab()[mfi]
      .getType();
}

///////////////////////////////////////////////////////////////////////////////
//                                                                    Modifiers

void IntegratorContext::ResetHierarchyConfiguration(
    std::shared_ptr<GriddingAlgorithm> gridding) {
  gridding_ = std::move(gridding);
  ResetHierarchyConfiguration();
}

void IntegratorContext::ResetHierarchyConfiguration(int first_level) {
  const int n_components =
      GetPatchHierarchy().GetDataDescription().n_state_components;
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
    const std::shared_ptr<::amrex::EBFArrayBoxFactory>& ebf =
        GetPatchHierarchy().GetEmbeddedBoundary(level);

    { // Redistribute reference states
      ::amrex::MultiFab refs(ba, dm, n_components, ghost_cell_width_,
                             ::amrex::MFInfo(), *ebf);
      if (data.reference_states) {
        refs.ParallelCopy(*data.reference_states, 0, 0, n_components,
                          ghost_cell_width_, ghost_cell_width_);
      }
      data.reference_states = std::move(refs);
    }

    data.boundary_fluxes = std::make_unique<::amrex::MultiCutFab>(
        ba, dm, n_cons_components, ghost_cell_width_,
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

void IntegratorContext::SetTimePoint(Duration dt, int level, Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  data_[l].time_point[d] = dt;
}

void IntegratorContext::SetCycles(std::ptrdiff_t cycles, int level,
                                  Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  data_[l].cycles[d] = cycles;
}

void IntegratorContext::FillGhostLayerTwoLevels(int fine, int coarse,
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
  ::amrex::FillPatchTwoLevels(scratch, ft[0], *GetPatchHierarchy().GetOptions().index_spaces[fine], cmf, ct, fmf, ft, 0, 0, nc, cgeom,
                              fgeom, coarse_condition, 0, fine_condition, 0,
                              ratio, mapper, bcr, 0, ::amrex::NullInterpHook(), ::amrex::NullInterpHook());
}

  void IntegratorContext::FillGhostLayerSingleLevel(int level, BoundaryCondition& bc, Direction dir) {
    ::amrex::MultiFab& scratch = GetScratch(level, dir);
    ::amrex::Vector<::amrex::BCRec> bcr(
                                        static_cast<std::size_t>(scratch.nComp()));
    const int nc = scratch.nComp();
    const ::amrex::Vector<::amrex::MultiFab*> smf{&GetData(level)};
    const ::amrex::Vector<double> stime{GetTimePoint(level, dir).count()};
    const ::amrex::Geometry& geom = GetGeometry(level);
    ::amrex::FillPatchSingleLevel(scratch, stime[0], smf, stime, 0, 0, nc, geom,
                                  bc, 0);
  }

void IntegratorContext::FillGhostLayerSingleLevel(int level, Direction dir) {
  BoundaryCondition& condition = GetBoundaryCondition(level);
  FillGhostLayerSingleLevel(level, condition, dir);
}

//namespace {
//constexpr std::ptrdiff_t ipow(int base, int exponent) {
//  std::ptrdiff_t prod{1};
//  while (exponent > 0) {
//    prod *= base;
//    exponent -= 1;
//  }
//  return prod;
//}
//} // namespace

void IntegratorContext::CoarsenConservatively(int fine_level, int coarse_level,
                                              Direction dir) {
  const int first =
      GetPatchHierarchy().GetDataDescription().first_cons_component;
  const int size = GetPatchHierarchy().GetDataDescription().n_cons_components;
  ::amrex::EB_average_down(GetScratch(fine_level, dir),
                           GetScratch(coarse_level, dir), first, size, 2);
}

void IntegratorContext::AccumulateCoarseFineFluxes(int level, Duration,
                                                   Direction dir) {
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

void IntegratorContext::ResetCoarseFineFluxes(int fine, int coarse,
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

void IntegratorContext::ApplyFluxCorrection(int fine, int coarse, Duration,
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

Duration IntegratorContext::ComputeStableDt(int level, Direction dir) {
  return method_.flux_method.ComputeStableDt(*this, level, dir);
}

void IntegratorContext::ComputeNumericFluxes(int level, Duration dt,
                                             Direction dir) {
  method_.flux_method.ComputeNumericFluxes(*this, level, dt, dir);
}

void IntegratorContext::CompleteFromCons(int level, Duration dt,
                                         Direction dir) {
  method_.reconstruction.CompleteFromCons(*this, level, dt, dir);
}

void IntegratorContext::UpdateConservatively(int level, Duration dt,
                                             Direction dir) {
  method_.time_integrator.UpdateConservatively(*this, level, dt, dir);
}

void IntegratorContext::PreAdvanceLevel(int level_num, Direction dir, Duration,
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
IntegratorContext::PostAdvanceLevel(int level_num, Direction dir, Duration dt,
                                    int) {
  SetCycles(GetCycles(level_num, dir) + 1, level_num, dir);
  SetTimePoint(GetTimePoint(level_num, dir) + dt, level_num, dir);
  return boost::outcome_v2::success();
}

void IntegratorContext::PreAdvanceHierarchy() {
  method_.flux_method.PreAdvanceHierarchy(*this);
}

void IntegratorContext::PostAdvanceHierarchy() {
  PatchHierarchy& hierarchy = GetPatchHierarchy();
  const int nlevels = hierarchy.GetNumberOfLevels();
  for (int level = 0; level < nlevels; ++level) {
    hierarchy.GetPatchLevel(level).time_point =
        GetTimePoint(level, Direction::X);
    hierarchy.GetPatchLevel(level).cycles = GetCycles(level, Direction::X);
  }
}

} // namespace fub::amrex::cutcell
