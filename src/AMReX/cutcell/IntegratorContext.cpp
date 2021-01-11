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

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/boundary_condition/BoundaryConditionRef.hpp"
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

IntegratorContext::LevelData&
IntegratorContext::LevelData::operator=(LevelData&& other) noexcept {
  if (other.coarse_fine.fineLevel() > 0) {
    // If we do not invoke clear in beforehand it will throw an error in AMReX
    coarse_fine.clear();
    coarse_fine.define(scratch.boxArray(), scratch.DistributionMap(),
                       other.coarse_fine.refRatio(),
                       other.coarse_fine.fineLevel(),
                       other.coarse_fine.nComp());
    coarse_fine.setVal(0.0);
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
    std::shared_ptr<GriddingAlgorithm> gridding, HyperbolicMethod hm)
    : IntegratorContext(std::move(gridding), hm,
                        hm.flux_method.GetStencilWidth() * 4,
                        hm.flux_method.GetStencilWidth() * 4) {}

IntegratorContext::IntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, HyperbolicMethod hm,
    int scratch_gcw, int flux_gcw)
    : scratch_ghost_cell_width_{scratch_gcw}, flux_ghost_cell_width_{flux_gcw},
      gridding_{std::move(gridding)}, data_{}, method_{std::move(hm)} {
  data_.reserve(
      static_cast<std::size_t>(GetPatchHierarchy().GetMaxNumberOfLevels()));
  // Allocate auxiliary data arrays for each refinement level in the hierarchy
  ResetHierarchyConfiguration();
  std::size_t n_levels =
      static_cast<std::size_t>(GetPatchHierarchy().GetNumberOfLevels());
  for (std::size_t i = 0; i < n_levels; ++i) {
    const auto level = static_cast<int>(i);
    data_[i].cycles =
        gridding_->GetPatchHierarchy().GetPatchLevel(level).cycles;
    data_[i].time_point =
        gridding_->GetPatchHierarchy().GetPatchLevel(level).time_point;
    data_[i].regrid_time_point = data_[i].time_point;
  }
}

IntegratorContext::IntegratorContext(const IntegratorContext& other)
    : scratch_ghost_cell_width_{other.scratch_ghost_cell_width_},
      flux_ghost_cell_width_{other.flux_ghost_cell_width_},
      gridding_{other.gridding_},
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

IntegratorContext& IntegratorContext::IntegratorContext::operator=(
    const IntegratorContext& other) {
  // We use the copy and move idiom to provide the strong exception guarantee.
  // If an exception occurs we do not change the original object.
  IntegratorContext tmp{other};
  return (*this = std::move(tmp));
}

///////////////////////////////////////////////////////////////////////////////
//                                                             Member Accessors

const AnyBoundaryCondition& IntegratorContext::GetBoundaryCondition() const {
  return gridding_->GetBoundaryCondition();
}

AnyBoundaryCondition& IntegratorContext::GetBoundaryCondition() {
  return gridding_->GetBoundaryCondition();
}

const std::shared_ptr<GriddingAlgorithm>&
IntegratorContext::GetGriddingAlgorithm() const noexcept {
  return gridding_;
}

const std::shared_ptr<CounterRegistry>&
IntegratorContext::GetCounterRegistry() const noexcept {
  return GetPatchHierarchy().GetCounterRegistry();
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

::amrex::MultiFab& IntegratorContext::GetScratch(int level) {
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].scratch;
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

Duration IntegratorContext::GetTimePoint(int level) const {
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].time_point;
}

std::ptrdiff_t IntegratorContext::GetCycles(int level) const {
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].cycles;
}

::amrex::EBFArrayBoxFactory& IntegratorContext::GetEmbeddedBoundary(int level) {
  return *GetPatchHierarchy().GetEmbeddedBoundary(level);
}

const ::amrex::EBFArrayBoxFactory&
IntegratorContext::GetEmbeddedBoundary(int level) const {
  return *GetPatchHierarchy().GetEmbeddedBoundary(level);
}

const HyperbolicMethod& IntegratorContext::GetHyperbolicMethod() const
    noexcept {
  return method_;
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

void IntegratorContext::CopyDataToScratch(int level_num) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::CopyDataToScratch");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::CopyDataToScratch({})", level_num));
  }
  gridding_->FillMultiFabFromLevel(GetScratch(level_num), level_num);
}

void IntegratorContext::CopyScratchToData(int level_num) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::CopyScratchToData");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::CopyScratchToData({})", level_num));
  }
  GetData(level_num).ParallelCopy(
      GetScratch(level_num),
      GetPatchHierarchy().GetGeometry(level_num).periodicity());
}

void IntegratorContext::ResetHierarchyConfiguration(
    std::shared_ptr<GriddingAlgorithm> gridding) {
  gridding_ = std::move(gridding);
  ResetHierarchyConfiguration();
}

void IntegratorContext::ResetHierarchyConfiguration(int first_level) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::ResetHierarchyConfiguration");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::ResetHierarchyConfiguration({})",
        first_level));
  }
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
      ::amrex::MultiFab refs(ba, dm, n_components, scratch_ghost_cell_width_,
                             ::amrex::MFInfo(), *ebf);
      if (data.reference_states) {
        refs.ParallelCopy(*data.reference_states, 0, 0, n_components,
                          scratch_ghost_cell_width_, scratch_ghost_cell_width_);
      }
      data.reference_states = std::move(refs);
    }

    data.boundary_fluxes = std::make_unique<::amrex::MultiCutFab>(
        ba, dm, n_cons_components, scratch_ghost_cell_width_,
        ebf->getMultiEBCellFlagFab());

    data.scratch.define(ba, dm, n_components, scratch_ghost_cell_width_);
    ::amrex::IntVect grow(scratch_ghost_cell_width_);

    for (std::size_t d = 0; d < static_cast<std::size_t>(AMREX_SPACEDIM); ++d) {
      ::amrex::IntVect fgrow(flux_ghost_cell_width_ + method_.flux_method.GetStencilWidth());
      fgrow[int(d)] = flux_ghost_cell_width_;
      const ::amrex::IntVect unit =
          ::amrex::IntVect::TheDimensionVector(int(d));
      data.fluxes[d].define(::amrex::convert(ba, unit), dm, n_cons_components,
                            fgrow);
      data.fluxes[d].setVal(0.0);
      data.stabilized_fluxes[d].define(::amrex::convert(ba, unit), dm,
                                       n_cons_components, fgrow);
      data.stabilized_fluxes[d].setVal(0.0);
      data.shielded_left_fluxes[d].define(::amrex::convert(ba, unit), dm,
                                          n_cons_components, fgrow);
      data.shielded_left_fluxes[d].setVal(0.0);
      data.shielded_right_fluxes[d].define(::amrex::convert(ba, unit), dm,
                                           n_cons_components, fgrow);
      data.shielded_right_fluxes[d].setVal(0.0);
      data.doubly_shielded_fluxes[d].define(::amrex::convert(ba, unit), dm,
                                            n_cons_components, fgrow);
      data.doubly_shielded_fluxes[d].setVal(0.0);
    }

    data.eb_factory = ebf;
    if (level > 0) {
      const ::amrex::IntVect ref_ratio = GetRatioToCoarserLevel(level);
      data.coarse_fine.clear();
      data.coarse_fine.define(ba, dm, ref_ratio, level, n_cons_components);
      data.coarse_fine.setVal(0.0);
    }
  }
  for (int level_num = first_level; level_num < static_cast<int>(data_.size());
       ++level_num) {
    CopyDataToScratch(level_num);
  }
}

void IntegratorContext::SetTimePoint(Duration dt, int level) {
  const std::size_t l = static_cast<std::size_t>(level);
  data_[l].time_point = dt;
}

void IntegratorContext::SetCycles(std::ptrdiff_t cycles, int level) {
  const std::size_t l = static_cast<std::size_t>(level);
  data_[l].cycles = cycles;
}

void IntegratorContext::ApplyBoundaryCondition(int level, Direction dir) {
  AnyBoundaryCondition& boundary_condition = GetBoundaryCondition();
  ApplyBoundaryCondition(level, dir, boundary_condition);
}

void IntegratorContext::ApplyBoundaryCondition(int level, Direction dir,
                                               AnyBoundaryCondition& bc) {
  Timer timer = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::ApplyBoundaryCondition");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::ApplyBoundaryCondition({})", level));
  }
  ::amrex::MultiFab& scratch = GetScratch(level);
  GriddingAlgorithm& grid = *GetGriddingAlgorithm();
  bc.FillBoundary(scratch, grid, level, dir);
}

void IntegratorContext::FillGhostLayerTwoLevels(int fine,
                                                AnyBoundaryCondition& fbc,
                                                int coarse,
                                                AnyBoundaryCondition& cbc) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::FillGhostLayerTwoLevels");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::FillGhostLayerTwoLevels({})", fine));
  }
  FUB_ASSERT(coarse >= 0 && fine > coarse);
  ::amrex::MultiFab& scratch = GetScratch(fine);
  ::amrex::Vector<::amrex::BCRec> bcr(
      static_cast<std::size_t>(scratch.nComp()));
  const int nc = scratch.nComp();
  const ::amrex::Vector<::amrex::MultiFab*> cmf{&GetScratch(coarse)};
  const ::amrex::Vector<::amrex::MultiFab*> fmf{&GetScratch(fine)};
  const ::amrex::Vector<double> ct{GetTimePoint(coarse).count()};
  const ::amrex::Vector<double> ft{GetTimePoint(fine).count()};
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  const ::amrex::Geometry& fgeom = GetGeometry(fine);
  const ::amrex::IntVect ratio = 2 * ::amrex::IntVect::TheUnitVector();
  ::amrex::Interpolater* mapper = &::amrex::pc_interp;
  std::size_t sfine = static_cast<std::size_t>(fine);
  auto&& index_space = *GetPatchHierarchy().GetOptions().index_spaces[sfine];
  BoundaryConditionRef fine_boundary(fbc, *gridding_, fine);
  BoundaryConditionRef coarse_boundary(cbc, *gridding_, coarse);
  ::amrex::FillPatchTwoLevels(scratch, ft[0], index_space, cmf, ct, fmf, ft, 0,
                              0, nc, cgeom, fgeom, coarse_boundary, 0,
                              fine_boundary, 0, ratio, mapper, bcr, 0,
                              ::amrex::NullInterpHook<::amrex::FArrayBox>(),
                              ::amrex::NullInterpHook<::amrex::FArrayBox>());
}

void IntegratorContext::FillGhostLayerTwoLevels(int fine, int coarse) {
  FillGhostLayerTwoLevels(fine, GetBoundaryCondition(), coarse,
                          GetBoundaryCondition());
}

void IntegratorContext::FillGhostLayerSingleLevel(int level,
                                                  AnyBoundaryCondition& bc) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::FillGhostLayerSingleLevel");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::FillGhostLayerSingleLevel({})", level));
  }
  ::amrex::MultiFab& scratch = GetScratch(level);
  ::amrex::Vector<::amrex::BCRec> bcr(
      static_cast<std::size_t>(scratch.nComp()));
  const int nc = scratch.nComp();
  const ::amrex::Vector<::amrex::MultiFab*> smf{&scratch};
  const ::amrex::Vector<double> stime{GetTimePoint(level).count()};
  const ::amrex::Geometry& geom = GetGeometry(level);
  BoundaryConditionRef boundary(bc, *gridding_, level);
  ::amrex::FillPatchSingleLevel(scratch, stime[0], smf, stime, 0, 0, nc, geom,
                                boundary, 0);
}

void IntegratorContext::FillGhostLayerSingleLevel(int level) {
  AnyBoundaryCondition& condition = GetBoundaryCondition();
  FillGhostLayerSingleLevel(level, condition);
}

void IntegratorContext::CoarsenConservatively(int fine_level,
                                              int coarse_level) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::CoarsenConservatively");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("cutcell::IntegratorContext::CoarsenConservatively({}, {})",
                    fine_level, coarse_level));
  }
  const int first =
      GetPatchHierarchy().GetDataDescription().first_cons_component;
  const int size = GetPatchHierarchy().GetDataDescription().n_cons_components;
  ::amrex::EB_average_down(GetScratch(fine_level), GetScratch(coarse_level),
                           first, size, 2);
}

void IntegratorContext::AccumulateCoarseFineFluxes(int level, double scale,
                                                   Direction dir) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::AccumulateCoarseFineFluxes");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::AccumulateCoarseFineFluxes({})", level));
  }
  const ::amrex::MultiFab& fluxes = GetFluxes(level, dir);
  const int dir_v = static_cast<int>(dir);
  const int ncomp = fluxes.nComp();
  const ::amrex::Geometry& geom = GetGeometry(level);
  const double* dx = geom.CellSize();
  const double vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
  const double vol_dt_over_dx = vol * scale / dx[dir_v];
  if (level > 0) {
    const std::size_t slevel = static_cast<std::size_t>(level);
    ::amrex::FluxRegister& flux_register = data_[slevel].coarse_fine;
    flux_register.FineAdd(fluxes, dir_v, 0, 0, ncomp, vol_dt_over_dx);
  }
  if (LevelExists(level + 1)) {
    const std::size_t next_level = static_cast<std::size_t>(level + 1);
    ::amrex::FluxRegister& flux_register = data_[next_level].coarse_fine;
    flux_register.CrseInit(fluxes, dir_v, 0, 0, ncomp, -vol_dt_over_dx,
                           ::amrex::FluxRegister::ADD);
  }
}

void IntegratorContext::ResetCoarseFineFluxes(int fine, int coarse) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::ResetCoarseFineFluxes");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("cutcell::IntegratorContext::ResetCoarseFineFluxes({}, {})",
                    fine, coarse));
  }
  FUB_ASSERT(fine > 0);
  FUB_ASSERT(fine >= coarse);
  FUB_ASSERT(coarse >= 0);
  std::size_t fine_level = static_cast<std::size_t>(fine);
  ::amrex::FluxRegister& flux_register = data_[fine_level].coarse_fine;
  flux_register.ClearInternalBorders(GetGeometry(coarse));
  flux_register.setVal(0.0);
}

void IntegratorContext::ApplyFluxCorrection(
    int fine, int coarse, [[maybe_unused]] Duration time_step_size) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::ApplyFluxCorrection");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("cutcell::IntegratorContext::ApplyFluxCorrection({}, {})",
                    fine, coarse));
  }
  FUB_ASSERT(fine > 0);
  FUB_ASSERT(fine >= coarse);
  FUB_ASSERT(coarse >= 0);
  const std::size_t next_level = static_cast<std::size_t>(fine);
  const int ncomp = GetPatchHierarchy().GetDataDescription().n_cons_components;
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  ::amrex::MultiFab& scratch = GetScratch(coarse);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    ::amrex::FluxRegister& flux_register = data_[next_level].coarse_fine;
    flux_register.Reflux(scratch, dir, 1.0, 0, 0, ncomp, cgeom);
  }
}

Duration IntegratorContext::ComputeStableDt(int level, Direction dir) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::ComputeStableDt");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("cutcell::IntegratorContext::ComputeStableDt({})", level));
  }
  return method_.flux_method.ComputeStableDt(*this, level, dir);
}

void IntegratorContext::ComputeNumericFluxes(int level, Duration dt,
                                             Direction dir) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::ComputeNumericFluxes");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::ComputeNumericFluxes({})", level));
  }
  method_.flux_method.ComputeNumericFluxes(*this, level, dt, dir);
}

void IntegratorContext::CompleteFromCons(int level, Duration dt) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::CompleteFromCons");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("cutcell::IntegratorContext::CompleteFromCons({})", level));
  }
  method_.reconstruction.CompleteFromCons(*this, level, dt);
}

void IntegratorContext::UpdateConservatively(int level, Duration dt,
                                             Direction dir) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::UpdateConservatively");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::UpdateConservatively({})", level));
  }
  method_.time_integrator.UpdateConservatively(*this, level, dt, dir);
}

int IntegratorContext::PreAdvanceLevel(int level_num, Duration,
                                       std::pair<int, int> subcycle) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::PreAdvanceLevel");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::PreAdvanceLevel({})", level_num));
  }
  const std::size_t l = static_cast<std::size_t>(level_num);
  const int max_level = GetPatchHierarchy().GetMaxNumberOfLevels();
  int level_which_changed = max_level;
  if (subcycle.first == 0) {
    if (data_[l].regrid_time_point != data_[l].time_point) {
      level_which_changed = gridding_->RegridAllFinerlevels(level_num);
      for (std::size_t lvl = l; lvl < data_.size(); ++lvl) {
        data_[lvl].regrid_time_point = data_[lvl].time_point;
      }
      if (level_which_changed < max_level) {
        ResetHierarchyConfiguration(level_which_changed);
      }
    }
  }
  return level_which_changed;
}

Result<void, TimeStepTooLarge>
IntegratorContext::PostAdvanceLevel(int level_num, Duration dt,
                                    std::pair<int, int> subcycle) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::PostAdvanceLevel");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "cutcell::IntegratorContext::PostAdvanceLevel({})", level_num));
  }
  SetCycles(GetCycles(level_num) + 1, level_num);
  double timepoint = (GetTimePoint(level_num) + dt).count();
  ::MPI_Bcast(&timepoint, 1, MPI_DOUBLE, 0, GetMpiCommunicator());
  SetTimePoint(Duration(timepoint), level_num);

  if (subcycle.first == subcycle.second - 1) {
    PatchLevel& level = GetPatchHierarchy().GetPatchLevel(level_num);
    level.time_point = GetTimePoint(level_num);
    level.cycles = GetCycles(level_num);
  }

  if (feedback_) {
    feedback_(*this, level_num, dt, subcycle);
  }

  return boost::outcome_v2::success();
}

void IntegratorContext::PreAdvanceHierarchy() {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "cutcell::IntegratorContext::PreAdvanceHierarchy");
  method_.flux_method.PreAdvanceHierarchy(*this);
}

void IntegratorContext::PostAdvanceHierarchy() {
  PatchHierarchy& hierarchy = GetPatchHierarchy();
  int nlevels = hierarchy.GetNumberOfLevels();
  const Duration time_point = GetTimePoint();
  for (int level = 0; level < nlevels; ++level) {
    hierarchy.GetPatchLevel(level).time_point = time_point;
    hierarchy.GetPatchLevel(level).cycles = GetCycles(level);
  }
}

} // namespace fub::amrex::cutcell
