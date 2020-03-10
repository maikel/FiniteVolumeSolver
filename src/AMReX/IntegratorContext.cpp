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

#include "fub/AMReX/IntegratorContext.hpp"

#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/core/algorithm.hpp"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxReg_C.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>

#include <fmt/format.h>

namespace fub::amrex {

////////////////////////////////////////////////////////////////////////////////
//                                   IntegratorContext::LevelData

////////////////////////////////////////////////////////////////////////////////
// Move Assignment Operator

IntegratorContext::LevelData& IntegratorContext::LevelData::
operator=(LevelData&& other) noexcept {
  if (other.coarse_fine.fineLevel() > 0) {
    // If we do not invoke clear in beforehand it will throw an error in AMReX
    coarse_fine.clear();
    coarse_fine.define(scratch.boxArray(), scratch.DistributionMap(),
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
//                                              IntegratorContext

////////////////////////////////////////////////////////////////////////////////
// Constructor and Assignment Operators

// By default we assume Strang-Splitting and subcycling from AMR.
// Each adds a factor of two to the ghost cell width requirements on coarse fine
// boundaries.
IntegratorContext::IntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, HyperbolicMethod hm)
    : cell_ghost_cell_width_{hm.flux_method.GetStencilWidth() * 2 * 2},
      face_ghost_cell_width_{cell_ghost_cell_width_ -
                             hm.flux_method.GetStencilWidth()},
      gridding_{std::move(gridding)}, data_{}, method_{std::move(hm)} {
  data_.reserve(
      static_cast<std::size_t>(GetPatchHierarchy().GetMaxNumberOfLevels()));
  // Allocate auxiliary data arrays for each refinement level in the hierarchy
  ResetHierarchyConfiguration();
  std::size_t n_levels =
      static_cast<std::size_t>(GetPatchHierarchy().GetNumberOfLevels());
  for (std::size_t i = 0; i < n_levels; ++i) {
    data_[i].cycles = gridding_->GetPatchHierarchy().GetPatchLevel(i).cycles;
    data_[i].time_point =
        gridding_->GetPatchHierarchy().GetPatchLevel(i).time_point;
    data_[i].regrid_time_point = data_[i].time_point;
  }
}

IntegratorContext::IntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, HyperbolicMethod nm,
    int cell_gcw, int face_gcw)
    : cell_ghost_cell_width_{cell_gcw}, face_ghost_cell_width_{face_gcw},
      gridding_{std::move(gridding)}, data_{}, method_{std::move(nm)} {
  data_.reserve(
      static_cast<std::size_t>(GetPatchHierarchy().GetMaxNumberOfLevels()));
  // Allocate auxiliary data arrays for each refinement level in the hierarchy
  ResetHierarchyConfiguration();
  std::size_t n_levels =
      static_cast<std::size_t>(GetPatchHierarchy().GetNumberOfLevels());
  for (std::size_t i = 0; i < n_levels; ++i) {
    data_[i].cycles = gridding_->GetPatchHierarchy().GetPatchLevel(i).cycles;
    data_[i].time_point =
        gridding_->GetPatchHierarchy().GetPatchLevel(i).time_point;
    data_[i].regrid_time_point = data_[i].time_point;
  }
}

IntegratorContext::IntegratorContext(const IntegratorContext& other)
    : cell_ghost_cell_width_{other.cell_ghost_cell_width_},
      face_ghost_cell_width_{other.face_ghost_cell_width_},
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

IntegratorContext& IntegratorContext::IntegratorContext::
operator=(const IntegratorContext& other) {
  // We use the copy and move idiom to provide the strong exception guarantee.
  // If an exception occurs we do not change the original object.
  IntegratorContext tmp{other};
  return (*this = std::move(tmp));
}

///////////////////////////////////////////////////////////////////////////////
//                                                            Member Accessors

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

const std::shared_ptr<CounterRegistry>&
IntegratorContext::GetCounterRegistry() const noexcept {
  return GetPatchHierarchy().GetCounterRegistry();
}

::amrex::MultiFab& IntegratorContext::GetData(int level) {
  return GetPatchHierarchy().GetPatchLevel(level).data;
}

const ::amrex::MultiFab& IntegratorContext::GetData(int level) const {
  return GetPatchHierarchy().GetPatchLevel(level).data;
}

::amrex::MultiFab& IntegratorContext::GetScratch(int level) {
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].scratch;
}

const ::amrex::MultiFab& IntegratorContext::GetScratch(int level) const {
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].scratch;
}

::amrex::MultiFab& IntegratorContext::GetFluxes(int level, Direction dir) {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].fluxes[d];
}

const ::amrex::MultiFab& IntegratorContext::GetFluxes(int level,
                                                      Direction dir) const {
  const std::size_t d = static_cast<std::size_t>(dir);
  const std::size_t l = static_cast<std::size_t>(level);
  return data_[l].fluxes[d];
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

///////////////////////////////////////////////////////////////////////////////
//                                                                    Observers

bool IntegratorContext::LevelExists(int level) const noexcept {
  return 0 <= level && level < GetPatchHierarchy().GetNumberOfLevels();
}

int IntegratorContext::GetRatioToCoarserLevel(int level, Direction dir) const
    noexcept {
  return GetPatchHierarchy().GetRatioToCoarserLevel(level, dir);
}

::amrex::IntVect IntegratorContext::GetRatioToCoarserLevel(int level) const
    noexcept {
  return GetPatchHierarchy().GetRatioToCoarserLevel(level);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                    Modifiers

void IntegratorContext::ResetHierarchyConfiguration(
    std::shared_ptr<GriddingAlgorithm> gridding) {
  gridding_ = std::move(gridding);
  ResetHierarchyConfiguration();
}

void IntegratorContext::ResetHierarchyConfiguration(int first_level) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "IntegratorContext::ResetHierarchyConfiguration");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "IntegratorContext::ResetHierarchyConfiguration({})", first_level));
  }
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
    ::amrex::IntVect grow = ::amrex::IntVect::TheZeroVector();
    for (int i = 0; i < GetPatchHierarchy().GetDataDescription().dimension;
         ++i) {
      grow[i] = cell_ghost_cell_width_;
    }
    data.scratch.define(ba, dm, n_comp, grow);
    for (std::size_t d = 0; d < static_cast<std::size_t>(AMREX_SPACEDIM); ++d) {
      const ::amrex::BoxArray fba =
          ::amrex::convert(ba, ::amrex::IntVect::TheDimensionVector(int(d)));
      ::amrex::IntVect fgrow = grow;
      fgrow[int(d)] = face_ghost_cell_width_;
      data.fluxes[d].define(fba, dm, n_cons_components, fgrow);
      data.fluxes[d].setVal(0.0);
    }
    if (level > 0) {
      const ::amrex::IntVect ref_ratio =
          GetPatchHierarchy().GetRatioToCoarserLevel(level);
      data.coarse_fine.clear();
      data.coarse_fine.define(ba, dm, ref_ratio, level, n_cons_components);
    }
  }
  for (std::size_t level_num = first_level; level_num < data_.size();
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
  BoundaryCondition& boundary_condition = GetBoundaryCondition(level);
  ApplyBoundaryCondition(level, dir, boundary_condition);
}

void IntegratorContext::ApplyBoundaryCondition(int level, Direction dir,
                                               BoundaryCondition& bc) {
  Timer timer = GetCounterRegistry()->get_timer(
      "IntegratorContext::ApplyBoundaryCondition");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::ApplyBoundaryCondition({})", level));
  }
  ::amrex::MultiFab& scratch = GetScratch(level);
  Duration time_point = GetTimePoint(level);
  const ::amrex::Geometry& geometry = GetGeometry(level);
  GriddingAlgorithm& grid = *GetGriddingAlgorithm();
  bc.FillBoundary(scratch, geometry, time_point, grid, dir);
}

void IntegratorContext::FillGhostLayerTwoLevels(
    int fine, BoundaryCondition& fine_condition, int coarse,
    BoundaryCondition& coarse_condition) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "IntegratorContext::FillGhostLayerTwoLevels");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::FillGhostLayerTwoLevels({})", fine));
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
  const ::amrex::IntVect ratio =
      GetGriddingAlgorithm()->GetPatchHierarchy().GetOptions().refine_ratio;
  ::amrex::Interpolater* mapper = &::amrex::pc_interp;
#ifdef AMREX_USE_EB
  auto index_space = GetPatchHierarchy().GetIndexSpaces();
  ::amrex::FillPatchTwoLevels(
      scratch, ft[0], *index_space[fine], cmf, ct, fmf, ft, 0, 0, nc, cgeom,
      fgeom, coarse_condition, 0, fine_condition, 0, ratio, mapper, bcr, 0,
      ::amrex::NullInterpHook(), ::amrex::NullInterpHook());
#else
  ::amrex::FillPatchTwoLevels(scratch, ft[0], cmf, ct, fmf, ft, 0, 0, nc, cgeom,
                              fgeom, coarse_condition, 0, fine_condition, 0,
                              ratio, mapper, bcr, 0);
#endif
}

void IntegratorContext::FillGhostLayerTwoLevels(int fine, int coarse) {
  return FillGhostLayerTwoLevels(fine, GetBoundaryCondition(fine), coarse,
                                 GetBoundaryCondition(coarse));
}

void IntegratorContext::FillGhostLayerSingleLevel(int level,
                                                  BoundaryCondition& bc) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "IntegratorContext::FillGhostLayerSingleLevel");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::FillGhostLayerSingleLevel({})", level));
  }
  ::amrex::MultiFab& scratch = GetScratch(level);
  ::amrex::Vector<::amrex::BCRec> bcr(
      static_cast<std::size_t>(scratch.nComp()));
  const int nc = scratch.nComp();
  const ::amrex::Vector<::amrex::MultiFab*> smf{&scratch};
  const ::amrex::Vector<double> stime{GetTimePoint(level).count()};
  const ::amrex::Geometry& geom = GetGeometry(level);
  ::amrex::FillPatchSingleLevel(scratch, stime[0], smf, stime, 0, 0, nc, geom,
                                bc, 0);
}

void IntegratorContext::FillGhostLayerSingleLevel(int level) {
  FillGhostLayerSingleLevel(level, GetBoundaryCondition(level));
}

void IntegratorContext::CoarsenConservatively(int fine_level,
                                              int coarse_level) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "IntegratorContext::CoarsenConservatively");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::CoarsenConservatively({}, {})",
                    fine_level, coarse_level));
  }
  const int first =
      GetPatchHierarchy().GetDataDescription().first_cons_component;
  const int size = GetPatchHierarchy().GetDataDescription().n_cons_components;
  const ::amrex::MultiFab& fine = GetScratch(fine_level);
  ::amrex::MultiFab& coarse = GetScratch(coarse_level);
  FUB_ASSERT(!fine.contains_nan());
  FUB_ASSERT(!coarse.contains_nan());
  ::amrex::average_down(fine, coarse, GetGeometry(fine_level),
                        GetGeometry(coarse_level), first, size,
                        GetRatioToCoarserLevel(fine_level));
}

void IntegratorContext::AccumulateCoarseFineFluxes(int level, double scale,
                                                   Direction dir) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "IntegratorContext::AccumulateCoarseFineFluxes");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "IntegratorContext::AccumulateCoarseFineFluxes({})", level));
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
      "IntegratorContext::ResetCoarseFineFluxes");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "IntegratorContext::ResetCoarseFineFluxes({}, {})", fine, coarse));
  }
  std::size_t fine_level = static_cast<std::size_t>(fine);
  ::amrex::FluxRegister& flux_register = data_[fine_level].coarse_fine;
  flux_register.ClearInternalBorders(GetGeometry(coarse));
  flux_register.setVal(0.0);
}

void IntegratorContext::ApplyFluxCorrection(
    int fine, int coarse, [[maybe_unused]] Duration time_step_size) {
  Timer timer1 =
      GetCounterRegistry()->get_timer("IntegratorContext::ApplyFluxCorrection");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(fmt::format(
        "IntegratorContext::ApplyFluxCorrection({}, {})", fine, coarse));
  }
  const std::size_t next_level = static_cast<std::size_t>(fine);
  const int ncomp = GetPatchHierarchy().GetDataDescription().n_cons_components;
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  ::amrex::MultiFab& scratch = GetScratch(coarse);
  for (int dir = 0; dir < Rank(); ++dir) {
    ::amrex::FluxRegister& flux_register = data_[next_level].coarse_fine;
    flux_register.Reflux(scratch, dir, 1.0, 0, 0, ncomp, cgeom);
  }
}

Duration IntegratorContext::ComputeStableDt(int level, Direction dir) {
  Timer timer1 =
      GetCounterRegistry()->get_timer("IntegratorContext::ComputeStableDt");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::ComputeStableDt({})", level));
  }
  return method_.flux_method.ComputeStableDt(*this, level, dir);
}

int IntegratorContext::Rank() const noexcept {
  return GetPatchHierarchy().GetDataDescription().dimension;
}

void IntegratorContext::ComputeNumericFluxes(int level, Duration dt,
                                             Direction dir) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "IntegratorContext::ComputeNumericFluxes");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::ComputeNumericFluxes({})", level));
  }
  method_.flux_method.ComputeNumericFluxes(*this, level, dt, dir);
}

void IntegratorContext::CompleteFromCons(int level, Duration dt) {
  Timer timer1 =
      GetCounterRegistry()->get_timer("IntegratorContext::CompleteFromCons");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::CompleteFromCons({})", level));
  }
  method_.reconstruction.CompleteFromCons(*this, level, dt);
}

void IntegratorContext::UpdateConservatively(int level, Duration dt,
                                             Direction dir) {
  Timer timer1 = GetCounterRegistry()->get_timer(
      "IntegratorContext::UpdateConservatively");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::UpdateConservatively({})", level));
  }
  method_.time_integrator.UpdateConservatively(*this, level, dt, dir);
}

void IntegratorContext::PreAdvanceLevel(int level_num, Duration,
                                        std::pair<int, int> subcycle) {
  Timer timer =
      GetCounterRegistry()->get_timer("IntegratorContext::PreAdvanceLevel");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::PreAdvanceLevel({})", level_num));
  }
  const std::size_t l = static_cast<std::size_t>(level_num);
  if (subcycle.first == 0) {
    if (data_[l].regrid_time_point != data_[l].time_point) {
      gridding_->RegridAllFinerlevels(level_num);
      for (std::size_t lvl = l; lvl < data_.size(); ++lvl) {
        data_[lvl].regrid_time_point = data_[lvl].time_point;
      }
      if (LevelExists(level_num + 1)) {
        ResetHierarchyConfiguration(level_num + 1);
        ResetCoarseFineFluxes(level_num + 1, level_num);
      }
    }
  }
}

void IntegratorContext::CopyDataToScratch(int level_num) {
  Timer timer1 =
      GetCounterRegistry()->get_timer("IntegratorContext::CopyDataToScratch");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::CopyDataToScratch({})", level_num));
  }
  gridding_->FillMultiFabFromLevel(GetScratch(level_num), level_num);
}

void IntegratorContext::CopyScratchToData(int level_num) {
  Timer timer1 =
      GetCounterRegistry()->get_timer("IntegratorContext::CopyScratchToData");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::CopyScratchToData({})", level_num));
  }
  GetData(level_num).ParallelCopy(
      GetScratch(level_num),
      GetPatchHierarchy().GetGeometry(level_num).periodicity());
}

Result<void, TimeStepTooLarge>
IntegratorContext::PostAdvanceLevel(int level_num, Duration dt,
                                    std::pair<int, int> subcycle) {
  Timer timer1 =
      GetCounterRegistry()->get_timer("IntegratorContext::PostAdvanceLevel");
  Timer timer_per_level{};
  if (count_per_level) {
    timer_per_level = GetCounterRegistry()->get_timer(
        fmt::format("IntegratorContext::PostAdvanceLevel({})", level_num));
  }
  SetCycles(GetCycles(level_num) + 1, level_num);
  double timepoint = (GetTimePoint(level_num) + dt).count();
  ::MPI_Bcast(&timepoint, 1, MPI_DOUBLE, 0, GetMpiCommunicator());
  SetTimePoint(Duration(timepoint), level_num);

  // If this is the last AMR subcycle we copy the scratch back to data
  // This essentially means, that the hyperbolic part on this level is
  // done and if any additional operator wants to act on this level it can
  // do so on the original data.
  if (subcycle.first == subcycle.second - 1) {
    PatchLevel& level = GetPatchHierarchy().GetPatchLevel(level_num);
    level.time_point = GetTimePoint(level_num);
    level.cycles = GetCycles(level_num);
  }

  return boost::outcome_v2::success();
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

} // namespace fub::amrex
