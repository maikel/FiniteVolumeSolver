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
#include "fub/AMReX/ViewFArrayBox.hpp"

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
      grow[i] = ghost_cell_width_;
    }
    data.scratch.define(ba, dm, n_comp, grow);
    // We need to allocate data for each dimensions because of
    // amrex::FluxBoundary::OverwriteFlux.
    for (std::size_t d = 0; d < static_cast<std::size_t>(AMREX_SPACEDIM); ++d) {
      const ::amrex::BoxArray fba =
          ::amrex::convert(ba, ::amrex::IntVect::TheDimensionVector(int(d)));
      ::amrex::IntVect fgrow = grow;
      fgrow[int(d)] = 1;
      //          ::amrex::IntVect::TheDimensionVector(int(d));
      data.fluxes[d].define(fba, dm, n_cons_components, fgrow);
    }
    if (level > 0) {
      const ::amrex::IntVect ref_ratio =
          GetPatchHierarchy().GetRatioToCoarserLevel(level);
      data.coarse_fine.clear();
      data.coarse_fine.define(ba, dm, ref_ratio, level, n_cons_components);
    }
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

void IntegratorContext::FillGhostLayerTwoLevels(
    int fine, BoundaryCondition& fine_condition, int coarse,
    BoundaryCondition& coarse_condition) {
  FUB_ASSERT(coarse >= 0 && fine > coarse);
  ::amrex::MultiFab& scratch = GetScratch(fine);
  ::amrex::Vector<::amrex::BCRec> bcr(
      static_cast<std::size_t>(scratch.nComp()));
  const int nc = scratch.nComp();
  const ::amrex::Vector<::amrex::MultiFab*> cmf{&GetData(coarse)};
  const ::amrex::Vector<::amrex::MultiFab*> fmf{&GetData(fine)};
  const ::amrex::Vector<double> ct{GetTimePoint(coarse).count()};
  const ::amrex::Vector<double> ft{GetTimePoint(fine).count()};
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  const ::amrex::Geometry& fgeom = GetGeometry(fine);
  const ::amrex::IntVect ratio =
      GetGriddingAlgorithm()->GetPatchHierarchy().GetOptions().refine_ratio;
  ::amrex::Interpolater* mapper = &::amrex::pc_interp;
  ::amrex::FillPatchTwoLevels(scratch, ft[0], cmf, ct, fmf, ft, 0, 0, nc, cgeom,
                              fgeom, coarse_condition, 0, fine_condition, 0,
                              ratio, mapper, bcr, 0);
}

void IntegratorContext::FillGhostLayerTwoLevels(int fine, int coarse) {
  return FillGhostLayerTwoLevels(fine, GetBoundaryCondition(fine), coarse,
                                 GetBoundaryCondition(coarse));
}

void IntegratorContext::FillGhostLayerSingleLevel(int level,
                                                  BoundaryCondition& bc) {
  ::amrex::MultiFab& scratch = GetScratch(level);
  ::amrex::Vector<::amrex::BCRec> bcr(
      static_cast<std::size_t>(scratch.nComp()));
  const int nc = scratch.nComp();
  const ::amrex::Vector<::amrex::MultiFab*> smf{&GetData(level)};
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

void IntegratorContext::AccumulateCoarseFineFluxes(int level, double time_scale,
                                                   Direction dir) {
  if (level > 0) {
    const ::amrex::MultiFab& fluxes = GetFluxes(level, dir);
    const int dim = GetPatchHierarchy().GetDataDescription().dimension;
    const double space_scale_inv =
        static_cast<double>(ipow(GetRatioToCoarserLevel(level, dir), dim));
    const double total_scale = time_scale / space_scale_inv;
    data_[static_cast<std::size_t>(level)].coarse_fine.FineAdd(
        fluxes, int(dir), 0, 0, fluxes.nComp(), total_scale);
  }
}

void IntegratorContext::ResetCoarseFineFluxes(int fine, int coarse) {
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    std::size_t sfine = static_cast<std::size_t>(fine);
    data_[sfine].coarse_fine.ClearInternalBorders(GetGeometry(coarse));
    const ::amrex::MultiFab& flux = GetFluxes(coarse, Direction(dir));
    const ::amrex::BoxArray& boxes = flux.boxArray();
    const ::amrex::DistributionMapping& distribution = flux.DistributionMap();
    const int ncomp = flux.nComp();
    ::amrex::MultiFab zero(boxes, distribution, ncomp, 0);
    zero.setVal(0.0);
    data_[sfine].coarse_fine.CrseInit(zero, dir, 0, 0, ncomp);
  }
}

void IntegratorContext::ApplyFluxCorrection(int fine, int coarse, Duration) {
  const std::size_t sfine = static_cast<std::size_t>(fine);
  const int ncomp = GetPatchHierarchy().GetDataDescription().n_cons_components;
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  std::array<::amrex::MultiFab*, AMREX_SPACEDIM> crse_fluxes{AMREX_D_DECL(
      &GetFluxes(coarse, Direction::X), &GetFluxes(coarse, Direction::Y),
      &GetFluxes(coarse, Direction::Z))};
  data_[sfine].coarse_fine.OverwriteFlux(crse_fluxes, 1.0, 0, 0, ncomp, cgeom);
}

Duration IntegratorContext::ComputeStableDt(int level, Direction dir) {
  return method_.flux_method.ComputeStableDt(*this, level, dir);
}

int IntegratorContext::Rank() const noexcept {
  return GetPatchHierarchy().GetDataDescription().dimension;
}

void IntegratorContext::ComputeNumericFluxes(int level, Duration dt,
                                             Direction dir) {
  method_.flux_method.ComputeNumericFluxes(*this, level, dt, dir);
}

void IntegratorContext::CompleteFromCons(int level, Duration dt) {
  method_.reconstruction.CompleteFromCons(*this, level, dt);
}

void IntegratorContext::UpdateConservatively(int level, Duration dt,
                                             Direction dir) {
  method_.time_integrator.UpdateConservatively(*this, level, dt, dir);
}

void IntegratorContext::PreAdvanceLevel(int level_num, Duration, int subcycle) {
  const std::size_t l = static_cast<std::size_t>(level_num);
  if (subcycle == 0 && data_[l].regrid_time_point != data_[l].time_point) {
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

Result<void, TimeStepTooLarge>
IntegratorContext::PostAdvanceLevel(int level_num, Duration dt, int) {
  SetCycles(GetCycles(level_num) + 1, level_num);
  SetTimePoint(GetTimePoint(level_num) + dt, level_num);
  PatchLevel& level = GetPatchHierarchy().GetPatchLevel(level_num);
  level.time_point = GetTimePoint(level_num);
  level.cycles = GetCycles(level_num);
  return boost::outcome_v2::success();
}

void IntegratorContext::PostAdvanceHierarchy() {
  PatchHierarchy& hierarchy = GetPatchHierarchy();
  int nlevels = hierarchy.GetNumberOfLevels();
  for (int level = 0; level < nlevels; ++level) {
    double timepoint = GetTimePoint(level).count();
    ::MPI_Bcast(&timepoint, 1, MPI_DOUBLE, 0, GetMpiCommunicator());
    int cycles = GetCycles(level);
    ::MPI_Bcast(&cycles, 1, MPI_INT, 0, GetMpiCommunicator());
    hierarchy.GetPatchLevel(level).time_point = Duration(timepoint);
    hierarchy.GetPatchLevel(level).cycles = cycles;
  }
}

} // namespace fub::amrex
