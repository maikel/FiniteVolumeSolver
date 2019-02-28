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

#include "fub/grid/AMReX/HyperbolicSplitIntegratorContext.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"
#include "fub/grid/AMReX/utility.hpp"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxReg_C.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>

#include <fmt/format.h>

namespace fub {
namespace amrex {
namespace {
dynamic_mdspan<double, AMREX_SPACEDIM + 1>
GetMdSpan(::amrex::MultiFab& multi_fab, PatchHandle patch) {
  ::amrex::FArrayBox& array = multi_fab[patch.global_index];
  return MakeMdSpan(array);
}
} // namespace

dynamic_mdspan<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetData(PatchHandle patch) {
  return GetMdSpan(GetPatchHierarchy()->GetPatchLevel(patch.level).data, patch);
}

::amrex::MultiFab& HyperbolicSplitIntegratorContext::GetData(int level) {
  return GetPatchHierarchy()->GetPatchLevel(level).data;
}

dynamic_mdspan<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetScratch(PatchHandle patch, Direction dir) {
  const int d = int(dir);
  return GetMdSpan(data_[patch.level].scratch[d], patch);
}

::amrex::MultiFab& HyperbolicSplitIntegratorContext::GetScratch(int level,
                                                                Direction dir) {
  return data_[level].scratch[int(dir)];
}

dynamic_mdspan<double, AMREX_SPACEDIM + 1>
HyperbolicSplitIntegratorContext::GetFluxes(PatchHandle patch, Direction dir) {
  const int d = int(dir);
  return GetMdSpan(data_[patch.level].fluxes[d], patch);
}

::amrex::MultiFab& HyperbolicSplitIntegratorContext::GetFluxes(int level,
                                                               Direction dir) {
  return data_[level].fluxes[int(dir)];
}

const ::amrex::Geometry&
HyperbolicSplitIntegratorContext::GetGeometry(int level) const {
  return GetPatchHierarchy()->GetGeometry(level);
}

Duration HyperbolicSplitIntegratorContext::GetTimePoint(int level,
                                                        Direction dir) const {
  return data_[level].time_point[int(dir)];
}

std::ptrdiff_t
HyperbolicSplitIntegratorContext::GetCycles(int level, Direction dir) const {
  return data_[level].cycles[int(dir)];
}

void HyperbolicSplitIntegratorContext::SetTimePoint(Duration dt, int level,
                                                    Direction dir) {
  data_[level].time_point[int(dir)] = dt;
}

void HyperbolicSplitIntegratorContext::SetCycles(std::ptrdiff_t cycles,
                                                 int level, Direction dir) {
  data_[level].cycles[int(dir)] = cycles;
}

double HyperbolicSplitIntegratorContext::GetDx(PatchHandle patch,
                                               Direction dir) const {
  const int d = int(dir);
  return GetPatchHierarchy()->GetGeometry(patch.level).CellSize(d);
}

HyperbolicSplitIntegratorContext::HyperbolicSplitIntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, int gcw)
    : ghost_cell_width_{gcw + 1}, gridding_{std::move(gridding)},
      data_(GetPatchHierarchy()->GetMaxNumberOfLevels()) {
  ResetHierarchyConfiguration();
}

void HyperbolicSplitIntegratorContext::ResetHierarchyConfiguration(
    int first_level) {
  const int n_cons_components =
      GetPatchHierarchy()->GetDataDescription().n_cons_components;
  const int n_levels = GetPatchHierarchy()->GetMaxNumberOfLevels();
  for (int level = first_level; level < n_levels; ++level) {
    LevelData& data = data_[level];
    const ::amrex::MultiFab& base =
        GetPatchHierarchy()->GetPatchLevel(level).data;
    const ::amrex::BoxArray& ba = base.boxArray();
    const ::amrex::DistributionMapping& dm = base.DistributionMap();
    const int n_comp = base.nComp();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      const ::amrex::IntVect unit = ::amrex::IntVect::TheDimensionVector(d);
      data.scratch[d].define(ba, dm, n_comp, ghost_cell_width_ * unit);
      data.fluxes[d].define(::amrex::convert(ba, unit), dm, n_cons_components,
                            (ghost_cell_width_ - 1) * unit);
    }
    if (level > 0) {
      const ::amrex::IntVect ref_ratio = 2 * ::amrex::IntVect::TheUnitVector();
      data.coarse_fine.clear();
      data.coarse_fine.define(ba, dm, ref_ratio, level, n_cons_components);
    }
  }
}

namespace {
template <typename Function>
struct AdaptBoundaryCondition : public ::amrex::PhysBCFunctBase {
  AdaptBoundaryCondition(Function f, const ::amrex::Geometry& geom, int level)
      : function_{f}, geom_{geom}, level_num_{level} {}

  void FillBoundary(::amrex::MultiFab& mf, int, int, double time_point,
                    int) override {
    if (geom_.isAllPeriodic())
      return;

    //! create a grown domain box containing valid + periodic cells
    const ::amrex::Box& domain = geom_.Domain();
    ::amrex::Box gdomain = ::amrex::convert(domain, mf.boxArray().ixType());
    const ::amrex::IntVect& ngrow = mf.nGrowVect();
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      if (geom_.isPeriodic(i)) {
        gdomain.grow(i, ngrow[i]);
      }
    }
    for (::amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
      ::amrex::Box box = mfi.tilebox();
      if (!gdomain.contains(box)) {
        PatchHandle patch{level_num_, mfi.index()};
        auto lower = [&box](int dir) {
          ::amrex::IntVect small = box.smallEnd();
          return small.shift(dir, -1);
        };
        auto upper = [&box](int dir) {
          ::amrex::IntVect big = box.bigEnd();
          return big.shift(dir, +1);
        };
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          if (ngrow[dir] && !gdomain.contains(lower(dir))) {
            function_(patch, Location{Direction(dir), 0}, Duration(time_point));
          }
          if (ngrow[dir] && !gdomain.contains(upper(dir))) {
            function_(patch, Location{Direction(dir), 1}, Duration(time_point));
          }
        }
      }
    }
  }

  Function function_;
  ::amrex::Geometry geom_;
  int level_num_;
};
} // namespace

void HyperbolicSplitIntegratorContext::FillGhostLayerTwoLevels(
    int fine, int coarse, Direction dir, BoundaryCondition boundary) {
  FUB_ASSERT(coarse >= 0 && fine > coarse);
  ::amrex::Vector<::amrex::BCRec> bcr(2 * AMREX_SPACEDIM); /* Fill it */
  ::amrex::MultiFab& scratch = GetScratch(fine, dir);
  const int nc = scratch.nComp();
  const ::amrex::Vector<::amrex::MultiFab*> cmf{&GetData(coarse)};
  const ::amrex::Vector<::amrex::MultiFab*> fmf{&GetData(fine)};
  const ::amrex::Vector<double> ct{GetTimePoint(coarse, dir).count()};
  const ::amrex::Vector<double> ft{GetTimePoint(fine, dir).count()};
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  const ::amrex::Geometry& fgeom = GetGeometry(fine);
  const ::amrex::IntVect ratio = 2 * ::amrex::IntVect::TheUnitVector();
  ::amrex::Interpolater* mapper = &::amrex::pc_interp;
  AdaptBoundaryCondition fine_condition(boundary, fgeom, fine);
  AdaptBoundaryCondition coarse_condition(boundary, cgeom, coarse);
  ::amrex::FillPatchTwoLevels(scratch, ft[0], cmf, ct, fmf, ft, 0, 0, nc, cgeom,
                              fgeom, coarse_condition, 0, fine_condition, 0,
                              ratio, mapper, bcr, 0);
}

void HyperbolicSplitIntegratorContext::FillGhostLayerSingleLevel(
    int level, Direction dir, BoundaryCondition boundary) {
  ::amrex::Vector<::amrex::BCRec> bcr(2 * AMREX_SPACEDIM); /* Fill it */
  ::amrex::MultiFab& scratch = GetScratch(level, dir);
  const int nc = scratch.nComp();
  const ::amrex::Vector<::amrex::MultiFab*> smf{&GetData(level)};
  const ::amrex::Vector<double> stime{GetTimePoint(level, dir).count()};
  const ::amrex::Geometry& geom = GetGeometry(level);
  AdaptBoundaryCondition condition(boundary, geom, level);
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
                                                                  Duration dt) {
  if (level > 0) {
    const ::amrex::MultiFab& fluxes = GetFluxes(level, dir);
  const double scale =
      1.0 / ipow(GetRatioToCoarserLevel(level), AMREX_SPACEDIM);
    data_[level].coarse_fine.FineAdd(fluxes, int(dir), 0, 0, fluxes.nComp(),
                                     scale);
  }
}

void HyperbolicSplitIntegratorContext::ResetCoarseFineFluxes(int fine,
                                                             int coarse,
                                                             Direction dir) {
  data_[fine].coarse_fine.ClearInternalBorders(GetGeometry(coarse));
  const ::amrex::MultiFab& flux = GetFluxes(coarse, dir);
  const ::amrex::BoxArray& boxes = flux.boxArray();
  const ::amrex::DistributionMapping& distribution = flux.DistributionMap();
  const int ncomp = flux.nComp();
  ::amrex::MultiFab zero(boxes, distribution, ncomp, 0);
  zero.setVal(0.0);
  data_[fine].coarse_fine.CrseInit(zero, int(dir), 0, 0, ncomp);
}

void HyperbolicSplitIntegratorContext::ApplyFluxCorrection(int fine, int coarse,
                                                           Duration dt,
                                                           Direction dir) {
  const int ncomp = GetPatchHierarchy()->GetDataDescription().n_cons_components;
  const ::amrex::Geometry& cgeom = GetGeometry(coarse);
  std::array<::amrex::MultiFab*, AMREX_SPACEDIM> crse_fluxes{AMREX_D_DECL(
      &GetFluxes(coarse, Direction::X), &GetFluxes(coarse, Direction::Y),
      &GetFluxes(coarse, Direction::Z))};
  data_[fine].coarse_fine.OverwriteFlux(crse_fluxes, 1.0, 0, 0, ncomp, cgeom);
}

CartesianCoordinates HyperbolicSplitIntegratorContext::GetCartesianCoordinates(
    PatchHandle patch) const {
  const ::amrex::Geometry& geom = GetGeometry(patch.level);
  const ::amrex::Box& box = GetPatchHierarchy()
                                ->GetPatchLevel(patch.level)
                                .data[patch.global_index]
                                .box();
  return ::fub::amrex::GetCartesianCoordinates(geom, box);
}

MPI_Comm HyperbolicSplitIntegratorContext::GetMpiCommunicator() const noexcept {
  return ::amrex::ParallelContext::CommunicatorAll();
}

void HyperbolicSplitIntegratorContext::PreAdvanceLevel(int level_num,
                                                       Direction dir,
                                                       Duration dt,
                                                       int subcycle) {
  const int d = int(dir);
  if (subcycle == 0 && level_num > 0 &&
      data_[level_num].regrid_time_point[d] != data_[level_num].time_point[d]) {
    gridding_->RegridAllFinerlevels(level_num - 1);
    for (int lvl = level_num; lvl < data_.size(); ++lvl) {
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

} // namespace amrex
} // namespace fub