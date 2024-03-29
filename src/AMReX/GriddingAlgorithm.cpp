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

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/boundary_condition/BoundaryConditionRef.hpp"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_ParmParse.H>

namespace fub::amrex {

GriddingAlgorithm::GriddingAlgorithm(const GriddingAlgorithm& other)
    : AmrCore(other.hierarchy_.GetGeometry(0),
              static_cast<const ::amrex::AmrInfo&>(other)),
      hierarchy_{other.hierarchy_},
      initial_data_{other.initial_data_}, tagging_{other.tagging_},
      boundary_condition_(other.boundary_condition_) {
  AmrMesh::finest_level = other.finest_level;
  AmrMesh::geom = other.geom;
  AmrMesh::dmap = other.dmap;
  AmrMesh::grids = other.grids;
}

GriddingAlgorithm&
GriddingAlgorithm::operator=(const GriddingAlgorithm& other) {
  GriddingAlgorithm tmp{other};
  return *this = std::move(tmp);
}

GriddingAlgorithm::GriddingAlgorithm(PatchHierarchy hier,
                                     AnyInitialData initial_data,
                                     AnyTaggingMethod tagging,
                                     AnyBoundaryCondition bc)
    : AmrCore(
          &hier.GetGridGeometry().coordinates, hier.GetMaxNumberOfLevels() - 1,
          ::amrex::Vector<int>(hier.GetGridGeometry().cell_dimensions.begin(),
                               hier.GetGridGeometry().cell_dimensions.end()),
          -1,
          ::amrex::Vector<::amrex::IntVect>(
              static_cast<std::size_t>(hier.GetMaxNumberOfLevels()),
              hier.GetRatioToCoarserLevel(
                  hier.GetOptions().max_number_of_levels))),
      hierarchy_(std::move(hier)), initial_data_(std::move(initial_data)),
      tagging_(std::move(tagging)), boundary_condition_(std::move(bc)) {
  const PatchHierarchyOptions& options = hierarchy_.GetOptions();
  AmrMesh::SetMaxGridSize(options.max_grid_size);
  AmrMesh::SetBlockingFactor(options.blocking_factor);
  AmrMesh::SetNProper(options.n_proper);
  AmrMesh::SetGridEff(options.grid_efficiency);
  if (hier.GetMaxNumberOfLevels() == 1) {
    AmrMesh::ref_ratio.resize(1, ::amrex::IntVect(1));
  }
  AmrMesh::verbose = options.verbose;
  AmrCore::verbose = options.verbose;
  AmrMesh::n_error_buf = ::amrex::Vector<::amrex::IntVect>(
      static_cast<std::size_t>(AmrMesh::n_error_buf.size()),
      options.n_error_buf);
  if (hierarchy_.GetNumberOfLevels() > 0) {
    for (int i = 0; i < hierarchy_.GetNumberOfLevels(); ++i) {
      const std::size_t ii = static_cast<std::size_t>(i);
      AmrMesh::geom[ii] = hierarchy_.GetGeometry(i);
      AmrMesh::dmap[ii] = hierarchy_.GetPatchLevel(i).distribution_mapping;
      AmrMesh::grids[ii] = hierarchy_.GetPatchLevel(i).box_array;
    }
    AmrMesh::finest_level = hierarchy_.GetNumberOfLevels() - 1;
  }
}

int GriddingAlgorithm::RegridAllFinerlevels(int which_level) {
  if (which_level < max_level) {
    auto timer = hierarchy_.GetCounterRegistry()->get_timer(
        "GriddingAlgorithm::RegridAllFinerLevels");
    const ::amrex::Vector<::amrex::BoxArray> before =
        ::amrex::AmrMesh::boxArray();
    AmrCore::regrid(which_level,
                    hierarchy_.GetPatchLevel(which_level).time_point.count());
    const ::amrex::Vector<::amrex::BoxArray> after =
        ::amrex::AmrMesh::boxArray();
    const int bsize = static_cast<int>(before.size());
    const int asize = static_cast<int>(after.size());
    int i = which_level + 1;
    for (; i < std::min(asize, bsize); ++i) {
      const std::size_t ii = static_cast<std::size_t>(i);
      if (before[ii] != after[ii]) {
        return i;
      }
    }
    if (i < bsize) {
      return i;
    }
  }
  return max_level + 1;
}

void GriddingAlgorithm::InitializeHierarchy(double level_time) {
  ::amrex::AmrCore::MakeNewGrids(level_time);
  const int n_levels = hierarchy_.GetNumberOfLevels();
  const int first = hierarchy_.GetDataDescription().first_cons_component;
  const int size = hierarchy_.GetDataDescription().n_cons_components;
  for (int level = n_levels - 1; level > 0; --level) {
    FUB_ASSERT(level > 0);
    ::amrex::average_down(hierarchy_.GetPatchLevel(level).data,
                          hierarchy_.GetPatchLevel(level - 1).data,
                          hierarchy_.GetGeometry(level),
                          hierarchy_.GetGeometry(level - 1), first, size,
                          hierarchy_.GetRatioToCoarserLevel(level));
    FUB_ASSERT(!hierarchy_.GetPatchLevel(level).data.contains_nan());
  }
  FUB_ASSERT(!hierarchy_.GetPatchLevel(0).data.contains_nan());
}

void GriddingAlgorithm::FillMultiFabFromLevel(::amrex::MultiFab& multifab,
                                              int level_number) {
  FillMultiFabFromLevel(multifab, level_number, boundary_condition_);
}

void GriddingAlgorithm::FillMultiFabFromLevel(::amrex::MultiFab& multifab,
                                              int level_number,
                                              AnyBoundaryCondition& bc) {
  PatchLevel& level = hierarchy_.GetPatchLevel(level_number);
  const int n_comps = level.data.nComp();
  ::amrex::Vector<::amrex::BCRec> bcr(static_cast<std::size_t>(n_comps));
  if (level_number == 0) {
    const ::amrex::Geometry& geom = hierarchy_.GetGeometry(level_number);
    BoundaryConditionRef boundary(bc, *this, level_number);
    const ::amrex::Vector<::amrex::MultiFab*> smf{&level.data};
    const ::amrex::Vector<double> stime{level.time_point.count()};
    ::amrex::FillPatchSingleLevel(multifab, level.time_point.count(), smf,
                                  stime, 0, 0, n_comps, geom, boundary, 0);
  } else {
    PatchLevel& coarse_level = hierarchy_.GetPatchLevel(level_number - 1);
    const ::amrex::Vector<::amrex::MultiFab*> cmf{&coarse_level.data};
    const ::amrex::Vector<::amrex::MultiFab*> fmf{&level.data};
    const ::amrex::Vector<double> ct{coarse_level.time_point.count()};
    const ::amrex::Vector<double> ft{level.time_point.count()};
    const ::amrex::Geometry& cgeom = hierarchy_.GetGeometry(level_number - 1);
    const ::amrex::Geometry& fgeom = hierarchy_.GetGeometry(level_number);
    const ::amrex::IntVect ratio =
        hierarchy_.GetRatioToCoarserLevel(level_number);
    ::amrex::Interpolater* mapper = &::amrex::pc_interp;
    const int fine = level_number;
    const int coarse = level_number - 1;
    BoundaryConditionRef fine_boundary(bc, *this, fine);
    BoundaryConditionRef coarse_boundary(bc, *this, coarse);
#ifdef AMREX_USE_EB
    span<const ::amrex::EB2::IndexSpace*> index_space =
        hierarchy_.GetIndexSpaces();
    ::amrex::FillPatchTwoLevels(
        multifab, level.time_point.count(), *index_space[fine], cmf, ct, fmf,
        ft, 0, 0, n_comps, cgeom, fgeom, coarse_boundary, 0, fine_boundary, 0,
        ratio, mapper, bcr, 0, ::amrex::NullInterpHook<::amrex::FArrayBox>(),
        ::amrex::NullInterpHook<::amrex::FArrayBox>());
#else
    ::amrex::FillPatchTwoLevels(multifab, level.time_point.count(), cmf, ct,
                                fmf, ft, 0, 0, n_comps, cgeom, fgeom,
                                coarse_boundary, 0, fine_boundary, 0, ratio,
                                mapper, bcr, 0);
#endif
  }
}

void GriddingAlgorithm::ErrorEst(int level, ::amrex::TagBoxArray& tags,
                                 double tp, int /* ngrow */) {
  auto timer =
      hierarchy_.GetCounterRegistry()->get_timer("GriddingAlgorithm::ErrorEst");
  const Duration time_point(tp);
  tagging_.TagCellsForRefinement(tags, *this, level, time_point);
}

void GriddingAlgorithm::MakeNewLevelFromScratch(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  auto timer = hierarchy_.GetCounterRegistry()->get_timer(
      "GriddingAlgorithm::MakeNewLevelFromScratch");
  // Allocate level data.
  const DataDescription& desc = hierarchy_.GetDataDescription();
  if (hierarchy_.GetNumberOfLevels() == level) {
    hierarchy_.PushBack(PatchLevel(level, Duration(time_point), box_array,
                                   distribution_mapping, desc,
                                   hierarchy_.GetMFInfo()));
  } else {
    PatchLevel patch_level(level, Duration(time_point), box_array,
                           distribution_mapping, desc, hierarchy_.GetMFInfo());
    hierarchy_.GetPatchLevel(level) = std::move(patch_level);
  }
  PatchLevel& patch_level = hierarchy_.GetPatchLevel(level);
  initial_data_.InitializeData(patch_level, *this, level,
                               static_cast<Duration>(time_point));
}

void GriddingAlgorithm::MakeNewLevelFromCoarse(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  auto timer = hierarchy_.GetCounterRegistry()->get_timer(
      "GriddingAlgorithm::MakeNewLevelFromCoarse");
  FUB_ASSERT(0 < level && level <= hierarchy_.GetNumberOfLevels());
  const PatchLevel& coarse_level = hierarchy_.GetPatchLevel(level - 1);
  const int n_comps = hierarchy_.GetDataDescription().n_state_components;
  PatchLevel fine_level(level, Duration(time_point), box_array,
                        distribution_mapping, hierarchy_.GetDataDescription(),
                        hierarchy_.GetMFInfo());
  const int cons_start = hierarchy_.GetDataDescription().first_cons_component;
  const int n_cons_components =
      hierarchy_.GetDataDescription().n_cons_components;
  ::amrex::Vector<::amrex::BCRec> bcr(static_cast<std::size_t>(n_comps));
  BoundaryConditionRef fine_boundary(boundary_condition_, *this, level);
  BoundaryConditionRef coarse_boundary(boundary_condition_, *this, level - 1);
  ::amrex::InterpFromCoarseLevel(
      fine_level.data, time_point, coarse_level.data, cons_start, cons_start,
      n_cons_components, hierarchy_.GetGeometry(level - 1),
      hierarchy_.GetGeometry(level), fine_boundary, 0, coarse_boundary, 0,
      {AMREX_D_DECL(2, 2, 2)}, &::amrex::pc_interp, bcr, 0);
  if (level == hierarchy_.GetNumberOfLevels()) {
    hierarchy_.PushBack(std::move(fine_level));
  } else {
    fine_level.data.ParallelCopy(hierarchy_.GetPatchLevel(level).data);
    fine_level.cycles = hierarchy_.GetPatchLevel(level).cycles;
    fine_level.time_point = hierarchy_.GetPatchLevel(level).time_point;
    hierarchy_.GetPatchLevel(level) = std::move(fine_level);
  }
}

void GriddingAlgorithm::RemakeLevel(
    int level_number, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  auto timer = hierarchy_.GetCounterRegistry()->get_timer(
      "GriddingAlgorithm::RemakeLevel");
  const DataDescription desc = hierarchy_.GetDataDescription();
  PatchLevel new_level(level_number, Duration(time_point), box_array,
                       distribution_mapping, desc, hierarchy_.GetMFInfo());
  FillMultiFabFromLevel(new_level.data, level_number);
  hierarchy_.GetPatchLevel(level_number) = std::move(new_level);
}

void GriddingAlgorithm::ClearLevel([[maybe_unused]] int level) {
  FUB_ASSERT(0 < level && level + 1 == hierarchy_.GetNumberOfLevels());
  hierarchy_.PopBack();
}

const AnyBoundaryCondition&
GriddingAlgorithm::GetBoundaryCondition() const noexcept {
  return boundary_condition_;
}

AnyBoundaryCondition& GriddingAlgorithm::GetBoundaryCondition() noexcept {
  return boundary_condition_;
}

const AnyInitialData& GriddingAlgorithm::GetInitialCondition() const noexcept {
  return initial_data_;
}

const AnyTaggingMethod& GriddingAlgorithm::GetTagging() const noexcept {
  return tagging_;
}

} // namespace fub::amrex
