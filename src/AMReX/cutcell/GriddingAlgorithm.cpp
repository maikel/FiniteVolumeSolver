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

#include "fub/grid/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/grid/AMReX/utility.hpp"

#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

namespace fub::amrex::cutcell {

GriddingAlgorithm::GriddingAlgorithm(std::shared_ptr<PatchHierarchy> hier,
                                     InitialData data, Tagging tagging,
                                     BoundaryCondition boundary)
    : AmrCore(
          &hier->GetGridGeometry().coordinates,
          hier->GetOptions().max_number_of_levels - 1,
          ::amrex::Vector<int>(hier->GetGridGeometry().cell_dimensions.begin(),
                               hier->GetGridGeometry().cell_dimensions.end())),
      hierarchy_{std::move(hier)}, initial_condition_{std::move(data)},
      boundary_condition_{std::move(boundary)}, tagging_{std::move(tagging)} {
  n_error_buf = ::amrex::Vector<int>(
      static_cast<std::size_t>(hierarchy_->GetOptions().max_number_of_levels),
      2);
  n_proper = 2;
}

const std::shared_ptr<PatchHierarchy>&
GriddingAlgorithm::GetPatchHierarchy() const noexcept {
  return hierarchy_;
}

bool GriddingAlgorithm::RegridAllFinerlevels(int which_level) {
  if (which_level < max_level) {
    const int before = AmrMesh::finest_level;
    AmrCore::regrid(which_level,
                    hierarchy_->GetPatchLevel(which_level).time_point.count());
    const int after = AmrMesh::finest_level;
    return before != after;
  }
  return false;
}

void GriddingAlgorithm::InitializeHierarchy(double level_time) {
  ::amrex::AmrCore::MakeNewGrids(level_time);
  const int n_levels = hierarchy_->GetNumberOfLevels();
  const int first = hierarchy_->GetDataDescription().first_cons_component;
  const int size = hierarchy_->GetDataDescription().n_cons_components;
  for (int level = n_levels - 1; level > 0; --level) {
    FUB_ASSERT(level > 0);
    ::amrex::EB_average_down(hierarchy_->GetPatchLevel(level).data,
                             hierarchy_->GetPatchLevel(level - 1).data, first,
                             size, 2);
  }
}

void GriddingAlgorithm::FillMultiFabFromLevel(::amrex::MultiFab& multifab,
                                              int level_number) {
  PatchLevel& level = hierarchy_->GetPatchLevel(level_number);
  const int n_comps = level.data.nComp();
  // TODO decide for BoundaryCondition interface
  ::amrex::Vector<::amrex::BCRec> bcr(static_cast<std::size_t>(n_comps));
  if (level_number == 0) {
    const ::amrex::Geometry& geom = hierarchy_->GetGeometry(level_number);
    const ::amrex::Vector<::amrex::MultiFab*> smf{&level.data};
    const ::amrex::Vector<double> stime{level.time_point.count()};
    ::fub::amrex::BoundaryCondition boundary(boundary_condition_, geom,
                                             level_number);
    ::amrex::FillPatchSingleLevel(multifab, level.time_point.count(), smf,
                                  stime, 0, 0, n_comps, geom, boundary, 0);
  } else {
    PatchLevel& coarse_level = hierarchy_->GetPatchLevel(level_number - 1);
    const ::amrex::Vector<::amrex::MultiFab*> cmf{&coarse_level.data};
    const ::amrex::Vector<::amrex::MultiFab*> fmf{&level.data};
    const ::amrex::Vector<double> ct{coarse_level.time_point.count()};
    const ::amrex::Vector<double> ft{level.time_point.count()};
    const ::amrex::Geometry& cgeom = hierarchy_->GetGeometry(level_number - 1);
    const ::amrex::Geometry& fgeom = hierarchy_->GetGeometry(level_number);
    const ::amrex::IntVect ratio = 2 * ::amrex::IntVect::TheUnitVector();
    ::amrex::Interpolater* mapper = &::amrex::pc_interp;
    ::fub::amrex::BoundaryCondition fine_boundary(boundary_condition_, fgeom,
                                                  level_number);
    ::fub::amrex::BoundaryCondition coarse_boundary(boundary_condition_, cgeom,
                                                    level_number - 1);
    ::amrex::FillPatchTwoLevels(multifab, level.time_point.count(), cmf, ct,
                                fmf, ft, 0, 0, n_comps, cgeom, fgeom,
                                coarse_boundary, 0, fine_boundary, 0, ratio,
                                mapper, bcr, 0);
  }
}

void GriddingAlgorithm::ErrorEst(int level, ::amrex::TagBoxArray& tags, double,
                                 int /* ngrow */) {
  PatchLevel& old_level = hierarchy_->GetPatchLevel(level);
  ::amrex::MultiFab& data = old_level.data;
  const int ncomp = data.nComp();
  const ::amrex::BoxArray& ba = data.boxArray();
  const ::amrex::DistributionMapping& dm = data.DistributionMap();
  ::amrex::MultiFab scratch(ba, dm, ncomp, 2);
  FillMultiFabFromLevel(scratch, level);

  for (::amrex::MFIter mfi(ba, dm); mfi.isValid(); ++mfi) {
    PatchDataView<const double, Rank + 1> data =
        MakePatchDataView(scratch[mfi]);
    PatchDataView<char, Rank> tag = MakePatchDataView(tags[mfi], 0);
    PatchHandle patch{level, &mfi};
    tagging_.TagCellsForRefinement(tag, data, patch);
  }
}

const GriddingAlgorithm::BoundaryCondition&
GriddingAlgorithm::GetBoundaryCondition() const noexcept {
  return boundary_condition_;
}

void GriddingAlgorithm::MakeNewLevelFromScratch(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  // Allocate level data.
  {
    const int n_comps = hierarchy_->GetDataDescription().n_state_components;
    std::shared_ptr<::amrex::EBFArrayBoxFactory> eb_factory =
        ::amrex::makeEBFabFactory(
            hierarchy_->GetOptions()
                .index_spaces[static_cast<std::size_t>(level)],
            hierarchy_->GetGeometry(level), box_array, distribution_mapping,
            {4, 4, 4}, ::amrex::EBSupport::full);
    hierarchy_->GetPatchLevel(level) =
        PatchLevel(level, Duration(time_point), box_array, distribution_mapping,
                   n_comps, std::move(eb_factory));
  }
  // Initialize Data with stored intial data condition.
  for (::amrex::MFIter mfi(box_array, distribution_mapping); mfi.isValid();
       ++mfi) {
    PatchHandle patch{level, &mfi};
    PatchDataView<double, Rank + 1> data =
        MakePatchDataView(hierarchy_->GetPatchLevel(level).data[mfi]);
    initial_condition_.InitializeData(data, patch);
  }
}

void GriddingAlgorithm::MakeNewLevelFromCoarse(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  FUB_ASSERT(level > 0);
  const PatchLevel& coarse_level = hierarchy_->GetPatchLevel(level - 1);
  const int n_comps = hierarchy_->GetDataDescription().n_state_components;
  const ::amrex::Geometry& geom = hierarchy_->GetGeometry(level);
  std::shared_ptr<::amrex::EBFArrayBoxFactory> factory =
      ::amrex::makeEBFabFactory(
          hierarchy_->GetOptions()
              .index_spaces[static_cast<std::size_t>(level)],
          geom, box_array, distribution_mapping, {4, 4, 4},
          ::amrex::EBSupport::full);
  PatchLevel fine_level(level, Duration(time_point), box_array,
                        distribution_mapping, n_comps, std::move(factory));
  const int cons_start = hierarchy_->GetDataDescription().first_cons_component;
  const int n_cons_components =
      hierarchy_->GetDataDescription().n_cons_components;
  ::amrex::Vector<::amrex::BCRec> bcr(static_cast<std::size_t>(n_comps));
  ::fub::amrex::BoundaryCondition fine_boundary(
      boundary_condition_, hierarchy_->GetGeometry(level), level);
  ::fub::amrex::BoundaryCondition coarse_boundary(
      boundary_condition_, hierarchy_->GetGeometry(level - 1), level - 1);

  ::amrex::InterpFromCoarseLevel(
      fine_level.data, time_point, coarse_level.data, cons_start, cons_start,
      n_cons_components, hierarchy_->GetGeometry(level - 1),
      hierarchy_->GetGeometry(level), fine_boundary, 0, coarse_boundary, 0,
      {AMREX_D_DECL(2, 2, 2)}, &::amrex::pc_interp, bcr, 0);
  hierarchy_->GetPatchLevel(level) = std::move(fine_level);
}

void GriddingAlgorithm::RemakeLevel(
    int level_number, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  const int n_comps = hierarchy_->GetDataDescription().n_state_components;
  const ::amrex::Geometry& geom = hierarchy_->GetGeometry(level_number);
  std::shared_ptr<::amrex::EBFArrayBoxFactory> factory =
      ::amrex::makeEBFabFactory(
          hierarchy_->GetOptions()
              .index_spaces[static_cast<std::size_t>(level_number)],
          geom, box_array, distribution_mapping, {4, 4, 4},
          ::amrex::EBSupport::full);
  PatchLevel new_level(level_number, Duration(time_point), box_array,
                       distribution_mapping, n_comps, std::move(factory));
  FillMultiFabFromLevel(new_level.data, level_number);
  hierarchy_->GetPatchLevel(level_number) = std::move(new_level);
}

void GriddingAlgorithm::ClearLevel(int level) {
  hierarchy_->GetPatchLevel(level) = PatchLevel{};
}

} // namespace fub::amrex::cutcell
