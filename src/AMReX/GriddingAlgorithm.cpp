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

#include "fub/grid/AMReX/GriddingAlgorithm.hpp"
#include "fub/grid/AMReX/utility.hpp"

#include <AMReX_FillPatchUtil.H>

namespace fub::amrex {
GriddingAlgorithm::GriddingAlgorithm(std::shared_ptr<PatchHierarchy> hier,
                                     InitialData initial_data, Tagging tagging)
    : AmrCore(
          &hier->GetGridGeometry().coordinates,
          hier->GetOptions().max_number_of_levels - 1,
          ::amrex::Vector<int>(hier->GetGridGeometry().cell_dimensions.begin(),
                               hier->GetGridGeometry().cell_dimensions.end())),
      hierarchy_{std::move(hier)},
      initial_data_{std::move(initial_data)}, tagging_{std::move(tagging)} {
  n_error_buf =
      ::amrex::Vector<int>(hierarchy_->GetOptions().max_number_of_levels, 2);
  n_proper = 2;
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
    ::amrex::average_down(hierarchy_->GetPatchLevel(level).data,
                          hierarchy_->GetPatchLevel(level - 1).data,
                          hierarchy_->GetGeometry(level),
                          hierarchy_->GetGeometry(level - 1), first, size, 2);
  }
}

void GriddingAlgorithm::FillMultiFabFromLevel(::amrex::MultiFab& multifab,
                                              int level_number) {
  PatchLevel& level = hierarchy_->GetPatchLevel(level_number);
  const int n_comps = level.data.nComp();
  // TODO decide for BoundaryCondition interface
  ::amrex::Vector<::amrex::BCRec> bcr(n_comps);
  if (level_number == 0) {
    const ::amrex::Geometry& geom = hierarchy_->GetGeometry(level_number);
    auto no_condition = MakePhysBCFunct(geom, bcr, [](auto&&...) {});
    const ::amrex::Vector<::amrex::MultiFab*> smf{&level.data};
    const ::amrex::Vector<double> stime{level.time_point.count()};
    ::amrex::FillPatchSingleLevel(multifab, level.time_point.count(), smf,
                                  stime, 0, 0, n_comps, geom, no_condition, 0);
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
    auto fine_boundary = MakePhysBCFunct(fgeom, bcr, [](auto&&...) {});
    auto coarse_boundary = MakePhysBCFunct(cgeom, bcr, [](auto&&...) {});
    ::amrex::FillPatchTwoLevels(multifab, level.time_point.count(),
                                cmf, ct, fmf, ft, 0, 0, n_comps, cgeom, fgeom,
                                coarse_boundary, 0, fine_boundary, 0, ratio,
                                mapper, bcr, 0);
  }
}

void GriddingAlgorithm::ErrorEst(int level, ::amrex::TagBoxArray& tags,
                                 double, int /* ngrow */) {
  PatchLevel& old_level = hierarchy_->GetPatchLevel(level);
  ::amrex::MultiFab& data = old_level.data;
  const int ngrow = tags.nGrow();
  const int ncomp = data.nComp();
  ::amrex::MultiFab scratch(data.boxArray(), data.DistributionMap(), ncomp,
                            ngrow);
  FillMultiFabFromLevel(scratch, level);
  ForEachFAB(
      [&](::amrex::TagBox& tag_box, const ::amrex::FArrayBox& fab) {
        mdspan<const double, Rank + 1> data = MakeMdSpan(fab);
        mdspan<char, Rank> tags = MakeMdSpan(tag_box, 0);
        const ::amrex::Box box = fab.box();
        const ::amrex::Geometry& geom = hierarchy_->GetGeometry(level);
        CartesianCoordinates coords = GetCartesianCoordinates(geom, box);
        tagging_.TagCellsForRefinement(tags, data, coords);
      },
      tags, scratch);
}

void GriddingAlgorithm::MakeNewLevelFromScratch(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  const int n_comps = hierarchy_->GetDataDescription().n_state_components;
  hierarchy_->GetPatchLevel(level) = PatchLevel(
      level, Duration(time_point), box_array, distribution_mapping, n_comps);
  const ::amrex::Geometry& geom = hierarchy_->GetGeometry(level);
  FUB_ASSERT(::amrex::BoxArray(geom.Domain()).contains(box_array));
  ForEachFAB(
      [&](::amrex::FArrayBox& fab) {
        mdspan<double, Rank + 1> data = MakeMdSpan(fab);
        CartesianCoordinates coords = GetCartesianCoordinates(geom, fab.box());
        initial_data_.InitializeData(data, coords);
      },
      hierarchy_->GetPatchLevel(level).data);
}

void GriddingAlgorithm::MakeNewLevelFromCoarse(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  FUB_ASSERT(level > 0);
  const PatchLevel& coarse_level = hierarchy_->GetPatchLevel(level - 1);
  const int n_comps = hierarchy_->GetDataDescription().n_state_components;
  PatchLevel fine_level(level, Duration(time_point), box_array,
                        distribution_mapping, n_comps);
  const int cons_start = hierarchy_->GetDataDescription().first_cons_component;
  const int n_cons_components =
      hierarchy_->GetDataDescription().n_cons_components;
  ::amrex::Vector<::amrex::BCRec> bcr(n_comps);
  auto fine_boundary =
      MakePhysBCFunct(hierarchy_->GetGeometry(level), bcr, [](auto&&...) {});
  auto coarse_boundary = MakePhysBCFunct(hierarchy_->GetGeometry(level - 1),
                                         bcr, [](auto&&...) {});
  ::amrex::InterpFromCoarseLevel(
      fine_level.data, time_point, coarse_level.data, cons_start, cons_start,
      n_cons_components, hierarchy_->GetGeometry(level - 1),
      hierarchy_->GetGeometry(level), fine_boundary, 0, coarse_boundary, 0,
      {2, 2}, &::amrex::pc_interp, bcr, 0);
  hierarchy_->GetPatchLevel(level) = std::move(fine_level);
}

void GriddingAlgorithm::RemakeLevel(
    int level_number, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  const int n_comps = hierarchy_->GetDataDescription().n_state_components;
  PatchLevel new_level(level_number, Duration(time_point), box_array,
                       distribution_mapping, n_comps);
  FillMultiFabFromLevel(new_level.data, level_number);
  hierarchy_->GetPatchLevel(level_number) = std::move(new_level);
}

void GriddingAlgorithm::ClearLevel(int level) {
  hierarchy_->GetPatchLevel(level) = PatchLevel{};
}

} // namespace fub::amrex