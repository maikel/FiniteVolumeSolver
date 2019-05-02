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
#include <AMReX_ParmParse.H>

namespace fub::amrex {

int PrepareParmParseAndReturnNumberOfRefinementLevels(
    const PatchHierarchy& hier) {
  ::amrex::ParmParse pp("amr");
  const int dim = hier.GetDataDescription().dimension;
  FUB_ASSERT(dim >= 1);
  pp.add("blocking_factor_x", 8);
  pp.add("blocking_factor_y", dim >= 2 ? 8 : 1);
  pp.add("blocking_factor_z", dim >= 3 ? 8 : 1);
  return hier.GetMaxNumberOfLevels() - 1;
}

GriddingAlgorithm::GriddingAlgorithm(const GriddingAlgorithm& other)
    : AmrCore(
          &other.hierarchy_.GetGridGeometry().coordinates,
          PrepareParmParseAndReturnNumberOfRefinementLevels(other.hierarchy_),
          ::amrex::Vector<int>(
              other.hierarchy_.GetGridGeometry().cell_dimensions.begin(),
              other.hierarchy_.GetGridGeometry().cell_dimensions.end()),
          -1,
          ::amrex::Vector<::amrex::IntVect>(
              static_cast<std::size_t>(other.hierarchy_.GetMaxNumberOfLevels()),
              other.hierarchy_.GetRatioToCoarserLevel(
                  other.hierarchy_.GetOptions().max_number_of_levels - 1))),
      hierarchy_{other.hierarchy_},
      initial_data_{other.initial_data_}, tagging_{other.tagging_},
      boundary_condition_(other.boundary_condition_) {
  AmrMesh::verbose = other.verbose;
  AmrMesh::max_level = other.max_level;
  AmrMesh::ref_ratio = other.ref_ratio;
  AmrMesh::finest_level = other.finest_level;
  AmrMesh::n_error_buf = other.n_error_buf;
  AmrMesh::blocking_factor = other.blocking_factor;
  AmrMesh::max_grid_size = other.max_grid_size;
  AmrMesh::grid_eff = other.grid_eff;
  AmrMesh::n_proper = other.n_proper;
  AmrMesh::use_fixed_coarse_grids = other.use_fixed_coarse_grids;
  AmrMesh::use_fixed_upto_level = other.use_fixed_upto_level;
  AmrMesh::refine_grid_layout = other.refine_grid_layout;
  AmrMesh::check_input = other.check_input;
  AmrMesh::iterate_on_new_grids = other.iterate_on_new_grids;
  AmrMesh::use_new_chop = other.use_new_chop;
  const std::size_t size = static_cast<std::size_t>(AmrMesh::max_level + 1);
  AmrMesh::geom.resize(size);
  AmrMesh::dmap.resize(size);
  AmrMesh::grids.resize(size);
  for (int i = 0; i <= AmrMesh::finest_level; ++i) {
    const std::size_t ii = static_cast<std::size_t>(i);
    AmrMesh::geom[ii] = hierarchy_.GetGeometry(i);
    AmrMesh::dmap[ii] = hierarchy_.GetPatchLevel(i).distribution_mapping;
    AmrMesh::grids[ii] = hierarchy_.GetPatchLevel(i).box_array;
  }
}

GriddingAlgorithm& GriddingAlgorithm::
operator=(const GriddingAlgorithm& other) {
  GriddingAlgorithm tmp{other};
  return *this = std::move(tmp);
}

GriddingAlgorithm::GriddingAlgorithm(GriddingAlgorithm&& other) noexcept
    : AmrCore(
          &other.hierarchy_.GetGridGeometry().coordinates,
          PrepareParmParseAndReturnNumberOfRefinementLevels(other.hierarchy_),
          ::amrex::Vector<int>(
              other.hierarchy_.GetGridGeometry().cell_dimensions.begin(),
              other.hierarchy_.GetGridGeometry().cell_dimensions.end()),
          -1,
          ::amrex::Vector<::amrex::IntVect>(
              static_cast<std::size_t>(other.hierarchy_.GetMaxNumberOfLevels()),
              other.hierarchy_.GetRatioToCoarserLevel(
                  other.hierarchy_.GetOptions().max_number_of_levels - 1))),
      hierarchy_{std::move(other.hierarchy_)},
      initial_data_{std::move(other.initial_data_)}, tagging_{std::move(
                                                         other.tagging_)},
      boundary_condition_(std::move(other.boundary_condition_)) {
  AmrMesh::verbose = std::move(other.verbose);
  AmrMesh::max_level = std::move(other.max_level);
  AmrMesh::ref_ratio = std::move(other.ref_ratio);
  AmrMesh::finest_level = std::move(other.finest_level);
  AmrMesh::n_error_buf = std::move(other.n_error_buf);
  AmrMesh::blocking_factor = std::move(other.blocking_factor);
  AmrMesh::max_grid_size = std::move(other.max_grid_size);
  AmrMesh::grid_eff = std::move(other.grid_eff);
  AmrMesh::n_proper = std::move(other.n_proper);
  AmrMesh::use_fixed_coarse_grids = std::move(other.use_fixed_coarse_grids);
  AmrMesh::use_fixed_upto_level = std::move(other.use_fixed_upto_level);
  AmrMesh::refine_grid_layout = std::move(other.refine_grid_layout);
  AmrMesh::check_input = std::move(other.check_input);
  AmrMesh::iterate_on_new_grids = std::move(other.iterate_on_new_grids);
  AmrMesh::use_new_chop = std::move(other.use_new_chop);
  AmrMesh::geom = std::move(other.geom);
  AmrMesh::dmap = std::move(other.dmap);
  AmrMesh::grids = std::move(other.grids);
}

GriddingAlgorithm& GriddingAlgorithm::
operator=(GriddingAlgorithm&& other) noexcept {
  AmrMesh::verbose = std::move(other.verbose);
  AmrMesh::max_level = std::move(other.max_level);
  AmrMesh::ref_ratio = std::move(other.ref_ratio);
  AmrMesh::finest_level = std::move(other.finest_level);
  AmrMesh::n_error_buf = std::move(other.n_error_buf);
  AmrMesh::blocking_factor = std::move(other.blocking_factor);
  AmrMesh::max_grid_size = std::move(other.max_grid_size);
  AmrMesh::grid_eff = std::move(other.grid_eff);
  AmrMesh::n_proper = std::move(other.n_proper);
  AmrMesh::use_fixed_coarse_grids = std::move(other.use_fixed_coarse_grids);
  AmrMesh::use_fixed_upto_level = std::move(other.use_fixed_upto_level);
  AmrMesh::refine_grid_layout = std::move(other.refine_grid_layout);
  AmrMesh::check_input = std::move(other.check_input);
  AmrMesh::iterate_on_new_grids = std::move(other.iterate_on_new_grids);
  AmrMesh::use_new_chop = std::move(other.use_new_chop);
  AmrMesh::geom = std::move(other.geom);
  AmrMesh::dmap = std::move(other.dmap);
  AmrMesh::grids = std::move(other.grids);
  hierarchy_ = std::move(other.hierarchy_);
  initial_data_ = std::move(other.initial_data_);
  tagging_ = std::move(other.tagging_);
  boundary_condition_ = std::move(other.boundary_condition_);
  return *this;
}

GriddingAlgorithm::GriddingAlgorithm(PatchHierarchy hier,
                                     InitialData initial_data, Tagging tagging)
    : AmrCore(
          &hier.GetGridGeometry().coordinates,
          PrepareParmParseAndReturnNumberOfRefinementLevels(hier),
          ::amrex::Vector<int>(hier.GetGridGeometry().cell_dimensions.begin(),
                               hier.GetGridGeometry().cell_dimensions.end()),
          -1,
          ::amrex::Vector<::amrex::IntVect>(
              static_cast<std::size_t>(hier.GetMaxNumberOfLevels()),
              hier.GetRatioToCoarserLevel(
                  hier.GetOptions().max_number_of_levels - 1))),
      hierarchy_{std::move(hier)},
      initial_data_{std::move(initial_data)}, tagging_{std::move(tagging)} {
}

GriddingAlgorithm::GriddingAlgorithm(PatchHierarchy hier,
                                     InitialData initial_data, Tagging tagging,
                                     BoundaryCondition boundary)
    : AmrCore(
          &hier.GetGridGeometry().coordinates,
          PrepareParmParseAndReturnNumberOfRefinementLevels(hier),
          ::amrex::Vector<int>(hier.GetGridGeometry().cell_dimensions.begin(),
                               hier.GetGridGeometry().cell_dimensions.end()),
          -1,
          ::amrex::Vector<::amrex::IntVect>(
              static_cast<std::size_t>(hier.GetMaxNumberOfLevels()),
              hier.GetRatioToCoarserLevel(
                  hier.GetOptions().max_number_of_levels - 1))),
      hierarchy_{std::move(hier)}, initial_data_{std::move(initial_data)},
      tagging_{std::move(tagging)}, boundary_condition_{std::move(boundary)} {
}

bool GriddingAlgorithm::RegridAllFinerlevels(int which_level) {
  if (which_level < max_level) {
    const int before = AmrMesh::finest_level;
    AmrCore::regrid(which_level,
                    hierarchy_.GetPatchLevel(which_level).time_point.count());
    const int after = AmrMesh::finest_level;
    return before != after;
  }
  return false;
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
  PatchLevel& level = hierarchy_.GetPatchLevel(level_number);
  const int n_comps = level.data.nComp();
  // TODO decide for BoundaryCondition interface
  ::amrex::Vector<::amrex::BCRec> bcr(static_cast<std::size_t>(n_comps));
  if (level_number == 0) {
    const ::amrex::Geometry& geom = hierarchy_.GetGeometry(level_number);
    ::fub::amrex::BoundaryCondition boundary(boundary_condition_, geom,
                                             level_number, GetPatchHierarchy());
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
    ::fub::amrex::BoundaryCondition fine_boundary(boundary_condition_, fgeom,
                                                  level_number, GetPatchHierarchy());
    ::fub::amrex::BoundaryCondition coarse_boundary(boundary_condition_, cgeom,
                                                    level_number - 1, GetPatchHierarchy());
    ::amrex::FillPatchTwoLevels(multifab, level.time_point.count(), cmf, ct,
                                fmf, ft, 0, 0, n_comps, cgeom, fgeom,
                                coarse_boundary, 0, fine_boundary, 0, ratio,
                                mapper, bcr, 0);
  }
}

void GriddingAlgorithm::ErrorEst(int level, ::amrex::TagBoxArray& tags, double,
                                 int /* ngrow */) {
  PatchLevel& old_level = hierarchy_.GetPatchLevel(level);
  ::amrex::MultiFab& data = old_level.data;
  ::amrex::IntVect ngrow = ::amrex::IntVect::TheUnitVector();
  const int dim = hierarchy_.GetDataDescription().dimension;
  if (dim < AMREX_SPACEDIM) {
    std::fill_n(&ngrow[dim], AMREX_SPACEDIM - dim, 0);
  }
  const int ncomp = data.nComp();
  const ::amrex::BoxArray& ba = data.boxArray();
  const ::amrex::DistributionMapping& dm = data.DistributionMap();
  ::amrex::MultiFab scratch(ba, dm, ncomp, 2 * ngrow);
  FillMultiFabFromLevel(scratch, level);
  for (::amrex::MFIter mfi(ba, dm); mfi.isValid(); ++mfi) {
    PatchDataView<const double, Rank + 1> data =
        MakePatchDataView(scratch[mfi]);
    PatchDataView<char, Rank> tag = MakePatchDataView(tags[mfi], 0);
    PatchHandle patch{level, &mfi};
    tagging_.TagCellsForRefinement(tag, data, hierarchy_, patch);
  }
}

void GriddingAlgorithm::MakeNewLevelFromScratch(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  // Allocate level data.
  {
    const int n_comps = hierarchy_.GetDataDescription().n_state_components;
    hierarchy_.GetPatchLevel(level) = PatchLevel(
        level, Duration(time_point), box_array, distribution_mapping, n_comps);
  }
  // Initialize Data with stored intial data condition.
  for (::amrex::MFIter mfi(box_array, distribution_mapping); mfi.isValid();
       ++mfi) {
    PatchHandle patch{level, &mfi};
    PatchDataView<double, Rank + 1> data =
        MakePatchDataView(hierarchy_.GetPatchLevel(level).data[mfi]);
    initial_data_.InitializeData(data, hierarchy_, patch);
  }
}

void GriddingAlgorithm::MakeNewLevelFromCoarse(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  FUB_ASSERT(level > 0);
  const PatchLevel& coarse_level = hierarchy_.GetPatchLevel(level - 1);
  const int n_comps = hierarchy_.GetDataDescription().n_state_components;
  PatchLevel fine_level(level, Duration(time_point), box_array,
                        distribution_mapping, n_comps);
  const int cons_start = hierarchy_.GetDataDescription().first_cons_component;
  const int n_cons_components =
      hierarchy_.GetDataDescription().n_cons_components;
  ::amrex::Vector<::amrex::BCRec> bcr(static_cast<std::size_t>(n_comps));
  ::fub::amrex::BoundaryCondition fine_boundary(
      boundary_condition_, hierarchy_.GetGeometry(level), level, GetPatchHierarchy());
  ::fub::amrex::BoundaryCondition coarse_boundary(
      boundary_condition_, hierarchy_.GetGeometry(level - 1), level - 1, GetPatchHierarchy());
  ::amrex::InterpFromCoarseLevel(
      fine_level.data, time_point, coarse_level.data, cons_start, cons_start,
      n_cons_components, hierarchy_.GetGeometry(level - 1),
      hierarchy_.GetGeometry(level), fine_boundary, 0, coarse_boundary, 0,
      {AMREX_D_DECL(2, 2, 2)}, &::amrex::pc_interp, bcr, 0);
  hierarchy_.GetPatchLevel(level) = std::move(fine_level);
}

void GriddingAlgorithm::RemakeLevel(
    int level_number, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  const int n_comps = hierarchy_.GetDataDescription().n_state_components;
  PatchLevel new_level(level_number, Duration(time_point), box_array,
                       distribution_mapping, n_comps);
  FillMultiFabFromLevel(new_level.data, level_number);
  hierarchy_.GetPatchLevel(level_number) = std::move(new_level);
}

void GriddingAlgorithm::ClearLevel(int level) {
  hierarchy_.GetPatchLevel(level) = PatchLevel{};
}

const GriddingAlgorithm::BoundaryCondition&
GriddingAlgorithm::GetBoundaryCondition() const noexcept {
  return boundary_condition_;
}

const InitialData& GriddingAlgorithm::GetInitialCondition() const noexcept {
  return initial_data_;
}

const Tagging& GriddingAlgorithm::GetTagging() const noexcept {
  return tagging_;
}

} // namespace fub::amrex
