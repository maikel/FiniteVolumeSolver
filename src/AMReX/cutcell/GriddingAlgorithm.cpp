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

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"

#include "fub/AMReX/ForEachFab.hpp"

#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

namespace fub::amrex::cutcell {

GriddingAlgorithm::GriddingAlgorithm(const GriddingAlgorithm& other)
    : AmrCore(other.hierarchy_.GetGeometry(0),
              static_cast<const ::amrex::AmrInfo&>(other)),
      hierarchy_{other.hierarchy_},
      initial_condition_{other.initial_condition_}, tagging_{other.tagging_},
      boundary_condition_(other.boundary_condition_) {
  const std::size_t size = static_cast<std::size_t>(AmrInfo::max_level + 1);
  AmrMesh::geom.resize(size);
  AmrMesh::dmap.resize(size);
  AmrMesh::grids.resize(size);
  for (int i = 0; i <= AmrMesh::finest_level; ++i) {
    const std::size_t ii = static_cast<std::size_t>(i);
    AmrMesh::geom[ii] = hierarchy_.GetGeometry(i);
    AmrMesh::dmap[ii] = hierarchy_.GetPatchLevel(i).distribution_mapping;
    AmrMesh::grids[ii] = hierarchy_.GetPatchLevel(i).box_array;
  }
  for (int level = 0; level < hierarchy_.GetMaxNumberOfLevels(); ++level) {
    boundary_condition_[static_cast<std::size_t>(level)].geometry =
        hierarchy_.GetGeometry(level);
    boundary_condition_[static_cast<std::size_t>(level)].parent = this;
  }
}

GriddingAlgorithm&
GriddingAlgorithm::operator=(const GriddingAlgorithm& other) {
  GriddingAlgorithm tmp{other};
  return *this = std::move(tmp);
}

GriddingAlgorithm::GriddingAlgorithm(GriddingAlgorithm&& other) noexcept
    : AmrCore(std::move(other)), hierarchy_{std::move(other.hierarchy_)},
      initial_condition_{std::move(other.initial_condition_)},
      tagging_{std::move(other.tagging_)},
      boundary_condition_(std::move(other.boundary_condition_)) {
  for (int level = 0; level < hierarchy_.GetMaxNumberOfLevels(); ++level) {
    boundary_condition_[static_cast<std::size_t>(level)].geometry =
        hierarchy_.GetGeometry(level);
    boundary_condition_[static_cast<std::size_t>(level)].parent = this;
  }
}

GriddingAlgorithm&
GriddingAlgorithm::operator=(GriddingAlgorithm&& other) noexcept {
  hierarchy_ = std::move(other.hierarchy_);
  initial_condition_ = std::move(other.initial_condition_);
  tagging_ = std::move(other.tagging_);
  boundary_condition_ = std::move(other.boundary_condition_);
  *static_cast<::amrex::AmrCore*>(this) =
      static_cast<::amrex::AmrCore&&>(other);
  for (int level = 0; level < hierarchy_.GetMaxNumberOfLevels(); ++level) {
    boundary_condition_[static_cast<std::size_t>(level)].geometry =
        hierarchy_.GetGeometry(level);
    boundary_condition_[static_cast<std::size_t>(level)].parent = this;
  }
  return *this;
}

GriddingAlgorithm::GriddingAlgorithm(PatchHierarchy hier,
                                     AnyInitialData data,
                                     Tagging tagging,
                                     AnyBoundaryCondition boundary)
    : AmrCore(
          &hier.GetGridGeometry().coordinates,
          hier.GetOptions().max_number_of_levels - 1,
          ::amrex::Vector<int>(hier.GetGridGeometry().cell_dimensions.begin(),
                               hier.GetGridGeometry().cell_dimensions.end())),
      hierarchy_{std::move(hier)},
      initial_condition_{std::move(data)}, tagging_{std::move(tagging)},
      boundary_condition_(std::size_t(hier.GetMaxNumberOfLevels()),
                          std::move(boundary)) {
  const PatchHierarchyOptions& options = hierarchy_.GetOptions();
  AmrMesh::SetMaxGridSize(options.max_grid_size);
  AmrMesh::SetBlockingFactor(options.blocking_factor);
  AmrMesh::SetNProper(options.n_proper);
  AmrMesh::SetGridEff(options.grid_efficiency);
  AmrMesh::verbose = options.verbose;
  AmrCore::verbose = options.verbose;
  AmrMesh::n_error_buf = ::amrex::Vector<::amrex::IntVect>(
      AmrMesh::n_error_buf.size(), options.n_error_buf);
  if (hier.GetNumberOfLevels() > 0) {
    for (int i = 0; i < hier.GetNumberOfLevels(); ++i) {
      const std::size_t ii = static_cast<std::size_t>(i);
      AmrMesh::geom[ii] = hierarchy_.GetGeometry(i);
      AmrMesh::dmap[ii] = hierarchy_.GetPatchLevel(i).distribution_mapping;
      AmrMesh::grids[ii] = hierarchy_.GetPatchLevel(i).box_array;
    }
  }
  for (int level = 0; level < hierarchy_.GetMaxNumberOfLevels(); ++level) {
    boundary_condition_[static_cast<std::size_t>(level)].geometry =
        hierarchy_.GetGeometry(level);
    boundary_condition_[static_cast<std::size_t>(level)].parent = this;
  }
}

PatchHierarchy& GriddingAlgorithm::GetPatchHierarchy() noexcept {
  return hierarchy_;
}

const PatchHierarchy& GriddingAlgorithm::GetPatchHierarchy() const noexcept {
  return hierarchy_;
}

int GriddingAlgorithm::RegridAllFinerlevels(int which_level) {
  if (which_level < max_level) {
    auto timer = hierarchy_.GetCounterRegistry()->get_timer(
        "cutcell::GriddingAlgorithm::RegridAllFinerLevels");
    const ::amrex::Vector<::amrex::BoxArray> before =
        ::amrex::AmrMesh::boxArray();
    AmrCore::regrid(which_level,
                    hierarchy_.GetPatchLevel(which_level).time_point.count());
    const ::amrex::Vector<::amrex::BoxArray> after =
        ::amrex::AmrMesh::boxArray();
    FUB_ASSERT(before.size() == after.size());
    for (int i = which_level + 1; i < before.size(); ++i) {
      if (before[i] != after[i]) {
        return i;
      }
    }
  }
  return max_level + 1;
}

void GriddingAlgorithm::InitializeHierarchy(double level_time) {
  auto timer = hierarchy_.GetCounterRegistry()->get_timer(
      "cutcell::GriddingAlgorithm::InitializeHierarchy");
  ::amrex::AmrCore::MakeNewGrids(level_time);
  const int n_levels = hierarchy_.GetNumberOfLevels();
  const int first = hierarchy_.GetDataDescription().first_cons_component;
  const int size = hierarchy_.GetDataDescription().n_cons_components;
  for (int level = n_levels - 1; level > 0; --level) {
    FUB_ASSERT(level > 0);
    ::amrex::EB_average_down(hierarchy_.GetPatchLevel(level).data,
                             hierarchy_.GetPatchLevel(level - 1).data, first,
                             size, 2);
  }
}

void GriddingAlgorithm::FillMultiFabFromLevel(::amrex::MultiFab& multifab,
                                              int level_number) {
  PatchLevel& level = hierarchy_.GetPatchLevel(level_number);
  const int n_comps = level.data.nComp();
  // TODO decide for BoundaryCondition interface
  ::amrex::Vector<::amrex::BCRec> bcr(static_cast<std::size_t>(n_comps));
  if (level_number == 0) {
    const ::amrex::Geometry& geom = hierarchy_.GetGeometry(level_number);
    const ::amrex::Vector<::amrex::MultiFab*> smf{&level.data};
    const ::amrex::Vector<double> stime{level.time_point.count()};
    AnyBoundaryCondition& boundary = boundary_condition_[size_t(level_number)];
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
    const ::amrex::IntVect ratio = 2 * ::amrex::IntVect::TheUnitVector();
    ::amrex::Interpolater* mapper = &::amrex::pc_interp;
    const std::size_t fine = std::size_t(level_number);
    const std::size_t coarse = std::size_t(level_number - 1);
    AnyBoundaryCondition& fine_boundary = boundary_condition_[fine];
    AnyBoundaryCondition& coarse_boundary = boundary_condition_[coarse];
    ::amrex::FillPatchTwoLevels(
        multifab, level.time_point.count(),
        *hierarchy_.GetOptions().index_spaces[fine], cmf, ct, fmf, ft, 0, 0,
        n_comps, cgeom, fgeom, coarse_boundary, 0, fine_boundary, 0, ratio,
        mapper, bcr, 0, ::amrex::NullInterpHook<::amrex::FArrayBox>(),
        ::amrex::NullInterpHook<::amrex::FArrayBox>());
  }
}

void GriddingAlgorithm::ErrorEst(int level, ::amrex::TagBoxArray& tags,
                                 double tp, int /* ngrow */) {
  auto timer = hierarchy_.GetCounterRegistry()->get_timer(
      "cutcell::GriddingAlgorithm::ErrorEst");
  tagging_.TagCellsForRefinement(tags, Duration(tp), level, *this);
}

void GriddingAlgorithm::SetBoundaryCondition(int level,
                                             AnyBoundaryCondition&& condition) {
  boundary_condition_[std::size_t(level)] = std::move(condition);
}

void GriddingAlgorithm::SetBoundaryCondition(
    int level, const AnyBoundaryCondition& condition) {
  boundary_condition_[std::size_t(level)] = condition;
}

const AnyBoundaryCondition&
GriddingAlgorithm::GetBoundaryCondition(int level) const noexcept {
  return boundary_condition_[std::size_t(level)];
}

AnyBoundaryCondition&
GriddingAlgorithm::GetBoundaryCondition(int level) noexcept {
  return boundary_condition_[std::size_t(level)];
}

::amrex::DistributionMapping GriddingAlgorithm::LoadBalance(
    int level, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) const {
  auto timer = hierarchy_.GetCounterRegistry()->get_timer(
      "cutcell::GriddingAlgorithm::LoadBalance");
  std::unique_ptr<::amrex::EBFArrayBoxFactory> eb_factory =
      ::amrex::makeEBFabFactory(
          hierarchy_.GetOptions().index_spaces[static_cast<std::size_t>(level)],
          hierarchy_.GetGeometry(level), box_array, distribution_mapping,
          {0, 0, 0}, ::amrex::EBSupport::basic);
  const ::amrex::FabArray<::amrex::EBCellFlagFab>& flags =
      eb_factory->getMultiEBCellFlagFab();
  ::amrex::MultiFab weigths(box_array, distribution_mapping, 1, 0);
  weigths.setVal(1.0);
  for (::amrex::MFIter mfi(weigths); mfi.isValid(); ++mfi) {
    if (flags[mfi].getType() == ::amrex::FabType::singlevalued) {
      weigths[mfi].setVal(hierarchy_.GetOptions().cutcell_load_balance_weight);
    }
  }
  EB_set_covered(weigths, 0.001);
  return ::amrex::DistributionMapping::makeKnapSack(weigths);
}

void GriddingAlgorithm::MakeNewLevelFromScratch(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_map) {
  auto timer = hierarchy_.GetCounterRegistry()->get_timer(
      "cutcell::GriddingAlgorithm::MakeNewLevelFromScratch");
  // Allocate level data.
  const ::amrex::DistributionMapping balanced_distribution_map =
      LoadBalance(level, box_array, distribution_map);
  {
    const int n_comps = hierarchy_.GetDataDescription().n_state_components;
    const int ngrow = hierarchy_.GetOptions().ngrow_eb_level_set;
    std::shared_ptr<::amrex::EBFArrayBoxFactory> eb_factory =
        ::amrex::makeEBFabFactory(
            hierarchy_.GetOptions()
                .index_spaces[static_cast<std::size_t>(level)],
            hierarchy_.GetGeometry(level), box_array, balanced_distribution_map,
            {ngrow, ngrow, ngrow}, ::amrex::EBSupport::full);
    hierarchy_.GetPatchLevel(level) = PatchLevel(
        level, Duration(time_point), box_array, balanced_distribution_map,
        n_comps, std::move(eb_factory), ngrow - 1);
  }

  PatchLevel& patch_level = hierarchy_.GetPatchLevel(level);
  initial_condition_.InitializeData(patch_level, *this, level,
                                    static_cast<Duration>(time_point));
  SetDistributionMap(level, balanced_distribution_map);
}

void GriddingAlgorithm::MakeNewLevelFromCoarse(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_map) {
  auto timer = hierarchy_.GetCounterRegistry()->get_timer(
      "cutcell::GriddingAlgorithm::MakeNewLevelFromCoarse");
  FUB_ASSERT(level > 0);
  const ::amrex::DistributionMapping balanced_distribution_map =
      LoadBalance(level, box_array, distribution_map);
  const PatchLevel& coarse_level = hierarchy_.GetPatchLevel(level - 1);
  const int n_comps = hierarchy_.GetDataDescription().n_state_components;
  const int ngrow = hierarchy_.GetOptions().ngrow_eb_level_set;
  const ::amrex::Geometry& geom = hierarchy_.GetGeometry(level);
  std::shared_ptr<::amrex::EBFArrayBoxFactory> factory =
      ::amrex::makeEBFabFactory(
          hierarchy_.GetOptions().index_spaces[static_cast<std::size_t>(level)],
          geom, box_array, balanced_distribution_map, {ngrow, ngrow, ngrow},
          ::amrex::EBSupport::full);
  PatchLevel fine_level(level, Duration(time_point), box_array,
                        balanced_distribution_map, n_comps, std::move(factory),
                        ngrow - 1);
  const int cons_start = hierarchy_.GetDataDescription().first_cons_component;
  const int n_cons_components =
      hierarchy_.GetDataDescription().n_cons_components;
  ::amrex::Vector<::amrex::BCRec> bcr(static_cast<std::size_t>(n_comps));
  const std::size_t fine = std::size_t(level);
  const std::size_t coarse = std::size_t(level - 1);
  AnyBoundaryCondition& fine_boundary = boundary_condition_[fine];
  AnyBoundaryCondition& coarse_boundary = boundary_condition_[coarse];
  ::amrex::InterpFromCoarseLevel(
      fine_level.data, time_point, coarse_level.data, cons_start, cons_start,
      n_cons_components, hierarchy_.GetGeometry(level - 1),
      hierarchy_.GetGeometry(level), fine_boundary, 0, coarse_boundary, 0,
      {AMREX_D_DECL(2, 2, 2)}, &::amrex::pc_interp, bcr, 0);
  hierarchy_.GetPatchLevel(level) = std::move(fine_level);
  SetDistributionMap(level, balanced_distribution_map);
}

void GriddingAlgorithm::RemakeLevel(
    int level_number, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_map) {
  auto timer = hierarchy_.GetCounterRegistry()->get_timer(
      "cutcell::GriddingAlgorithm::RemakeLevel");
  const ::amrex::DistributionMapping balanced_distribution_map =
      LoadBalance(level_number, box_array, distribution_map);
  const int n_comps = hierarchy_.GetDataDescription().n_state_components;
  const int ngrow = hierarchy_.GetOptions().ngrow_eb_level_set;
  const ::amrex::Geometry& geom = hierarchy_.GetGeometry(level_number);
  std::shared_ptr<::amrex::EBFArrayBoxFactory> factory =
      ::amrex::makeEBFabFactory(
          hierarchy_.GetOptions()
              .index_spaces[static_cast<std::size_t>(level_number)],
          geom, box_array, balanced_distribution_map, {ngrow, ngrow, ngrow},
          ::amrex::EBSupport::full);
  PatchLevel new_level(level_number, Duration(time_point), box_array,
                       balanced_distribution_map, n_comps, std::move(factory),
                       ngrow - 1);
  FillMultiFabFromLevel(new_level.data, level_number);
  hierarchy_.GetPatchLevel(level_number) = std::move(new_level);
  SetDistributionMap(level_number, balanced_distribution_map);
}

void GriddingAlgorithm::ClearLevel(int level) {
  hierarchy_.GetPatchLevel(level) = PatchLevel{};
}

void GriddingAlgorithm::PostProcessBaseGrids(::amrex::BoxArray& ba) const {
  auto timer = hierarchy_.GetCounterRegistry()->get_timer(
      "cutcell::GriddingAlgorithm::PostProcessBaseGrids");
  if (hierarchy_.GetOptions().remove_covered_grids) {
    ::amrex::DistributionMapping dm(ba);
    std::unique_ptr<::amrex::EBFArrayBoxFactory> eb_factory =
        ::amrex::makeEBFabFactory(hierarchy_.GetOptions().index_spaces[0],
                                  hierarchy_.GetGeometry(0), ba, dm, {0, 0, 0},
                                  ::amrex::EBSupport::basic);
    const ::amrex::FabArray<::amrex::EBCellFlagFab>& flags =
        eb_factory->getMultiEBCellFlagFab();
    ::amrex::Vector<::amrex::Box> not_covered{};
    ForEachFab(ba, dm, [&](const ::amrex::MFIter& mfi) {
      if (flags[mfi].getType() != ::amrex::FabType::covered) {
        not_covered.push_back(mfi.validbox());
      }
    });
    ::amrex::AllGatherBoxes(not_covered);
    ba = ::amrex::BoxArray{::amrex::BoxList(std::move(not_covered))};
  }
}

} // namespace fub::amrex::cutcell
