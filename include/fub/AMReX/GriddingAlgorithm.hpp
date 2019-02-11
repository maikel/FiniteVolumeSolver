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

#ifndef FUB_AMREX_GRIDDING_ALGORITHM_HPP
#define FUB_AMREX_GRIDDING_ALGORITHM_HPP

#include "fub/AMReX/FArrayBox.hpp"
#include "fub/AMReX/PatchHierarchy.hpp"

#include <memory>

#include "AMReX_AmrCore.H"

namespace fub {
class TaggingStrategy;

namespace amrex {

template <typename PatchHierarchy, typename InitialData, typename Tagging>
class GriddingAlgorithm : private ::amrex::AmrCore {
public:
  GriddingAlgorithm(PatchHierarchy& hier, InitialData initial_data,
                    Tagging tagging);

  PatchHierarchy& patch_hierarchy() noexcept { return *hierarchy_; }
  const PatchHierarchy& patch_hierarchy() const noexcept { return *hierarchy_; }

  /// Make a finer level on top of `which_level`.
  bool Regrid(int which_level) {
    const int before = AmrMesh::finest_level;
    AmrCore::regrid(which_level,
                    hierarchy_->patch_level(which_level).time_point);
    const int after = AmrMesh::finest_level;
    return before != after;
  }

  void MakeCoarsestLevel(double level_time) {
    ::amrex::AmrCore::MakeNewGrids(level_time);
  }

  void MakeFinerLevel(std::ptrdiff_t level_cycle, double level_time);

  void SetTagging(TaggingStrategy* tagging);

private:
  void ErrorEst(int level, ::amrex::TagBoxArray& tags, double time_point,
                int ngrow) override {
    const ::amrex::Geometry& level_geom = hierarchy_->geometry(level);
    ForEachFAB(
        [&](::amrex::TagBox& tags, auto&& arrays) {
          auto mdspans = ViewDataOnFAB(arrays);
          auto coords = GetCartesianCoordinates(level_geom, GetBox(arrays));
          tagging_.TagCellsForRefinement(MakeMdSpan(tags), mdspans, coords);
        },
        tags, hierarchy_->patch_level(level).data);
  }

  void MakeNewLevelFromScratch(
      int level, double time_point, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override {
    hierarchy_->MakeLevelFromScratch(level, time_point, box_array,
                                     distribution_mapping);
    const ::amrex::Geometry& level_geom = hierarchy_->geometry(level);
    ForEachFAB(
        [&](auto&& arrays) {
          auto mdspans = ViewDataOnFAB(arrays);
          auto coords = GetCartesianCoordinates(level_geom, GetBox(arrays));
          initial_data_.InitializeData(mdspans, coords);
        },
        hierarchy_->patch_level(level).data);
  }

  void MakeNewLevelFromCoarse(
      int level, double time_point, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override {
    hierarchy_->MakeLevelFromCoarse(level, time_point, box_array,
                                    distribution_mapping);
  }

  void RemakeLevel(
      int level, double time_point, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override {
    hierarchy_->RemakeLevel(level, time_point, box_array, distribution_mapping);
  }

  void ClearLevel(int level) override { hierarchy_->ClearLevel(level); }

  PatchHierarchy* hierarchy_;
  InitialData initial_data_;
  Tagging tagging_;
};

template <typename PatchHierarchy, typename InitialData, typename Tagging>
GriddingAlgorithm<PatchHierarchy, InitialData, Tagging>::GriddingAlgorithm(
    PatchHierarchy& hier, InitialData initial_data, Tagging tagging)
    : AmrCore(&hier.geometry().coordinates, hier.options().max_refinement_level,
              ::amrex::Vector<int>(hier.geometry().cell_dimensions.begin(),
                                   hier.geometry().cell_dimensions.end())),
      hierarchy_{std::addressof(hier)},
      initial_data_{initial_data}, tagging_{tagging} {}

// template <typename PatchHierarchy, typename InitialData>
// void InitializeHierarchy(GriddingAlgorithm<PatchHierarchy>& gridding,
//                          InitialData initial_data) {
//   ForEachFAB(
//       [&](::amrex::FArrayBox& data) {
//         const ::amrex::Geometry& geom =
//             gridding.patch_hierarchy().patch_level(0).geometry;
//         CartesianCoordinates coords = GetCartesianCoordinates(geom,
//         data.box()); initial_data.InitializeData(AsStateView(data), coords);
//       },
//       gridding.patch_hierarchy().patch_level(0).data);
//   int level = 0;
//   while (gridding.patch_hierarchy().options().max_refinement_level < level) {
//     if (!gridding.Regrid(level)) {
//       return;
//     }
//     level += 1;
//     ForEachFAB(
//         [&](::amrex::FArrayBox& data) {
//           const ::amrex::Geometry& geom =
//               gridding.patch_hierarchy().patch_level(level).geometry;
//           CartesianCoordinates coords =
//               GetCartesianCoordinates(geom, data.box());
//           initial_data.InitializeData(AsStateView(data), coords);
//         },
//         gridding.patch_hierarchy().patch_level(level));
//   }
// }

} // namespace amrex
} // namespace fub

#endif