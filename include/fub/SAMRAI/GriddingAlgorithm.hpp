// Copyright (c) 2018 Maikel Nadolski
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

#ifndef FUB_SAMRAI_GRIDDING_ALGORITHM_HPP
#define FUB_SAMRAI_GRIDDING_ALGORITHM_HPP

#include "fub/SAMRAI/PatchView.hpp"
#include "fub/SAMRAI/RegisterVariables.hpp"

#include <SAMRAI/mesh/CascadePartitioner.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/mesh/TileClustering.h>

namespace fub {
namespace samrai {
inline SAMRAI::hier::ComponentSelector
SelectComponents(span<const int> data_ids) {
  SAMRAI::hier::ComponentSelector selector;
  for (int id : data_ids) {
    selector.setFlag(id);
  }
  return selector;
}

template <typename Equation, typename InitialData, typename TaggingStrategy>
struct InitialTagAndInit : SAMRAI::mesh::TagAndInitializeStrategy {
  InitialTagAndInit(std::vector<int> data_ids,
                    SAMRAI::hier::ComponentSelector which_to_allocate,
                    InitialData initial_data, TaggingStrategy tagging)
      : TagAndInitializeStrategy("InitialTagAndInit"), data_ids_{std::move(
                                                           data_ids)},
        which_to_allocate_{std::move(which_to_allocate)},
        initial_data_{std::move(initial_data)}, tagging_{std::move(tagging)} {}

  bool getUserSuppliedRefineBoxes(SAMRAI::hier::BoxContainer&, int, int,
                                  double) override {
    return false;
  }

  void resetRefineBoxes(const SAMRAI::hier::BoxContainer&, int) override {}

  void initializeLevelData(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      int level_number, double init_data_time, bool, bool,
      const std::shared_ptr<SAMRAI::hier::PatchLevel>&,
      bool allocate_data) override {
    std::shared_ptr<SAMRAI::hier::PatchLevel> level =
        hierarchy->getPatchLevel(level_number);
    if (allocate_data) {
      level->allocatePatchData(which_to_allocate_, init_data_time);
    }
    for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : *level) {
      using StateView = StateView<double, Equation>;
      std::array<int, StateView::size()> ids;
      std::copy(data_ids_.begin(), data_ids_.end(), ids.begin());
      auto states = ViewCellDataOnPatch<StateView>(*patch, ids);
      initial_data_.InitializeData(states, GetCartesianCoordinates(*patch));
    }
  }

  void tagCellsForRefinement(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      int level_number, int /* regrid_cycle */, double /* regrid_time */,
      int tag_index, bool /* initial_time */, bool /* coarsest_sync_level */,
      bool can_be_refined, double /* regrid_start_time */) override {
    std::shared_ptr<SAMRAI::hier::PatchLevel> level =
        hierarchy->getPatchLevel(level_number);
    if (can_be_refined) {
      for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : *level) {
        using StateView = StateView<double, Equation>;
        auto cell_data = static_cast<SAMRAI::pdat::CellData<int>*>(
            patch->getPatchData(tag_index).get());
        auto tags = MakeMdSpan(boost::hana::type_c<dynamic_mdspan<int, 3>>,
                               cell_data->getArrayData());
        std::array<int, StateView::size()> ids;
        std::copy(data_ids_.begin(), data_ids_.end(), ids.begin());
        auto states = ViewCellDataOnPatch<StateView>(*patch, ids);
        auto coords = GetCartesianCoordinates(*patch);
        tagging_.TagCellsForRefinement(tags, states, coords);
      }
    }
  }

  void preprocessErrorEstimation(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      int level_number, int cycle, double regrid_time, double regrid_start_time,
      bool initial_time) override {}

  bool usesTimeIntegration(int cycle, double time) override { return false; }

  bool everUsesTimeIntegration() const override { return false; }

  bool coarsestLevelBoxesOK(const SAMRAI::hier::BoxContainer&) const override {
    return true;
  }

  int getErrorCoarsenRatio() const override { return 1; }

  void
  checkCoarsenRatios(const std::vector<SAMRAI::hier::IntVector>&) override {}

  bool refineUserBoxInputOnly(int cycle, double time) override { return false; }

  void processHierarchyBeforeAddingNewLevel(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      int level_number,
      const std::shared_ptr<SAMRAI::hier::BoxLevel>& new_box_level) override {}

  void processLevelBeforeRemoval(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      int level_number,
      const std::shared_ptr<SAMRAI::hier::PatchLevel>& old_level) override {}

  void resetHierarchyConfiguration(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&, int, int) override {
  }

  std::vector<int> data_ids_;
  SAMRAI::hier::ComponentSelector which_to_allocate_;
  InitialData initial_data_;
  TaggingStrategy tagging_;
};

template <typename Equation, typename InitialData, typename TaggingStrategy>
std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm>
GriddingAlgorithm(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hier,
                  const Equation&, span<const int> data_ids,
                  InitialData initial_data, TaggingStrategy tagging) {
  const SAMRAI::hier::ComponentSelector which_to_allocate =
      SelectComponents(data_ids);
  const SAMRAI::tbox::Dimension dim = hier->getDim();
  auto tag_and_init = std::make_shared<
      InitialTagAndInit<Equation, InitialData, TaggingStrategy>>(
      std::vector<int>(data_ids.begin(), data_ids.end()), which_to_allocate,
      std::move(initial_data), std::move(tagging));
  return std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
      hier, "GriddingAlgorithm", nullptr, tag_and_init,
      std::make_shared<SAMRAI::mesh::TileClustering>(dim),
      std::make_shared<SAMRAI::mesh::CascadePartitioner>(dim,
                                                         "CascadePartitioner"));
}

} // namespace samrai
} // namespace fub

#endif