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

#include "fub/SAMRAI/GriddingAlgorithm.hpp"

#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/CascadePartitioner.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/mesh/TileClustering.h>
#include <SAMRAI/pdat/CellDoubleConstantRefine.h>

namespace fub {
namespace samrai {
namespace {
SAMRAI::hier::ComponentSelector SelectComponents(span<const int> data_ids) {
  SAMRAI::hier::ComponentSelector selector;
  for (int id : data_ids) {
    selector.setFlag(id);
  }
  return selector;
}

std::vector<SAMRAI::pdat::CellData<double>*>
GetCellData(SAMRAI::hier::Patch& patch, span<const int> ids) {
  std::vector<SAMRAI::pdat::CellData<double>*> data;
  data.reserve(ids.size());
  for (int id : ids) {
    data.push_back(static_cast<SAMRAI::pdat::CellData<double>*>(
        patch.getPatchData(id).get()));
  }
  return data;
}

struct TagAndInit : SAMRAI::mesh::TagAndInitializeStrategy {
  TagAndInit(GriddingAlgorithm* parent,
             SAMRAI::hier::ComponentSelector which_to_allocate)
      : TagAndInitializeStrategy("TagAndInit"), parent_{parent},
        which_to_allocate_{std::move(which_to_allocate)} {}

  bool getUserSuppliedRefineBoxes(SAMRAI::hier::BoxContainer&, int, int,
                                  double) override {
    return false;
  }

  void resetRefineBoxes(const SAMRAI::hier::BoxContainer&, int) override {}

  void initializeLevelData(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      int level_number, double init_data_time, bool, bool initial_time,
      const std::shared_ptr<SAMRAI::hier::PatchLevel>& old_level,
      bool allocate_data) override {
    std::shared_ptr<SAMRAI::hier::PatchLevel> level =
        hierarchy->getPatchLevel(level_number);
    if (allocate_data) {
      level->allocatePatchData(which_to_allocate_, init_data_time);
    }
    if (initial_time) {
      parent_->GetInitialData().InitializeData(
          parent_->GetPatchHierarchy(), level_number, Duration(init_data_time));
    } else {
      // Fill from old_level and refine from coarser level where needed
      if (level_number > 0) {
        const int coarser_level_number = level_number - 1;
        parent_->refine_data_algorithm_
            ->createSchedule(level, old_level, coarser_level_number, hierarchy)
            ->fillData(init_data_time);
      }
    }
  }

  void tagCellsForRefinement(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      int level_number, int /* regrid_cycle */, double regrid_time,
      int tag_index, bool /* initial_time */, bool /* coarsest_sync_level */,
      bool can_be_refined, double /* regrid_start_time */) override {
    std::shared_ptr<SAMRAI::hier::PatchLevel> level =
        hierarchy->getPatchLevel(level_number);
    if (can_be_refined) {
      for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : *level) {
        SAMRAI::pdat::CellData<int>* tags =
            static_cast<SAMRAI::pdat::CellData<int>*>(
                patch->getPatchData(tag_index).get());
        tags->fill(0);
        parent_->tagging_.TagCellsForRefinement(*parent_, level_number, tag_index, Duration(regrid_time));
      }
    }
  }

  void preprocessErrorEstimation(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&, int, int, double,
      double, bool) override {}

  bool usesTimeIntegration(int, double) override { return false; }

  bool everUsesTimeIntegration() const override { return false; }

  bool coarsestLevelBoxesOK(const SAMRAI::hier::BoxContainer&) const override {
    return true;
  }

  int getErrorCoarsenRatio() const override { return 1; }

  void
  checkCoarsenRatios(const std::vector<SAMRAI::hier::IntVector>&) override {}

  bool refineUserBoxInputOnly(int, double) override { return false; }

  void processHierarchyBeforeAddingNewLevel(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&, int,
      const std::shared_ptr<SAMRAI::hier::BoxLevel>&) override {}

  void processLevelBeforeRemoval(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&, int,
      const std::shared_ptr<SAMRAI::hier::PatchLevel>&) override {}

  void resetHierarchyConfiguration(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&, int, int) override {
  }

  GriddingAlgorithm* parent_;
  SAMRAI::hier::ComponentSelector which_to_allocate_;
};
} // namespace

GriddingAlgorithm::GriddingAlgorithm(PatchHierarchy hier, InitialData init,
                                     Tagging tag, std::vector<int> buffer)
    : hierarchy_{hier}, initial_data_{std::move(init)}, tagging_{std::move(
                                                            tag)},
      id_set_{hier.GetDataDescription()}, tag_buffer_{std::move(buffer)} {
  const SAMRAI::hier::ComponentSelector which_to_allocate =
      SelectComponents(hier.GetDataDescription().data_ids);
  const SAMRAI::tbox::Dimension dim = hier.GetNative()->getDim();
  for (int level = 0; level < hier.GetNative()->getMaxNumberOfLevels();
       ++level) {
    hier.GetNative()->setRatioToCoarserLevel(SAMRAI::hier::IntVector(dim, 2),
                                             level);
  }

  refine_data_algorithm_ = std::make_shared<SAMRAI::xfer::RefineAlgorithm>();
  for (int id : id_set_.data_ids) {
    refine_data_algorithm_->registerRefine(
        id, id, id, std::make_shared<SAMRAI::pdat::CellDoubleConstantRefine>());
  }

  algorithm_ = std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
      hier.GetNative(), "GriddingAlgorithm", nullptr,
      std::make_shared<TagAndInit>(this, which_to_allocate),
      std::make_shared<SAMRAI::mesh::TileClustering>(dim),
      std::make_shared<SAMRAI::mesh::CascadePartitioner>(dim,
                                                         "CascadePartitioner"));
}

const PatchHierarchy& GriddingAlgorithm::GetPatchHierarchy() const noexcept {
  return hierarchy_;
}

PatchHierarchy& GriddingAlgorithm::GetPatchHierarchy() noexcept {
  return hierarchy_;
}

const InitialData& GriddingAlgorithm::GetInitialData() const noexcept {
  return initial_data_;
}

InitialData& GriddingAlgorithm::GetInitialData() noexcept {
  return initial_data_;
}

void GriddingAlgorithm::RegridAllFinerLevels(int level_num, int cycle,
                                             Duration time_point) {
  algorithm_->regridAllFinerLevels(level_num, tag_buffer_, cycle,
                                   time_point.count());
}

void GriddingAlgorithm::InitializeHierarchy(Duration initial_time_point,
                                            std::ptrdiff_t initial_cycle) {
  algorithm_->makeCoarsestLevel(initial_time_point.count());
  bool initial_time = true;
  for (int buffer : tag_buffer_) {
    algorithm_->makeFinerLevel(buffer, initial_time, initial_cycle,
                               initial_time_point.count());
  }
}

} // namespace samrai
} // namespace fub
