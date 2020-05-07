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
#include <SAMRAI/tbox/MemoryDatabase.h>

namespace fub {
namespace samrai {
namespace {
struct TagAndInit : SAMRAI::mesh::TagAndInitializeStrategy {
  TagAndInit(GriddingAlgorithm* parent, AnyInitialData* init,
             AnyTaggingMethod* tagging_method,
             std::shared_ptr<SAMRAI::xfer::RefineAlgorithm> refine,
             SAMRAI::hier::ComponentSelector which_to_allocate)
      : TagAndInitializeStrategy("TagAndInit"), parent_{parent},
        initial_data_(init), tagging_method_(tagging_method),
        refine_data_algorithm_{std::move(refine)}, which_to_allocate_{std::move(
                                                       which_to_allocate)} {}

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
      initial_data_->InitializeData(level, *parent_, level_number,
                                    Duration(init_data_time));
    } else {
      // Fill from old_level and refine from coarser level where needed
      if (level_number > 0) {
        const int coarser_level_number = level_number - 1;
        refine_data_algorithm_
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
        tagging_method_->TagCellsForRefinement(
            tag_index, *parent_, level_number, Duration(regrid_time));
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
  AnyInitialData* initial_data_;
  AnyTaggingMethod* tagging_method_;
  std::shared_ptr<SAMRAI::xfer::RefineAlgorithm> refine_data_algorithm_;
  SAMRAI::hier::ComponentSelector which_to_allocate_;
};
} // namespace

GriddingAlgorithm::GriddingAlgorithm(PatchHierarchy hier, AnyInitialData init,
                                     AnyTaggingMethod tag,
                                     std::vector<int> buffer)
    : hierarchy_{std::move(hier)}, initial_data_{std::move(init)},
      tagging_method_{std::move(tag)}, tag_buffer_{std::move(buffer)} {
  const SAMRAI::hier::ComponentSelector which_to_allocate =
      SelectComponents(hierarchy_.GetDataDescription().data_ids);
  const SAMRAI::tbox::Dimension dim = hierarchy_.GetNative()->getDim();
  for (int level = 1; level < hierarchy_.GetNative()->getMaxNumberOfLevels();
       ++level) {
    hierarchy_.GetNative()->setRatioToCoarserLevel(
        SAMRAI::hier::IntVector(dim, 2), level);
  }

  refine_data_algorithm_ = std::make_shared<SAMRAI::xfer::RefineAlgorithm>();
  for (int id : hierarchy_.GetDataDescription().data_ids) {
    refine_data_algorithm_->registerRefine(
        id, id, id, std::make_shared<SAMRAI::pdat::CellDoubleConstantRefine>());
  }

  std::shared_ptr cascade_options =
      std::make_shared<SAMRAI::tbox::MemoryDatabase>("CascadePartitioner");
  cascade_options->putVector("tile_size", std::vector{8, 8});

  algorithm_ = std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
      hierarchy_.GetNative(), MakeUniqueName(), nullptr,
      std::make_shared<TagAndInit>(this, &initial_data_, &tagging_method_,
                                   refine_data_algorithm_, which_to_allocate),
      std::make_shared<SAMRAI::mesh::TileClustering>(dim),
      std::make_shared<SAMRAI::mesh::CascadePartitioner>(dim, MakeUniqueName(),
                                                         cascade_options));
}

GriddingAlgorithm::GriddingAlgorithm(const GriddingAlgorithm& ga)
    : hierarchy_(ga.GetPatchHierarchy()), initial_data_(ga.GetInitialData()),
      tagging_method_(ga.GetTaggingMethod()),
      boundary_condition_(ga.GetBoundaryCondition()),
      tag_buffer_(ga.GetTagBuffer()) {
  const SAMRAI::hier::ComponentSelector which_to_allocate =
      SelectComponents(hierarchy_.GetDataDescription().data_ids);
  const SAMRAI::tbox::Dimension dim = hierarchy_.GetNative()->getDim();
  for (int level = 1; level < hierarchy_.GetNative()->getMaxNumberOfLevels();
       ++level) {
    hierarchy_.GetNative()->setRatioToCoarserLevel(
        SAMRAI::hier::IntVector(dim, 2), level);
  }
  refine_data_algorithm_ = std::make_shared<SAMRAI::xfer::RefineAlgorithm>();
  for (int id : hierarchy_.GetDataDescription().data_ids) {
    refine_data_algorithm_->registerRefine(
        id, id, id, std::make_shared<SAMRAI::pdat::CellDoubleConstantRefine>());
  }

  algorithm_ = std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
      hierarchy_.GetNative(), MakeUniqueName(), nullptr,
      std::make_shared<TagAndInit>(this, &initial_data_, &tagging_method_,
                                   refine_data_algorithm_, which_to_allocate),
      std::make_shared<SAMRAI::mesh::TileClustering>(dim),
      std::make_shared<SAMRAI::mesh::CascadePartitioner>(dim,
                                                         MakeUniqueName()));
}

const PatchHierarchy& GriddingAlgorithm::GetPatchHierarchy() const noexcept {
  return hierarchy_;
}

PatchHierarchy& GriddingAlgorithm::GetPatchHierarchy() noexcept {
  return hierarchy_;
}

const AnyInitialData& GriddingAlgorithm::GetInitialData() const noexcept {
  return initial_data_;
}

const AnyTaggingMethod& GriddingAlgorithm::GetTaggingMethod() const noexcept {
  return tagging_method_;
}

const std::vector<int>& GriddingAlgorithm::GetTagBuffer() const noexcept {
  return tag_buffer_;
}

Duration GriddingAlgorithm::GetTimePoint() const noexcept {
  return hierarchy_.GetTimePoint();
}

std::ptrdiff_t GriddingAlgorithm::GetCycles() const noexcept {
  return hierarchy_.GetCycles();
}

void GriddingAlgorithm::RegridAllFinerLevels(int level_num) {
  Duration time_point = hierarchy_.GetTimePoint(level_num);
  int cycles = hierarchy_.GetCycles(level_num);
  if (level_num < 0) {
    algorithm_->makeCoarsestLevel(time_point.count());
  } else {
    algorithm_->regridAllFinerLevels(level_num, tag_buffer_, cycles,
                                     time_point.count());
  }
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

const AnyBoundaryCondition& GriddingAlgorithm::GetBoundaryCondition() const
    noexcept {
  return boundary_condition_;
}

AnyBoundaryCondition& GriddingAlgorithm::GetBoundaryCondition() noexcept {
  return boundary_condition_;
}

} // namespace samrai
} // namespace fub
