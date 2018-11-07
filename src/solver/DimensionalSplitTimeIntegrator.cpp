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

#include "fub/solver/DimensionalSplitTimeIntegrator.hpp"

#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/mesh/CascadePartitioner.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/TileClustering.h"

#include "fub/core/assert.hpp"

namespace fub {
double DimensionalSplitTimeIntegrator::computeStableDt(
    const SAMRAI::hier::PatchHierarchy& hierarchy, double time_point,
    Direction dir) const {
  double time_step_size = std::numeric_limits<double>::infinity();
  forEachPatch(hierarchy, [&](const SAMRAI::hier::Patch& patch) {
    const double local_time_step_size =
        this->computeStableDtOnPatch(patch, time_point, dir);
    time_step_size = std::min(time_step_size, local_time_step_size);
  });
  return time_step_size;
}

void DimensionalSplitTimeIntegrator::allocatePatchDataOnPatchLevel(
    SAMRAI::hier::PatchLevel& level) const {
  SAMRAI::hier::ComponentSelector data_to_allocate;
  this->flagPatchDataIdsToAllocate(data_to_allocate);
  level.allocatePatchData(data_to_allocate);
}

void DimensionalSplitTimeIntegrator::advanceTime(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    double time_point, double time_step_size, Direction dir) const {
  forEachPatch(*hierarchy, [&](const SAMRAI::hier::Patch& patch) {
    this->advanceTimeOnPatch(patch, time_point, time_step_size, dir);
  });
}

namespace {
struct BoundaryConditionWrapper : SAMRAI::xfer::RefinePatchStrategy {
  const BoundaryCondition* boundary_condition;

  explicit BoundaryConditionWrapper(const BoundaryCondition& bc)
      : boundary_condition{&bc} {}

  void setPhysicalBoundaryConditions(
      SAMRAI::hier::Patch& patch, const double fill_time,
      const SAMRAI::hier::IntVector& ghost_width_to_fill) override {
    FUB_ASSERT(boundary_condition);
    boundary_condition->setPhysicalBoundaryConditions(patch, fill_time,
                                                      ghost_width_to_fill);
  }

  SAMRAI::hier::IntVector
  getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const override {
    return boundary_condition->getStencilWidth(dim);
  }

  void preprocessRefine(SAMRAI::hier::Patch&, const SAMRAI::hier::Patch&,
                        const SAMRAI::hier::Box&,
                        const SAMRAI::hier::IntVector&) override {}

  void postprocessRefine(SAMRAI::hier::Patch&, const SAMRAI::hier::Patch&,
                         const SAMRAI::hier::Box&,
                         const SAMRAI::hier::IntVector&) override {}
};
} // namespace

void DimensionalSplitTimeIntegrator::fillGhostLayer(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    const BoundaryCondition& boundary_condition, double time_point,
    Direction dir) const {
  std::shared_ptr<SAMRAI::xfer::RefineAlgorithm> algorithm =
      getFillGhostLayerRefineAlgorithm(dir);
  const int max_level_number = hierarchy->getNumberOfLevels();
  BoundaryConditionWrapper bc{boundary_condition};
  for (int level_number = 0; level_number < max_level_number; ++level_number) {
    std::shared_ptr<SAMRAI::hier::PatchLevel> level =
        hierarchy->getPatchLevel(level_number);
    algorithm->createSchedule(level, &bc)->fillData(time_point);
  }
}

void InitializePatchHierarchy(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    const DimensionalSplitTimeIntegrator& integrator,
    const InitialCondition& initial_condition) {
  struct InitializeStrategy : public SAMRAI::mesh::StandardTagAndInitStrategy {
    const DimensionalSplitTimeIntegrator* that;
    const InitialCondition* initial_condition;

    InitializeStrategy(const DimensionalSplitTimeIntegrator& solv,
                       const InitialCondition& init)
        : that{&solv}, initial_condition{&init} {}

    void initializeLevelData(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
        int level_number, double init_data_time, bool can_be_refined,
        bool initial_time,
        const std::shared_ptr<SAMRAI::hier::PatchLevel>& old_level,
        bool allocate_data) override {
      const std::shared_ptr<SAMRAI::hier::PatchLevel>& patch_level =
          hierarchy->getPatchLevel(level_number);
      if (allocate_data) {
        that->allocatePatchDataOnPatchLevel(*patch_level);
      }
      if (initial_condition) {
        for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : *patch_level) {
          initial_condition->InitializeDataOnPatch(*patch);
        }
      }
    }

    void resetHierarchyConfiguration(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
        int coarsest_level, int finest_level) override {}
  };
  InitializeStrategy initialize(integrator, initial_condition);
  const SAMRAI::tbox::Dimension& dim = hierarchy->getDim();
  SAMRAI::mesh::GriddingAlgorithm algorithm(
      hierarchy, "", nullptr,
      std::make_shared<SAMRAI::mesh::StandardTagAndInitialize>("", &initialize),
      std::make_shared<SAMRAI::mesh::TileClustering>(dim),
      std::make_shared<SAMRAI::mesh::CascadePartitioner>(dim, ""));
  algorithm.makeCoarsestLevel(0.0);
}

} // namespace fub