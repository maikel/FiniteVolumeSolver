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

#include "fub/solver/DimensionalSplitIntegrator.hpp"
#include "SAMRAI/hier/CoarseFineBoundary.h"

namespace fub {
std::shared_ptr<SAMRAI::hier::PatchHierarchy>
DimensionalSplitIntegrator::makePatchHierarchy(
    const PatchHierarchyOptions& options) const {
  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy =
      makeCartesianPatchHierarchy(options.indices, options.coordinates);
  const SAMRAI::tbox::Dimension& dim = hierarchy->getDim();
  hierarchy->setMaxNumberOfLevels(options.max_number_of_levels);
  const SAMRAI::hier::IntVector ratio_to_coarser(dim, 2);
  for (int level = 0; level < options.max_number_of_levels; ++level) {
    hierarchy->setRatioToCoarserLevel(ratio_to_coarser, level);
  }
  initializePatchHierarchy(hierarchy, 0.0);
  return hierarchy;
}

double DimensionalSplitIntegrator::estimateHierarchyTimeStepSize(
    const SAMRAI::hier::PatchHierarchy& hierarchy, double time_point) const {
  double time_step_size = std::numeric_limits<double>::infinity();
  forEachPatch(hierarchy, [&](const SAMRAI::hier::Patch& patch) {
    const double local_time_step_size =
        this->estimateStableTimeStepSize(patch, time_point);
    time_step_size = std::min(time_step_size, local_time_step_size);
  });
  return time_step_size;
}

namespace {
template <typename FixPatch>
void fixHierarchyConservationOnCoarseFineBoundary(
    const SAMRAI::hier::PatchHierarchy& hierarchy,
    const SAMRAI::hier::IntVector& gcw, FixPatch fix_patch) {
  const int max_level_number = hierarchy.getNumberOfLevels();
  for (int level_number = 0; level_number < max_level_number; ++level_number) {
    const SAMRAI::hier::PatchLevel& level =
        *hierarchy.getPatchLevel(level_number);
    SAMRAI::hier::CoarseFineBoundary coarse_fine_boundary(hierarchy,
                                                          level_number, gcw);
    for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
      const std::vector<SAMRAI::hier::BoundaryBox>& boundaries =
          coarse_fine_boundary.getBoundaries(patch->getGlobalId(), 1);
      for (const SAMRAI::hier::BoundaryBox& boundary : boundaries) {
        fix_patch(*patch, boundary);
      }
    }
  }
}
} // namespace

void DimensionalSplitIntegrator::integratePatchHierarchy(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    double time_point, double time_step_size, int dir) const {
  forEachPatch(*hierarchy, [&](const SAMRAI::hier::Patch& patch) {
    this->integratePatch(patch, time_point, time_step_size, dir);
  });
  coarsenFluxOnCoarseFineInteraces(hierarchy);
  fixHierarchyConservationOnCoarseFineBoundary(
      *hierarchy, SAMRAI::hier::IntVector::getZero(hierarchy->getDim()),
      [&](const SAMRAI::hier::Patch& patch,
          const SAMRAI::hier::BoundaryBox& boundary) {
        this->fixPatchConservationOnCoarseFineBoundary(patch, boundary,
                                                       time_step_size);
      });
  postIntegratePatchHierarchy(hierarchy, time_point, time_step_size, dir);
}

void DimensionalSplitIntegrator::fillGhostLayer(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    double time_point, int dir) const {
  doFillGhostLayer(hierarchy, time_point, dir);
}

} // namespace fub