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

#ifndef FUB_SOLVER_DIMENSIONAL_SPLIT_INTEGRATOR_HPP
#define FUB_SOLVER_DIMENSIONAL_SPLIT_INTEGRATOR_HPP

#include "fub/SAMRAI/utility.hpp"
#include "fub/solver/InitialCondition.hpp"
#include "fub/solver/BoundaryCondition.hpp"
#include "fub/solver/Direction.hpp"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

namespace fub {
class DimensionalSplitTimeIntegrator {
public:
  virtual ~DimensionalSplitTimeIntegrator() = default;

  void allocatePatchDataOnPatchLevel(SAMRAI::hier::PatchLevel& level) const;

  double computeStableDt(const SAMRAI::hier::PatchHierarchy& hierarchy,
                         double time_point) const;

  void
  fillGhostLayer(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                 const BoundaryCondition& boundary_condition, double time_point,
                 Direction dir) const;

  void
  advanceTime(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
              double time_point, double time_step_size, Direction dir) const;

private:
  virtual void flagPatchDataIdsToAllocate(
      SAMRAI::hier::ComponentSelector& which_to_allocate) const = 0;

  virtual double computeStableDtOnPatch(const SAMRAI::hier::Patch& patch,
                                        double time_point) const = 0;

  virtual void advanceTimeOnPatch(const SAMRAI::hier::Patch& patch,
                                  double time_point, double time_step_size,
                                  Direction dir) const = 0;

  virtual std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>
  getFillGhostLayerRefineAlgorithm(Direction dir) const = 0;

  virtual std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
  getFillGhostLayerCoarsenAlgorithm(Direction dir) const = 0;

  virtual void postprocessAdvanceLevelState(
      const std::shared_ptr<SAMRAI::hier::PatchLevel>& level) const {}
};

void initializePatchHierarchy(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    const DimensionalSplitTimeIntegrator& integrator,
    const InitialCondition& initial_condition);

} // namespace fub

#endif