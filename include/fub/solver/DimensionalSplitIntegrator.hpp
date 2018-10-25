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
#include "SAMRAI/hier/Patch.h"

namespace fub {
class DimensionalSplitIntegrator {
public:
  struct PatchHierarchyOptions {
    IndexRange indices;
    CoordinateRange coordinates;
    int max_number_of_levels;
  };

  virtual ~DimensionalSplitIntegrator() = default;

  std::shared_ptr<SAMRAI::hier::PatchHierarchy>
  makePatchHierarchy(const PatchHierarchyOptions& options) const;

  double
  estimateHierarchyTimeStepSize(const SAMRAI::hier::PatchHierarchy& hierarchy,
                                double time_point) const;

  void
  fillGhostLayer(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                 double time_point, int dir) const;

  void integratePatchHierarchy(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      double time_point, double time_step_size, int dir) const;

private:
  virtual void initializePatchHierarchy(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      double time_point) const = 0;

  virtual void doFillGhostLayer(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      double time_point, int dir) const = 0;

  virtual double estimateStableTimeStepSize(const SAMRAI::hier::Patch& patch,
                                            double time_point) const = 0;

  virtual void integratePatch(const SAMRAI::hier::Patch& patch,
                              double time_point, double time_step_size,
                              int dir) const = 0;

  virtual void coarsenFluxOnCoarseFineInteraces(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy) const = 0;

  virtual void fixPatchConservationOnCoarseFineBoundary(
      const SAMRAI::hier::Patch& patch,
      const SAMRAI::hier::BoundaryBox& boundary,
      double time_step_size) const = 0;

  virtual void postIntegratePatchHierarchy(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      double time_point, double time_step_size, int dir) const {};
};

} // namespace fub

#endif