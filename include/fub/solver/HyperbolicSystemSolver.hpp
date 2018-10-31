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

#ifndef FUB_SOLVER_HYPERBOLIC_SYSTEM_SOLVER_HPP
#define FUB_SOLVER_HYPERBOLIC_SYSTEM_SOLVER_HPP

#include "fub/solver/DimensionalSplitTimeIntegrator.hpp"
#include "fub/solver/SplittingMethod.hpp"

namespace fub {

class HyperbolicSystemSolver {
public:
  HyperbolicSystemSolver(const DimensionalSplitTimeIntegrator& integrator,
                         const SplittingMethod& splitting)
      : integrator_{&integrator}, splitting_method_{&splitting} {}

  double ComputeStableDt(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      const BoundaryCondition& boundary_condition, double time_point) const;

  void
  AdvanceTime(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
              const BoundaryCondition& boundary_condition, double time_point,
              double time_step_size) const;

  const DimensionalSplitTimeIntegrator& GetTimeIntegrator() const noexcept {
    return *integrator_;
  }

  const SplittingMethod& GetSplittingMethod() const noexcept {
    return *splitting_method_;
  }

private:
  const DimensionalSplitTimeIntegrator* integrator_;
  const SplittingMethod* splitting_method_;
};

} // namespace fub

#endif