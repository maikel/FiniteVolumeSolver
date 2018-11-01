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

#include "fub/solver/HyperbolicSystemSolver.hpp"

namespace fub {

double HyperbolicSystemSolver::computeStableDt(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    const BoundaryCondition& boundary_condition, double time_point) const {
  integrator_->fillGhostLayer(hierarchy, boundary_condition, time_point,
                              Direction::X);
  const double dt = integrator_->computeStableDt(*hierarchy, time_point);
  return dt / hierarchy->getDim().getValue();
}

void HyperbolicSystemSolver::advanceTime(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    const BoundaryCondition& boundary_condition, double time_point,
    double time_step_size) const {
  const int dim = hierarchy->getDim().getValue();
  if (dim == 1) {
    integrator_->fillGhostLayer(hierarchy, boundary_condition, time_point,
                                Direction::X);
    integrator_->advanceTime(hierarchy, time_point, time_step_size,
                             Direction::X);
  } else if (dim == 2) {
    splitting_method_->advanceTime(
        hierarchy, time_point, time_step_size,
        [&](const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            double time_point, double time_step_size) {
          integrator_->fillGhostLayer(hierarchy, boundary_condition, time_point,
                                      Direction::X);
          integrator_->advanceTime(hierarchy, time_point, time_step_size,
                                   Direction::X);
        },
        [&](const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            double time_point, double time_step_size) {
          integrator_->fillGhostLayer(hierarchy, boundary_condition, time_point,
                                      Direction::Y);
          integrator_->advanceTime(hierarchy, time_point, time_step_size,
                                   Direction::Y);
        });
  } else if (dim == 3) {
    splitting_method_->advanceTime(
        hierarchy, time_point, time_step_size,
        [&](const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            double time_point, double time_step_size) {
          integrator_->fillGhostLayer(hierarchy, boundary_condition, time_point,
                                      Direction::X);
          integrator_->advanceTime(hierarchy, time_point, time_step_size,
                                   Direction::X);
        },
        [&](const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            double time_point, double time_step_size) {
          integrator_->fillGhostLayer(hierarchy, boundary_condition, time_point,
                                      Direction::Y);
          integrator_->advanceTime(hierarchy, time_point, time_step_size,
                                   Direction::Y);
        },
        [&](const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            double time_point, double time_step_size) {
          integrator_->fillGhostLayer(hierarchy, boundary_condition, time_point,
                                      Direction::Z);
          integrator_->advanceTime(hierarchy, time_point, time_step_size,
                                   Direction::Z);
        });
  } else {
    throw std::runtime_error(
        "No Dimensional Splitting for dim > 3 implemented.");
  }
}

} // namespace fub