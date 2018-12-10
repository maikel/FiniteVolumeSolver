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

#include "fub/ideal_gas/KineticDriver.hpp"
#include "fub/ideal_gas/boundary_condition/ReflectiveCondition.hpp"

namespace fub {
namespace ideal_gas {

KineticDriver::KineticDriver(const std::shared_ptr<IdealGasKinetics>& equation,
                             const CoordinateRange& coordinates,
                             std::ptrdiff_t n_cells)
    : time_integrator_{equation}, source_term_{equation},
      boundary_condition_(std::make_shared<SplitBoundaryConditions>()),
      hierarchy_(MakeCartesianPatchHierarchy(
          IndexRange{
              SAMRAI::hier::Index(equation->GetDimension(), 0),
              SAMRAI::hier::Index(equation->GetDimension(), n_cells - 1)},
          coordinates)) {
  std::shared_ptr<const ideal_gas::ReflectiveCondition> cond =
      std::make_shared<ideal_gas::ReflectiveCondition>(time_integrator_);
  boundary_condition_->SetBoundaryCondition(cond, Direction::X, 0);
  boundary_condition_->SetBoundaryCondition(cond, Direction::X, 1);
}

void KineticDriver::InitializeHierarchy(const InitialCondition& init) {
  fub::InitializePatchHierarchy(hierarchy_, time_integrator_, init);
}

void KineticDriver::SetLeftBoundaryCondition(
    const std::shared_ptr<const SplitBoundaryCondition>& condition) {
  boundary_condition_->SetBoundaryCondition(condition, Direction::X, 0);
}

void KineticDriver::SetRightBoundaryCondition(
    const std::shared_ptr<const SplitBoundaryCondition>& condition) {
  boundary_condition_->SetBoundaryCondition(condition, Direction::X, 1);
}

double KineticDriver::ComputeStableDt() const {
  const double dt =
      solver_.computeStableDt(hierarchy_, *boundary_condition_, time_point_);
  return 0.5 * dt;
}

void KineticDriver::Advance(double dt) {
  double relt = 0.0;
  while (relt < dt) {
    const double stable_dt = ComputeStableDt();
    const double limit_dt = dt - relt;
    FUB_ASSERT(limit_dt > 0.0);
    const double actual_dt = std::min(stable_dt, limit_dt);
    solver_.advanceTime(hierarchy_, *boundary_condition_, time_point_,
                        actual_dt);
    time_point_ += actual_dt;
  }
}

} // namespace ideal_gas
} // namespace fub