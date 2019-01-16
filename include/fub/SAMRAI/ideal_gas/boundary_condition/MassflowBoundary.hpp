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

#ifndef FUB_EULER_BOUNDARY_CONDITION_MASSFLOW_CONDITION_HPP
#define FUB_EULER_BOUNDARY_CONDITION_MASSFLOW_CONDITION_HPP

#include "fub/core/mdspan.hpp"
#include "fub/SAMRAI/ideal_gas/DimensionalSplitTimeIntegrator.hpp"
#include "fub/SAMRAI/ideal_gas/FlameMasterKinetics.hpp"
#include "fub/SAMRAI/DimensionalSplitBoundaryCondition.hpp"

namespace fub {
namespace ideal_gas {

class MassflowBoundary : public DimensionalSplitBoundaryCondition {
public:
  explicit MassflowBoundary(double required_massflow, double surface_area,
                             const DimensionalSplitTimeIntegrator& i,
                             FlameMasterKinetics& kinetics)
      : integrator_{&i}, kinetics_{&kinetics}, required_massflow_{
                                                   required_massflow}, surface_area_{surface_area} {}

  void SetRequiredMassFlow(double m) {
    required_massflow_ = m;
  }

  double GetRequiredMassFlow(double m) const noexcept {
    return required_massflow_;
  }

  void SetPhysicalBoundaryCondition(const SAMRAI::hier::Patch& patch,
                                    const SAMRAI::hier::Box& fill_box,
                                    double fill_time, Direction dir,
                                    int side) const override;

  int GetStencilWidth() const override { return 1; }

private:
  const DimensionalSplitTimeIntegrator* integrator_;
  FlameMasterKinetics* kinetics_;
  double required_massflow_;
  double surface_area_;
};

} // namespace ideal_gas
} // namespace fub

#endif