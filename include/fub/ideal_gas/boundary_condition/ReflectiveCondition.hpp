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

#ifndef FUB_EULER_BOUNDARY_CONDITION_REFLECTIVE_CONDITION_HPP
#define FUB_EULER_BOUNDARY_CONDITION_REFLECTIVE_CONDITION_HPP

#include "fub/ideal_gas/ForwardEulerTimeIntegrator.hpp"
#include "fub/solver/BoundaryCondition.hpp"

namespace fub {
namespace ideal_gas {

class ReflectiveCondition : public BoundaryCondition {
public:
  explicit ReflectiveCondition(const ForwardEulerTimeIntegrator& i)
      : integrator_{&i} {}

  void setPhysicalBoundaryConditions(
      const SAMRAI::hier::Patch& patch, double fill_time,
      const SAMRAI::hier::IntVector& ghost_width_to_fill) const override;

  SAMRAI::hier::IntVector
  getStencilWidth(const SAMRAI::tbox::Dimension& dim) const override {
    return SAMRAI::hier::IntVector::getZero(dim);
  }

private:
  const ForwardEulerTimeIntegrator* integrator_;
};

} // namespace ideal_gas
} // namespace fub

#endif