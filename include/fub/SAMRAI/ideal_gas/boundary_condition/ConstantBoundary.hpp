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

#ifndef FUB_EULER_BOUNDARY_CONDITION_CONSTANT_CONDITION_HPP
#define FUB_EULER_BOUNDARY_CONDITION_CONSTANT_CONDITION_HPP

#include "fub/core/mdspan.hpp"
#include "fub/SAMRAI/ideal_gas/DimensionalSplitTimeIntegrator.hpp"
#include "fub/SAMRAI/DimensionalSplitBoundaryCondition.hpp"

namespace fub {
namespace ideal_gas {

class ConstantBoundary : public DimensionalSplitBoundaryCondition {
public:
  ConstantBoundary(std::vector<double> buffer,
                   const DimensionalSplitTimeIntegrator& i)
      : integrator_{&i}, buffer_(std::move(buffer)),
        states_(buffer_.data(), 1, buffer_.size()) {}

  ConstantBoundary(DynamicMdSpan<double, 2> states,
                   const DimensionalSplitTimeIntegrator& i)
      : integrator_{&i},
        buffer_(states.get_span().begin(), states.get_span().end()),
        states_(buffer_.data(), states.extents()) {}

  void SetPhysicalBoundaryCondition(const SAMRAI::hier::Patch& patch,
                                    const SAMRAI::hier::Box& fill_box,
                                    double fill_time, Direction dir,
                                    int side) const override;

  int GetStencilWidth() const override { return 0; }

  void SetStates(std::vector<double> states) {
    buffer_ = std::move(states);
    states_ = DynamicMdSpan<double, 2>(buffer_.data(), 1, buffer_.size());
  }

  void SetStates(span<double> states) {
    if (states.size() == buffer_.size()) {
      std::copy(states.begin(), states.end(), buffer_.begin());
    } else {
      buffer_ = std::vector<double>(states.begin(), states.end());
      states_ = DynamicMdSpan<double, 2>(buffer_.data(), 1, buffer_.size());
    }
  }

private:
  const DimensionalSplitTimeIntegrator* integrator_;
  std::vector<double> buffer_;
  DynamicMdSpan<double, 2> states_;
};

} // namespace ideal_gas
} // namespace fub

#endif