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

#ifndef FUB_EULER_BOUNDARY_CONDITION_ISENTROPIC_EXPANSION_CONDITION_HPP
#define FUB_EULER_BOUNDARY_CONDITION_ISENTROPIC_EXPANSION_CONDITION_HPP

#include "fub/core/span.hpp"
#include "fub/SAMRAI/ideal_gas/FlameMasterKinetics.hpp"
#include "fub/SAMRAI/ideal_gas/DimensionalSplitTimeIntegrator.hpp"
#include "fub/SAMRAI/DimensionalSplitBoundaryCondition.hpp"

#include <vector>

namespace fub {
namespace ideal_gas {

class IsentropicExpansionBoundary : public DimensionalSplitBoundaryCondition {
public:
  explicit IsentropicExpansionBoundary(double pressure,
                                       const DimensionalSplitTimeIntegrator& i,
                                       FlameMasterKinetics& kinetics)
      : integrator_{&i}, kinetics_{&kinetics}, pressure_{pressure} {}

  void SetPhysicalBoundaryCondition(const SAMRAI::hier::Patch& patch,
                                    const SAMRAI::hier::Box& fill_box,
                                    double fill_time, Direction dir,
                                    int side) const override;

  int GetStencilWidth() const override { return 0; }

  void SetMassFractions(span<const double> fractions) {
    mass_fractions_ = std::vector<double>(fractions.begin(), fractions.end());
  }

  const FlameMasterKinetics& GetEquation() const noexcept {
    return *kinetics_;
  }

  FlameMasterKinetics& GetEquation() noexcept {
    return *kinetics_;
  }

private:
  const DimensionalSplitTimeIntegrator* integrator_;
  FlameMasterKinetics* kinetics_;
  double pressure_;
  std::vector<double> mass_fractions_{};
};

} // namespace ideal_gas
} // namespace fub

#endif