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

#include "fub/SAMRAI/ideal_gas/boundary_condition/ConstantBoundary.hpp"

namespace fub {
namespace ideal_gas {

void ConstantBoundary::SetPhysicalBoundaryCondition(
    const SAMRAI::hier::Patch& patch, const SAMRAI::hier::Box& fill_box,
    [[maybe_unused]] double fill_time, Direction dir, int side) const {
  using Scratch = DimensionalSplitTimeIntegrator::Scratch;
  IdealGasEquation::CompletePatchData complete =
      integrator_->GetComletePatchData(patch, Scratch(dir));
  const int d = static_cast<int>(dir);
  const int n_species = complete.species.getDepth();
  using Variable = IdealGasEquation::Variable;
  const int s0 = static_cast<int>(Variable::species);
  auto state = [&](int i, Variable var) {
    return states_(i, static_cast<int>(var));
  };
  for (const SAMRAI::hier::Index& index : fill_box) {
    SAMRAI::pdat::CellIndex cell(index);
    const int i0 = side == 0 ? index[0] - fill_box.lower()[0]
                             : fill_box.upper()[0] - index[0];
    complete.density(cell) = state(i0, Variable::density);
    complete.momentum(cell) = state(i0, Variable::momentum);
    complete.energy(cell) = state(i0, Variable::energy);
    complete.pressure(cell) = state(i0, Variable::pressure);
    complete.temperature(cell) = state(i0, Variable::temperature);
    complete.speed_of_sound(cell) = state(i0, Variable::speed_of_sound);
    for (int s = 0; s < n_species; ++s) {
      complete.species(cell, s) = states_(i0, s0 + s);
    }
  }
}

} // namespace ideal_gas
} // namespace fub