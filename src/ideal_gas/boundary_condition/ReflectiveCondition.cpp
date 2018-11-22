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

#include "fub/ideal_gas/boundary_condition/ReflectiveCondition.hpp"
#include "fub/SAMRAI/utility.hpp"

#include "SAMRAI/pdat/CellIndex.h"

namespace fub {
namespace ideal_gas {
namespace {
void ReflectOnBoundaryBox_(
    const IdealGasEquation::CompleteStatePatchData& complete,
    const SAMRAI::hier::BoundaryBox& bbox, const SAMRAI::hier::Box& fill_box) {
  const int locidx = bbox.getLocationIndex();
  const int direction = locidx / 2;
  const int side = locidx % 2;
  const int sign = side == 0 ? +1 : -1;
  for (const SAMRAI::hier::Index& index : fill_box) {
    SAMRAI::pdat::CellIndex to(index);
    SAMRAI::pdat::CellIndex from(shift(index, Direction(direction), sign));
    complete.density(to) = complete.density(from);
    complete.energy(to) = complete.energy(from);
    complete.pressure(to) = complete.pressure(from);
    complete.temperature(to) = complete.temperature(from);
    complete.speed_of_sound(to) = complete.speed_of_sound(from);
    const int n_species = complete.species.getDepth();
    for (int s = 0; s < n_species; ++s) {
      complete.species(to, s) = complete.species(from, s);
    }
    const int dim = complete.momentum.getDepth();
    for (int d = 0; d < dim; ++d) {
      complete.momentum(to, d) = complete.momentum(from, d);
    }
    complete.momentum(to, direction) *= -1;
  }
}
} // namespace

void ReflectiveCondition::setPhysicalBoundaryConditions(
    const SAMRAI::hier::Patch& patch, double fill_time,
    const SAMRAI::hier::IntVector& ghost_width_to_fill) const {
  const SAMRAI::hier::PatchGeometry& geometry = *patch.getPatchGeometry();
  const std::vector<SAMRAI::hier::BoundaryBox>& boundaries =
      geometry.getCodimensionBoundaries(1);
  for (const SAMRAI::hier::BoundaryBox& bbox : boundaries) {
    SAMRAI::hier::Box fill_box =
        geometry.getBoundaryFillBox(bbox, patch.getBox(), ghost_width_to_fill);
    if (fill_box.size() > 0) {
      const int direction = bbox.getLocationIndex() / 2;
      using Scratch = ForwardEulerTimeIntegrator::Scratch;
      IdealGasEquation::CompleteStatePatchData complete =
          integrator_->getCompleteState(patch, Scratch(direction));
      ReflectOnBoundaryBox_(complete, bbox, fill_box);
    }
  }
}

} // namespace ideal_gas
} // namespace fub