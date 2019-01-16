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

#include "fub/SAMRAI/ideal_gas/InitializeFromData.hpp"
#include "fub/SAMRAI/utility.hpp"

namespace fub {
namespace ideal_gas {

InitializeFromData::InitializeFromData(
    const DimensionalSplitTimeIntegrator& integrator, DynamicMdSpan<double, 2> grid)
    : time_integrator_{&integrator}, data_{grid.data(),
                                           grid.data() + grid.size()},
      view_{data_.data(), grid.extents()} {}

void InitializeFromData::InitializeDataOnPatch(
    const SAMRAI::hier::Patch& patch) const {
  FUB_ASSERT(patch.getDim().getValue() == 1);
  SAMRAI::hier::Box box = patch.getBox();
  using PrimPatchData = IdealGasEquation::PrimPatchData;
  using CompletePatchData = IdealGasEquation::CompletePatchData;
  using Variable = DimensionalSplitTimeIntegrator::Variable;
  PrimPatchData prim =
      time_integrator_->GetEquation()->GetPrimPatchData(patch);
  auto at = [&](int cell, Variable var) -> double& {
    return view_(cell, static_cast<int>(var));
  };
  for (const SAMRAI::hier::Index& index : box) {
    SAMRAI::pdat::CellIndex cell(index);
    prim.temperature(cell) = at(index[0], Variable::temperature);
    prim.momentum(cell) = at(index[0], Variable::momentum);
    prim.pressure(cell) = at(index[0], Variable::pressure);
    const int n_species = time_integrator_->GetEquation()->GetNSpecies();
    auto species = [](int s) { return Variable(int(Variable::species) + s); };
    for (int s = 0; s < n_species; ++s) {
      prim.species(cell) = at(index[0], species(s));
    }
  }
  CompletePatchData complete =
      time_integrator_->GetEquation()->GetCompletePatchData(patch);
  time_integrator_->GetEquation()->FillFromPrim(complete, prim);
}

} // namespace ideal_gas
} // namespace fub