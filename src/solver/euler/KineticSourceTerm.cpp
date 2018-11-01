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

#include "fub/solver/euler/KineticSourceTerm.hpp"

namespace fub {
namespace euler {
namespace {
IdealGas::CompleteState getCompleteState_(const IdealGas& ideal_gas,
                                          const SAMRAI::hier::Patch& patch) {
  using Variable = IdealGas::Variable;
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = ideal_gas.getCellData(patch, Variable::density);
  SAMRAI::pdat::CellData<double>& momentum = ideal_gas.getCellData(patch, Variable::momentum);
  SAMRAI::pdat::CellData<double>& energy = ideal_gas.getCellData(patch, Variable::energy);
  SAMRAI::pdat::CellData<double>& pressure = ideal_gas.getCellData(patch, Variable::pressure);
  SAMRAI::pdat::CellData<double>& temperature = ideal_gas.getCellData(patch, Variable::temperature);
  SAMRAI::pdat::CellData<double>& speed_of_sound = ideal_gas.getCellData(patch, Variable::speed_of_sound);
  SAMRAI::pdat::CellData<double>& species = ideal_gas.getCellData(patch, Variable::species);
  // clang-format on
  return {density,     momentum,       energy, pressure,
          temperature, speed_of_sound, species};
}
} // namespace

void KineticSourceTerm::advanceTimeOnPatch(
    const std::shared_ptr<SAMRAI::hier::Patch>& patch, double /* time_point */,
    double time_step_size) const {
  IdealGas::CompleteState complete = getCompleteState_(*ideal_gas_, *patch);
  ideal_gas_->advanceSourceTerm(complete, time_step_size);
}

} // namespace euler
} // namespace fub