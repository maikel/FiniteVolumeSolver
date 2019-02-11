// Copyright (c) 2018-2019 Maikel Nadolski
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

#include "fub/SAMRAI/ideal_gas/MusclHancockMethod.hpp"
#include "fub/SAMRAI/ideal_gas/HlleGodunovMethod.hpp"

namespace fub {
namespace ideal_gas {
double
MusclHancockMethod::ComputeStableDtOnPatch(const CompletePatchData& state,
                                           const SAMRAI::hier::Patch& patch,
                                           Direction dir) const {
  return HlleGodunovMethod::ComputeStableDtOnPatch(state, patch, dir);
}

SAMRAI::hier::IntVector
MusclHancockMethod::GetStencilWidth(const SAMRAI::tbox::Dimension& dim) {
  return SAMRAI::hier::IntVector(dim, 2);
}

MusclHancockMethod::MusclHancockMethod(
    std::shared_ptr<const FlameMasterKinetics> equation)
    : equation_{std::move(equation)} {
  VariableDatabase& vardb = *SAMRAI::hier::VariableDatabase::getDatabase();
  std::shared_ptr<SAMRAI::hier::VariableContext> muscl_context =
      vardb.getContext("MusclHancock");
  std::transform(
      equation_->GetPatchDataIds().begin(), equation_->GetPatchDataIds().end(),
      half_time_ids.begin(), [&](int id) {
        std::shared_ptr<SAMRAI::hier::Variable> variable;
        vardb.mapIndexToVariable(id, variable);
        return vardb.registerVariableAndContext(
            variable, muscl_context,
            SAMRAI::hier::IntVector::getOne(equation_->GetDimension()));
      });
}

SAMRAI::pdat::CellData<double>&
MusclHancockMethod::GetHalfTimeCellData(const SAMRAI::hier::Patch& patch,
                                        Variable variable) const {
  return *static_cast<SAMRAI::pdat::CellData<double>*>(
      patch.getPatchData(half_time_ids_[int(variable)]).get());
}

IdealGasEquation::CompletePatchData MusclHancockMethod::GetHalfTimePatchData(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = GetHalfTimeCellData(patch, Variable::density);
  SAMRAI::pdat::CellData<double>& momentum = GetHalfTimeCellData(patch, Variable::momentum);
  SAMRAI::pdat::CellData<double>& energy = GetHalfTimeCellData(patch, Variable::energy);
  SAMRAI::pdat::CellData<double>& pressure = GetHalfTimeCellData(patch, Variable::pressure);
  SAMRAI::pdat::CellData<double>& temperature = GetHalfTimeCellData(patch, Variable::temperature);
  SAMRAI::pdat::CellData<double>& speed_of_sound = GetHalfTimeCellData(patch, Variable::speed_of_sound);
  SAMRAI::pdat::CellData<double>& species = GetHalfTimeCellData(patch, Variable::species);
  // clang-format on
  return {density,     momentum,       energy, pressure,
          temperature, speed_of_sound, species};
}

namespace {
void ReconstructStateAtHalfTime_(
    const IdealGasEquation::CompletePatchData& to,
    const IdealGasEquation::CompletePatchData& from, double time, Direction dir,
    const FlameMasterKinetics& kinetics) {
      
    }
}

void ComputeFluxesOnPatch(const FluxPatchData& fluxes,
                          const CompletePatchData& state,
                          const SAMRAI::hier::Patch& patch, double dt,
                          Direction dir) {
  const CompletePatchData state_at_half_time = GetHalfTimePatchData(patch);
  ReconstructStateAtHalfTime_(statestate_at_half_time, state, dt / 2, dir,
                              *equation_);
  return HlleGodunovMethod::ComputeFluxesOnPatch(fluxes, state_at_half_time,
                                                 patch, dt / 2, dir);
}

} // namespace ideal_gas
} // namespace fub