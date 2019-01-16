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

#ifndef FUB_IDEAL_GAS_INITIAL_DATA_AMBIENT_HPP
#define FUB_IDEAL_GAS_INITIAL_DATA_AMBIENT_HPP

#include "fub/SAMRAI/ideal_gas/IdealGasEquation.hpp"
#include <vector>

namespace fub {
namespace ideal_gas {

class Ambient : public fub::InitialCondition {
public:
  using CompletePatchData = IdealGasEquation::CompletePatchData;

  struct PrimState {
    double temperature;
    double pressure;
    double velocity;
    std::vector<double> species;
  };

  const IdealGasEquation* ideal_gas;
  PrimState state;

  Ambient(PrimState w, const IdealGasEquation& equation)
      : ideal_gas{&equation}, state(std::move(w)) {}

private:
  void FillPrimState(const SAMRAI::hier::Patch& patch,
                     const SAMRAI::hier::Index& index) const {
    IdealGasEquation::PrimPatchData prim =
        ideal_gas->GetPrimPatchData(patch);
    SAMRAI::pdat::CellIndex i(index);
    prim.momentum(i) = state.velocity;
    prim.pressure(i) = state.pressure;
    prim.temperature(i) = state.temperature;
    for (int s = 0; s < state.species.size(); ++s) {
      prim.species(i, s) = state.species[s];
    }
  }

  void InitializeDataOnPatch(const SAMRAI::hier::Patch& patch) const override {
    const SAMRAI::hier::Box& box = patch.getBox();
    for (const SAMRAI::hier::Index& index : box) {
      FillPrimState(patch, index);
    }
    ideal_gas->FillFromPrim(ideal_gas->GetCompletePatchData(patch),
                            ideal_gas->GetPrimPatchData(patch));
  }
};

} // namespace ideal_gas
} // namespace fub

#endif