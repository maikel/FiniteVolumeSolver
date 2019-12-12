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

#include "fub/ideal_gas/IdealGasEquation.hpp"

namespace fub {
namespace ideal_gas {

void IdealGasEquation::Flux(ConservativeArray& flux,
                            const StateArray& state) noexcept {
  flux.density = state.momentum.col(0);
  const Array1d velocity_x = state.momentum.col(0) / state.density;
  flux.momentum.col(0) = velocity_x * state.momentum(0) + state.pressure;
  flux.momentum.col(1) = velocity_x * state.momentum(1);
  flux.momentum.col(2) = velocity_x * state.momentum(2);
  flux.energy = velocity_x * (state.energy + state.pressure);
  for (int s = 0; s < state.species.cols(); ++s) {
    flux.species.col(s) = velocity_x * state.species.col(s);
  }
}

} // namespace ideal_gas
} // namespace fub