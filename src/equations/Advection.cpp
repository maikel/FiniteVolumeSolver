// Copyright (c) 2019 Maikel Nadolski
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

#include "fub/equations/Advection.hpp"

namespace fub {

void Advection::Flux(Cons& flux, const State& state, Direction dir) const
    noexcept {
  const int dir_v = int(dir);
  flux.mass = state.mass * velocity[dir_v];
}

void Advection::Reconstruct(State& state, const Cons& cons) const noexcept {
  state.mass = cons.mass;
}

void Advection::SolveRiemannProblem(State& state, const State& left,
                                    const State& right,
                                    Direction dir) {
  const int dir_v = int(dir);
  if (0.0 < velocity[dir_v]) {
    state = left;
  } else {
    state = right;
  }
}

} // namespace fub