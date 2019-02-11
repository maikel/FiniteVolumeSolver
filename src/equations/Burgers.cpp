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

#include "fub/equations/Burgers.hpp"

namespace fub {

void Burgers::Flux(Cons& flux, const State& state, Direction dir) const
    noexcept {
  flux.mass = 0.5 * state.mass * state.mass;
}

void Burgers::Reconstruct(State& state, const Cons& cons) const noexcept {
  state.mass = cons.mass;
}

void Burgers::SolveRiemannProblem(State& riemann_solution, const State& left,
                                  const State& right, Direction) {
  const Array1d s = (0.5 * (right.mass * right.mass - left.mass * left.mass)) /
                    (right.mass - left.mass);
  Array1d zero;
  zero.fill(0.0);
  const auto sol1 = (s < zero).select(left.mass, right.mass);
  const auto sol2 = (right.mass < zero)
          .select(right.mass, (left.mass > zero).select(left.mass, zero));
  riemann_solution.mass = (left.mass < right.mass).select(sol1, sol2);
}

} // namespace fub