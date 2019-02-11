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

#ifndef FUB_INITIAL_DATA_RIEMANN_PROBLEM_HPP
#define FUB_INITIAL_DATA_RIEMANN_PROBLEM_HPP

#include <array>

namespace fub {

template <typename State, typename Geometry> class RiemannProblem {
public:
  explicit RiemannProblem(const State& state1, const State& state2,
                          const Geometry& geom)
      : states_{state1, state2}, geometry_{geom} {}

  template <typename StateView, typename CoordMapping>
  void InitializeData(StateView states, const CoordMapping& x) {
    ForEachIndex(Extents(states), [&](auto index) {
      if (geometry_.ComputeDistanceTo(x(index)) < 0.0) {
        states(index) = states_[0];
      } else {
        states(index) = states_[1];
      }
    });
  }

private:
  std::array<State, 2> states_;
  Geometry geometry_;
};

} // namespace fub

#endif