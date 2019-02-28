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

#ifndef FUB_SOLVER_DIMENSIONAL_SPLIT_HYPERBOLIC_TIME_INTEGRATOR_HPP
#define FUB_SOLVER_DIMENSIONAL_SPLIT_HYPERBOLIC_TIME_INTEGRATOR_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"

namespace fub {

template <typename Eq> class HyperbolicSplitPatchIntegrator {
public:
  using Equation = Eq;

  HyperbolicSplitPatchIntegrator(const Equation& eq) : equation{eq} {}

  const Equation& GetEquation() const noexcept { return equation; }

  template <typename NextView, typename FluxView, typename PrevView>
  void UpdateConservatively(NextView next, FluxView fluxes, PrevView prev,
                            Duration dt, double dx,
                            Direction dir = Direction::X);

private:
  Equation equation;

  Cons<Equation> next_state;
  Cons<Equation> prev_state;
  Cons<Equation> flux_left;
  Cons<Equation> flux_right;
};

template <typename Equation>
template <typename NextView, typename FluxView, typename PrevView>
void HyperbolicSplitPatchIntegrator<Equation>::UpdateConservatively(
    NextView next, FluxView fluxes, PrevView prev, Duration dt, double dx,
    Direction dir) {
  const double lambda = dt.count() / dx;
  FUB_ASSERT(Extents(next) == Extents(prev));
  constexpr int Rank = Equation::Rank();
  ForEachIndex(Mapping(next), [&](const auto... is) {
    std::array<std::ptrdiff_t, Rank> index{is...};
    FUB_ASSERT(Extents(fluxes).extent(int(dir)) ==
               Extents(prev).extent(int(dir)) + 1);
    // Load fluxes left and right at index
    Load(flux_left, fluxes, index);
    Load(flux_right, fluxes, Shift(index, dir, 1));
    // Load state at index
    Load(prev_state, prev, index);
    // Do The computation
    ForEachComponent(
        [lambda](auto&& next, auto prev, auto flux_left, auto flux_right) {
          next = prev + lambda * (flux_left - flux_right);
        },
        next_state, prev_state, flux_left, flux_right);
    // Store the result at index
    Store(next, next_state, index);
  });
}

} // namespace fub

#endif