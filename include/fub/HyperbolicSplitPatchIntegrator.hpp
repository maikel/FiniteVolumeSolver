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
#include "fub/ForEach.hpp"

namespace fub {

template <typename Eq> class HyperbolicSplitPatchIntegrator {
public:
  using Equation = Eq;

  HyperbolicSplitPatchIntegrator(const Equation& eq) : equation{eq} {}

  const Equation& GetEquation() const noexcept { return equation; }

  template <typename NextView, typename FluxView, typename PrevView>
  void UpdateConservatively(NextView next, FluxView fluxes, PrevView prev, Direction dir,
                            Duration dt, double dx);

private:
  Equation equation;

  Conservative<Equation> next_state{equation};
  Conservative<Equation> prev_state{equation};
  Conservative<Equation> flux_left{equation};
  Conservative<Equation> flux_right{equation};
};

template <typename Equation>
template <typename NextView, typename FluxView, typename PrevView>
void HyperbolicSplitPatchIntegrator<Equation>::UpdateConservatively(
    NextView next, FluxView fluxes, PrevView prev, Direction dir, Duration dt, double dx) {
  const double lambda = dt.count() / dx;
  constexpr int Rank = Equation::Rank();
  ForEachIndex(Shrink(Box<0>(fluxes), dir, {0, 1}), [&](const auto... is) {
    const std::array<std::ptrdiff_t, Rank> face_left{is...};
    const std::array<std::ptrdiff_t, Rank> face_right = Shift(face_left, dir, 1);
    // Load fluxes left and right at index
    Load(flux_left, fluxes, face_left);
    Load(flux_right, fluxes, face_right);
    // Load state at index
    const std::array<std::ptrdiff_t, Rank> cell = face_left;
    Load(prev_state, prev, cell);
    // Do The computation
    ForEachVariable<Conservative<Equation>>(
        [lambda](auto&& next, auto prev, auto flux_left, auto flux_right) {
          next = prev + lambda * (flux_left - flux_right);
        },
        next_state, prev_state, flux_left, flux_right);
    // Store the result at index
    Store(next, next_state, cell);
  });
}

} // namespace fub

#endif
