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

#ifndef FUB_EQUATIONS_SHALLOW_WATER_HPP
#define FUB_EQUATIONS_SHALLOW_WATER_HPP

#include "fub/Eigen.hpp"
#include "fub/Equation.hpp"
#include "fub/StateFacade.hpp"

#include <array>

namespace fub {

struct ShallowWater {
  template <typename Heigth, typename Momentum> struct ConsData {
    using Equation = ShallowWater;
    BOOST_HANA_DEFINE_STRUCT(ConsData, (Heigth, heigth), (Momentum, momentum));
  };

  template <typename Heigth, typename Momentum> struct StateData {
    using Equation = ShallowWater;
    BOOST_HANA_DEFINE_STRUCT(StateData, (Heigth, heigth), (Momentum, momentum));
  };

  // clang-format off
  using State = StateFacade<StateData<
    Scalar,  // Heigth
    Array2d  // Momentum
  >>;
  using Cons = StateFacade<ConsData<
    Scalar,  // Heigth
    Array2d  // Momentum
  >>;
  // clang-format on

  /// Computes the linear transport flux in the specified direction.
  ///
  /// \note This method is mandatory.
  void Flux(Cons& flux, const State& state, Direction dir = Direction::X) const
      noexcept;

  /// This method is just the identity.
  /// The types State and Cons differ in general and need this kind of
  /// conversion method.
  ///
  /// \note This method is mandatory.
  void Reconstruct(State& state, const Cons& cons) const noexcept;

  /// Computes the exact solution to the 1-dim riemann problem at relative
  /// coordinates `(x, t)`
  ///
  /// \note Providing this method is OPTIONAL but enables the automatic usage
  /// with the GodunovMethod and other utility classes.
  void SolveRiemannProblem(State& state, const State& left, const State& right,
                           Direction dir = Direction::X);

  Array2d ComputeSignals(const State& left, const State& right, Direction dir = Direction::X);

  Scalar gravity_{Scalar::Constant(9.81)};
};

template <typename T, typename S>
auto AsCons(const StateFacade<ShallowWater::StateData<T, S>>& state) {
  return boost::hana::make<StateTag<ShallowWater::ConsData>>(
      Scalar::Map(state.heigth.data()), Array2d::Map(state.momentum.data()));
}

template <typename T, typename S>
auto AsCons(StateFacade<ShallowWater::StateData<T, S>>& state) {
  return boost::hana::make<StateTag<ShallowWater::ConsData>>(
      Scalar::Map(state.heigth.data()), Array2d::Map(state.momentum.data()));
}

} // namespace fub

#endif