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

#ifndef FUB_EQUATIONS_BURGERS_HPP
#define FUB_EQUATIONS_BURGERS_HPP

#include "fub/Eigen.hpp"
#include "fub/StateFacade.hpp"
#include "fub/Equation.hpp"

#include <array>

namespace fub {

struct Burgers {
  template <typename Mass> struct ConsData {
    using Equation = Burgers;
    BOOST_HANA_DEFINE_STRUCT(ConsData, (Mass, mass));
  };

  template <typename Mass> struct StateData {
    using Equation = Burgers;
    BOOST_HANA_DEFINE_STRUCT(StateData, (Mass, mass));
  };

  // clang-format off
  using State = StateFacade<StateData<
    Scalar  // Mass
  >>;
  using Cons = StateFacade<ConsData<
    Scalar  // Mass
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
  /// coordinates `(0, t)`
  ///
  /// \note Providing this method is OPTIONAL but enables the automatic usage
  /// with the GodunovMethod and other utility classes.
  void SolveRiemannProblem(State& state, const State& left, const State& right,
                           Direction dir = Direction::X);
};

template <typename T>
auto AsCons(const StateFacade<Burgers::StateData<T>>& state) {
  return boost::hana::make<StateTag<Burgers::ConsData>>(
      Scalar::Map(state.mass.data()));
}

template <typename T>
auto AsCons(StateFacade<Burgers::StateData<T>>& state) {
  return boost::hana::make<StateTag<Burgers::ConsData>>(
      Scalar::Map(state.mass.data()));
}

} // namespace fub

#endif