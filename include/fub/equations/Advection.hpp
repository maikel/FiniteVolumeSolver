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

#ifndef FUB_EQUATIONS_ADVECTION_HPP
#define FUB_EQUATIONS_ADVECTION_HPP

#include "fub/Eigen.hpp"
#include "fub/Equation.hpp"
#include "fub/StateFacade.hpp"

#include <array>

namespace fub {

template <int Rank> struct Advection {
  static constexpr int kRank = Rank;

  template <typename Mass> struct ConsData {
    using Equation = Advection;
    BOOST_HANA_DEFINE_STRUCT(ConsData, (Mass, mass));
  };

  template <typename Mass> struct StateData {
    using Equation = Advection;
    BOOST_HANA_DEFINE_STRUCT(StateData, (Mass, mass));
  };

  // clang-format off
  using State = StateFacade<StateData<
    double  // Mass
  >>;
  using Cons = StateFacade<ConsData<
    double  // Mass
  >>;
  // clang-format on

  using SimdState = SimdifyT<State>;
  using SimdCons = SimdifyT<Cons>;

  constexpr static int rank() { return Rank; }

  /// Initializes this equation with a fixed velocity for each direction
  /// specified by vel.
  explicit Advection(const std::array<double, Rank>& vel) : velocity{vel} {}

  /// Computes the linear transport flux in the specified direction.
  ///
  /// The output parameter will be filled with `v_{text{dir}} \cdot q`.
  ///
  /// \param[out] flux The conservative state which will store the results.
  /// \param[in] state The input state.
  /// \param[in] dir   The split direction of this flux .
  ///
  /// \note This method is mandatory.
  void Flux(SimdCons& flux, const SimdState& state,
            Direction dir = Direction::X) const noexcept;

  /// This method is just the identity.
  /// The types State and Cons differ in general and need this kind of
  /// conversion method.
  ///
  /// \note This method is mandatory.
  void Reconstruct(SimdState& state, const SimdCons& cons) const noexcept;

  /// Computes the exact solution to the 1-dim riemann problem along the
  /// relative coordinates `(0, t)`
  ///
  /// \note Providing this method is OPTIONAL but enables the automatic usage
  /// with the GodunovMethod and other utility classes.
  void SolveRiemannProblem(SimdState& state, const SimdState& left,
                           const SimdState& right,
                           Direction dir = Direction::X);

  std::array<double, Rank> velocity;
};

template <int Rank>
void Advection<Rank>::Flux(SimdCons& flux, const SimdState& state,
                           Direction dir) const noexcept {
  const int dir_v = int(dir);
  flux.mass = state.mass * velocity[dir_v];
}

template <int Rank>
void Advection<Rank>::Reconstruct(SimdState& state, const SimdCons& cons) const
    noexcept {
  state.mass = cons.mass;
}

template <int Rank>
void Advection<Rank>::SolveRiemannProblem(SimdState& state,
                                          const SimdState& left,
                                          const SimdState& right,
                                          Direction dir) {
  const int dir_v = int(dir);
  if (0.0 < velocity[dir_v]) {
    state = left;
  } else {
    state = right;
  }
}

template <std::size_t Rank>
Advection(const std::array<double, Rank>&)->Advection<Rank>;

} // namespace fub

#endif