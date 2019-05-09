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

// #include "fub/Equation.hpp"
#include "fub/Direction.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/State.hpp"
#include "fub/StateArray.hpp"
#include "fub/ext/Eigen.hpp"

#include <array>

namespace fub {
template <typename Mass> struct AdvectionVariables { Mass mass; };

template <typename... Xs> struct StateTraits<AdvectionVariables<Xs...>> {
  static constexpr auto names = std::make_tuple("Mass");
  static constexpr auto pointers_to_member =
      std::make_tuple(&AdvectionVariables<Xs...>::mass);
};

struct Advection2d {
  using ConservativeDepths = AdvectionVariables<ScalarDepth>;
  using CompleteDepths = ConservativeDepths;

  using Complete = ::fub::Complete<Advection2d>;
  using Conservative = ::fub::Conservative<Advection2d>;

  template <int N> using CompleteArray = ::fub::CompleteArray<Advection2d, N>;
  template <int N>
  using ConservativeArray = ::fub::ConservativeArray<Advection2d, N>;

  /// Constructs an equation object with velocity `v`.
  Advection2d(const std::array<double, 2>& v) noexcept : velocity{v} {}

  /// Returns 2, which is the space dimension for this equation.
  static constexpr int Rank() { return 2; }

  /// Computes the linear transport flux in the specified direction.
  ///
  /// \param[out] flux The conservative state which will store the results.
  /// \param[in] state The input state.
  /// \param[in] dir   The split direction of this flux.
  void Flux(Conservative& flux, const Complete& state, Direction dir) const
      noexcept;

  /// Computes the linear transport flux in the specified direction (Array
  /// version).
  ///
  /// \param[out] flux The conservative state which will store the results.
  /// \param[in] state The input state.
  /// \param[in] dir   The split direction of this flux.
  template <int N>
  void Flux(ConservativeArray<N>& flux, const CompleteArray<N>& state,
            Direction dir) const noexcept {
    const int d = static_cast<int>(dir);
    flux.mass = velocity[d] * state.mass;
  }

  /// This member variable stores the constant transport velocity.
  std::array<double, 2> velocity;
};

template <> class ExactRiemannSolver<Advection2d> {
public:
  using Complete = typename Advection2d::Complete;
  template <int N> using CompleteArray = typename Advection2d::CompleteArray<N>;

  ExactRiemannSolver(const Advection2d& equation) : equation_{equation} {}

  /// Returns either left or right, depending on the upwind velocity.
  void SolveRiemannProblem(Complete& state, const Complete& left,
                           const Complete& right, Direction dir);

  /// Returns either left or right, depending on the upwind velocity.
  template <int N>
  void SolveRiemannProblem(CompleteArray<N>& state,
                           const CompleteArray<N>& left,
                           const CompleteArray<N>& right, Direction dir) {
    if (equation_.velocity[static_cast<std::size_t>(dir)] > 0) {
      state = left;
    } else {
      state = right;
    }
  }

  /// Returns the upwind velocity in the specified direction.
  std::array<double, 1> ComputeSignals(const Complete&, const Complete&,
                                       Direction dir);

  /// Returns the upwind velocity in the specified direction.
  template <int N>
  std::array<Array<double, 1, N>, 1> ComputeSignals(const CompleteArray<N>&,
                                                    const CompleteArray<N>&,
                                                    Direction dir) {
    return {Array<double, 1, N>::Constant(
        equation_.velocity[static_cast<std::size_t>(dir)])};
  }

private:
  Advection2d equation_;
};

} // namespace fub

#endif
