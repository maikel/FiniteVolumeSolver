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

#include "fub/Equation.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/State.hpp"
#include "fub/ext/Eigen.hpp"

#include <array>

namespace fub {
template <typename Heigth, typename Momentum> struct ShallowWaterVariables {
  BOOST_HANA_DEFINE_STRUCT(ShallowWaterVariables, (Heigth, heigth),
                           (Momentum, momentum));
};

struct ShallowWater
    : VariableDescription<ShallowWaterVariables<Scalar, Vector2d>> {
  using Conservative = ::fub::Conservative<ShallowWater>;
  using Complete = ::fub::Complete<ShallowWater>;
  template <int N>
  using CompleteArray = ::fub::CompleteArray<ShallowWater, N>;
  template <int N>
  using ConservativeArray = ::fub::ConservativeArray<ShallowWater, N>;
  static constexpr int Rank() noexcept { return 2; }

  void Flux(Conservative& flux, const Complete& state,
            Direction dir = Direction::X) const noexcept;

  double gravity_{10.0};
};

template <> class ExactRiemannSolver<ShallowWater> {
public:
  using Complete = typename ShallowWater::Complete;

  ExactRiemannSolver(const ShallowWater& equation) : equation_{equation} {}

  /// Returns either left or right, depending on the upwind velocity.
  void SolveRiemannProblem(Complete& state, const Complete& left,
                           const Complete& right, Direction dir);

  /// Returns the upwind velocity in the specified direction.
  std::array<double, 2> ComputeSignals(const Complete&, const Complete&,
                                       Direction dir);

  Complete ComputeMiddleState(const Complete& left, const Complete& right,
                              Direction dir);

private:
  ShallowWater equation_;
};

struct ShallowWaterSignalVelocities {
  using Complete = typename ShallowWater::Complete;

  std::array<double, 2> operator()(const ShallowWater& equation,
                                   const Complete& left, const Complete& right,
                                   Direction dir);
};

} // namespace fub

#endif