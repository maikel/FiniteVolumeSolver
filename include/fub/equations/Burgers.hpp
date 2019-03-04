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

#include "fub/Direction.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/State.hpp"
#include "fub/VariableDescription.hpp"
#include "fub/ext/Eigen.hpp"
#include "fub/ext/hana.hpp"

#include <array>

namespace fub {
template <typename U> struct BurgersVariables {
  BOOST_HANA_DEFINE_STRUCT(BurgersVariables, (U, u));
};

struct Burgers1d : VariableDescription<BurgersVariables<Scalar>> {
  using Complete = ::fub::Complete<Burgers1d>;
  using Conservative = ::fub::Conservative<Burgers1d>;

  static constexpr int Rank() { return 1; }

  void Flux(Conservative& flux, const Complete& state,
            Direction dir = Direction::X) const noexcept;
};

template <> class ExactRiemannSolver<Burgers1d> {
public:
  using Complete = typename Burgers1d::Complete;

  ExactRiemannSolver(const Burgers1d&){}

  void SolveRiemannProblem(Complete& state, const Complete& left,
                           const Complete& right, Direction dir) const;

  std::array<double, 1> ComputeSignals(const Complete& left,
                                       const Complete& right,
                                       Direction dir) const;
};

} // namespace fub

#endif