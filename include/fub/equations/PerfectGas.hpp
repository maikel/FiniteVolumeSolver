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

#ifndef FUB_EQUATIONS_PERFECT_GAS_HPP
#define FUB_EQUATIONS_PERFECT_GAS_HPP

#include "fub/Equation.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/EinfeldtSignalVelocities.hpp"
#include "fub/ext/Eigen.hpp"

#include <array>

namespace fub {
template <typename Density, typename Momentum, typename Energy,
          typename Pressure, typename SpeedOfSound>
struct PerfectGasComplete {
  BOOST_HANA_DEFINE_STRUCT(PerfectGasComplete, (Density, density),
                           (Momentum, momentum), (Energy, energy),
                           (Pressure, pressure), (SpeedOfSound, speed_of_sound));
};

template <typename Density, typename Momentum, typename Energy>
struct PerfectGasConservative {
  BOOST_HANA_DEFINE_STRUCT(PerfectGasConservative, (Density, density),
                           (Momentum, momentum), (Energy, energy));
};

template <int Rank>
using PerfectGasConsShape =
    PerfectGasConservative<constant<1>, constant<Rank>, constant<1>>;

template <int Rank>
using PerfectGasCompleteShape =
    PerfectGasComplete<constant<1>, constant<Rank>, constant<1>, constant<1>,
                       constant<1>>;

template <int Rank> struct PerfectGas;

template <>
struct PerfectGas<1>
    : VariableDescription<PerfectGasConsShape<1>, PerfectGasCompleteShape<1>> {
  using Cons = ::fub::Cons<PerfectGas<1>>;
  using Complete = ::fub::Complete<PerfectGas<1>>;

  static constexpr int Rank() noexcept { return 1; }

  void Flux(Cons& flux, const Complete& state,
            Direction dir = Direction::X) const noexcept;

  void Reconstruct(Complete& complete, const Cons& cons) const noexcept;

  double gamma{1.4};
  double gamma_minus_1_inv{1.0 / (gamma - 1.0)};
};

template <>
struct PerfectGas<2>
    : VariableDescription<PerfectGasConsShape<2>, PerfectGasCompleteShape<2>> {

  using Cons = ::fub::Cons<PerfectGas<2>>;
  using Complete = ::fub::Complete<PerfectGas<2>>;

  static constexpr int Rank() noexcept { return 2; }

  void Flux(Cons& flux, const Complete& state,
            Direction dir = Direction::X) const noexcept;

  void Reconstruct(Complete& complete, const Cons& cons) const noexcept;

  double gamma{1.4};
  double gamma_minus_1_inv{1.0 / (gamma - 1.0)};
};

template <>
struct PerfectGas<3>
    : VariableDescription<PerfectGasConsShape<3>, PerfectGasCompleteShape<3>> {
  using Cons = ::fub::Cons<PerfectGas<3>>;
  using Complete = ::fub::Complete<PerfectGas<3>>;

  static constexpr int Rank() noexcept { return 3; }

  void Flux(Cons& flux, const Complete& state,
            Direction dir = Direction::X) const noexcept;

  void Reconstruct(Complete& complete, const Cons& cons) const noexcept;

  double gamma{1.4};
  double gamma_minus_1_inv{1.0 / (gamma - 1.0)};
};


static_assert(HasScalarFlux<PerfectGas<1>>::value);
static_assert(HasScalarFlux<PerfectGas<2>>::value);
static_assert(HasScalarFlux<PerfectGas<3>>::value);

static_assert(HasScalarReconstruction<PerfectGas<1>>::value);
static_assert(HasScalarReconstruction<PerfectGas<2>>::value);
static_assert(HasScalarReconstruction<PerfectGas<3>>::value);

/// \ingroup EinfeldtSignalVelocities
template <>
class EinfeldtSignalVelocities<PerfectGas<1>> {
public:
  using Complete = typename PerfectGas<1>::Complete;

  EinfeldtSignalVelocities(const PerfectGas<1>& equation)
      : equation_{equation} {}

  std::array<double, 2> ComputeSignals(const Complete& left,
                                       const Complete& right,
                                       Direction dir) const;

private:
  PerfectGas<1> equation_;
};

/// \ingroup EinfeldtSignalVelocities
template <>
class EinfeldtSignalVelocities<PerfectGas<2>> {
public:
  using Complete = typename PerfectGas<2>::Complete;

  EinfeldtSignalVelocities(const PerfectGas<2>& equation)
      : equation_{equation} {}

  std::array<double, 2> ComputeSignals(const Complete& left,
                                       const Complete& right,
                                       Direction dir) const;

private:
  PerfectGas<2> equation_;
};


/// \ingroup EinfeldtSignalVelocities
template <>
class EinfeldtSignalVelocities<PerfectGas<3>> {
public:
  using Complete = typename PerfectGas<3>::Complete;

  EinfeldtSignalVelocities(const PerfectGas<3>& equation)
      : equation_{equation} {}

  std::array<double, 2> ComputeSignals(const Complete& left,
                                       const Complete& right,
                                       Direction dir) const;

private:
  PerfectGas<3> equation_;
};


} // namespace fub

#endif