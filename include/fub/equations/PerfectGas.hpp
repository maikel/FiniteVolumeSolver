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

#include "fub/CompleteFromCons.hpp"
#include "fub/EinfeldtSignalVelocities.hpp"
#include "fub/Equation.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/ext/Eigen.hpp"

#include <array>

namespace fub {
template <typename Density, typename Momentum, typename Energy,
          typename Pressure, typename SpeedOfSound>
struct PerfectGasComplete {
  BOOST_HANA_DEFINE_STRUCT(PerfectGasComplete, (Density, density),
                           (Momentum, momentum), (Energy, energy),
                           (Pressure, pressure),
                           (SpeedOfSound, speed_of_sound));
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

template <int N>
struct PerfectGas
    : VariableDescription<PerfectGasConsShape<N>, PerfectGasCompleteShape<N>> {
  using Conservative = ::fub::Conservative<PerfectGas<N>>;
  using Complete = ::fub::Complete<PerfectGas<N>>;

  static constexpr int Rank() noexcept { return N; }

  void Flux(Conservative& flux, const Complete& state,
            [[maybe_unused]] Direction dir = Direction::X) const noexcept;

  double gamma{1.4};
  double gamma_minus_1_inv{1.0 / (gamma - 1.0)};
};

extern template struct PerfectGas<1>;
extern template struct PerfectGas<2>;
extern template struct PerfectGas<3>;

template <int Dim> struct CompleteFromConsImpl<PerfectGas<Dim>> {
  static void apply(const PerfectGas<Dim>& equation,
                    Complete<PerfectGas<Dim>>& complete,
                    const Conservative<PerfectGas<Dim>>& cons);

  static void apply(const PerfectGas<Dim>& equation,
                    Complete<PerfectGas<Dim>>& complete,
                    const Complete<PerfectGas<Dim>>& cons);
};

extern template struct CompleteFromConsImpl<PerfectGas<1>>;
extern template struct CompleteFromConsImpl<PerfectGas<2>>;
extern template struct CompleteFromConsImpl<PerfectGas<3>>;

template <int Dim> struct EinfeldtSignalVelocitiesImpl<PerfectGas<Dim>> {
  using Complete = typename PerfectGas<Dim>::Complete;

  static std::array<double, 2> apply(const PerfectGas<Dim>& equation,
                                     const Complete& left,
                                     const Complete& right,
                                     [[maybe_unused]] Direction dir) noexcept;
};

extern template struct EinfeldtSignalVelocitiesImpl<PerfectGas<1>>;
extern template struct EinfeldtSignalVelocitiesImpl<PerfectGas<2>>;
extern template struct EinfeldtSignalVelocitiesImpl<PerfectGas<3>>;

} // namespace fub

#endif