// Copyright (c) 2018 Maikel Nadolski
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

/// \defgroup Euler Euler Equations
/// This module defines classes and methods to be used to concrete solver
/// implementations wrt. the euler equations.

#ifndef FUB_IDEAL_GAS_IDEAL_GAS_EQUATION_HPP
#define FUB_IDEAL_GAS_IDEAL_GAS_EQUATION_HPP

#include "fub/core/assert.hpp"
#include "fub/core/mdspan.hpp"

#include "fub/Direction.hpp"
#include "fub/Eigen.hpp"
#include "fub/StateFacade.hpp"

namespace fub {
namespace ideal_gas {

struct IdealGasEquation {
  enum class Variable : int {
    density,
    momentum,
    energy,
    pressure,
    temperature,
    speed_of_sound,
    species,
    size ///< This is only a dummy value to get the number of variables
  };

  template <typename Density, typename Momentum, typename Energy,
            typename Pressure, typename Temperature, typename SpeedOfSound,
            typename Species>
  struct Complete {
    BOOST_HANA_DEFINE_STRUCT(Complete, (Density, density), (Momentum, momentum),
                             (Energy, energy), (Pressure, pressure),
                             (Temperature, temperature),
                             (SpeedOfSound, speed_of_sound),
                             (Species, species));
  };

  template <typename... Vars>
  using CompleteState = StateFacade<Complete<Vars...>>;

  template <typename Density, typename Momentum, typename Energy,
            typename Species>
  struct Conservative {
    BOOST_HANA_DEFINE_STRUCT(Conservative, (Density, density),
                             (Momentum, momentum), (Energy, energy),
                             (Species, species));
  };

  template <typename... Vars>
  using ConservativeState = StateFacade<Conservative<Vars...>>;

  template <typename Temperature, typename Velocity, typename Pressure,
            typename Species>
  struct Primitive {
    BOOST_HANA_DEFINE_STRUCT(Primitive, (Temperature, temperature),
                             (Velocity, velocity), (Pressure, pressure),
                             (Species, species));
  };

  template <typename... Vars>
  using PrimitiveState = StateFacade<Primitive<Vars...>>;

  using ConsArray = ConservativeState<Array1d, Array3d, Array1d, ArrayXd>;

  using StateArray = CompleteState<Array1d, Array3d, Array1d, Array1d, Array1d,
                                   Array1d, ArrayXd>;

  using PrimArray = PrimitiveState<Array1d, Array3d, Array1d, ArrayXd>;

  template <typename T>
  using ConsSpan =
      ConservativeState<basic_mdspan<T, DynamicExtents<3>, layout_left>,
                        basic_mdspan<T, DynamicExtents<4>, layout_left>,
                        basic_mdspan<T, DynamicExtents<3>, layout_left>,
                        basic_mdspan<T, DynamicExtents<4>, layout_left>>;

  template <typename T>
  using StateSpan =
      CompleteState<basic_mdspan<T, DynamicExtents<3>, layout_left>,
                    basic_mdspan<T, DynamicExtents<4>, layout_left>,
                    basic_mdspan<T, DynamicExtents<3>, layout_left>,
                    basic_mdspan<T, DynamicExtents<3>, layout_left>,
                    basic_mdspan<T, DynamicExtents<3>, layout_left>,
                    basic_mdspan<T, DynamicExtents<3>, layout_left>,
                    basic_mdspan<T, DynamicExtents<4>, layout_left>>;

  static void Flux(ConsArray& flux, const StateArray& state) noexcept;

  // static Array1d Flux(const StateArray& state, const Array1d& velocity,
  //                     Variable var, int n = 0) noexcept {
  //   switch (var) {
  //   case Variable::density:
  //     return state.momentum.col(0);
  //   case Variable::momentum:
  //     if (n == 0) {
  //       return state.momentum.col(0) * velocity + state.pressure;
  //     } else {
  //       return state.momentum.col(n) * velocity;
  //     }
  //   case Variable::energy:
  //     return velocity * (state.energy + state.pressure);
  //   default:
  //     return velocity * state.species.col(n);
  //   }
  // }
};

template <typename... Ts>
auto AsCons(const StateFacade<IdealGasEquation::Complete<Ts...>>& state) {
  return boost::hana::make<StateTag<IdealGasEquation::Conservative>>(
    Eigen::Map<const Array1d>(state.density.data()), 
    Eigen::Map<const Array3d>(state.momentum.data()),
    Eigen::Map<const Array1d>(state.energy.data()), 
    Eigen::Map<const ArrayXd>(state.species.data(), state.species.rows(), state.species.cols())
  );
}

} // namespace ideal_gas
} // namespace fub

#endif