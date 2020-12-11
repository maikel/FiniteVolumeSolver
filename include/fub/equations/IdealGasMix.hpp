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

#ifndef FUB_EQUATIONS_IDEAL_GAS_MIX_HPP
#define FUB_EQUATIONS_IDEAL_GAS_MIX_HPP

#include "fub/equations/ideal_gas_mix/FlameMasterReactor.hpp"

#include "fub/CompleteFromCons.hpp"
#include "fub/Equation.hpp"
#include "fub/State.hpp"
#include "fub/ext/Eigen.hpp"

#include <array>

namespace fub {

/// This is a template class for constructing conservative states for the
/// perfect gas equations.
template <typename Density, typename Momentum, typename Energy,
          typename Species>
struct IdealGasMixConservative {
  Density density;
  Momentum momentum;
  Energy energy;
  Species species;
};

template <int Rank>
using IdealGasConservativeShape =
    IdealGasMixConservative<ScalarDepth, VectorDepth<Rank>, ScalarDepth,
                            VectorDepth<-1>>;

// We "register" the conservative state with our framework.
// This enables us to name and iterate over all member variables in a given
// conservative state.
template <typename... Xs> struct StateTraits<IdealGasMixConservative<Xs...>> {
  static constexpr auto names =
      std::make_tuple("Density", "Momentum", "Energy", "Species");

  static constexpr auto pointers_to_member =
      std::make_tuple(&IdealGasMixConservative<Xs...>::density,
                      &IdealGasMixConservative<Xs...>::momentum,
                      &IdealGasMixConservative<Xs...>::energy,
                      &IdealGasMixConservative<Xs...>::species);

  template <int Rank> using Depths = IdealGasConservativeShape<Rank>;
};

template <typename Density, typename Momentum, typename Energy,
          typename Species, typename Pressure, typename SpeedOfSound,
          typename Temperature, typename HeatCapacityAtConstantPressure,
          typename HeatCapacityRatio>
struct IdealGasMixComplete
    : IdealGasMixConservative<Density, Momentum, Energy, Species> {
  Pressure pressure;
  SpeedOfSound speed_of_sound;
  Temperature temperature;
  HeatCapacityAtConstantPressure c_p;
  HeatCapacityRatio gamma;
};

template <int Rank>
using IdealGasMixCompleteShape =
    IdealGasMixComplete<ScalarDepth, VectorDepth<Rank>, ScalarDepth,
                        VectorDepth<-1>, ScalarDepth, ScalarDepth, ScalarDepth,
                        ScalarDepth, ScalarDepth>;

// We "register" the complete state with our framework.
// This enables us to name and iterate over all member variables in a given
// compete state.
template <typename... Xs> struct StateTraits<IdealGasMixComplete<Xs...>> {
  static constexpr auto names = std::make_tuple(
      "Density", "Momentum", "Energy", "Species", "Pressure", "SpeedOfSound",
      "Temperature", "HeatCapacityAtConstantPressure", "HeatCapacityRatio");
  static constexpr auto pointers_to_member = std::make_tuple(
      &IdealGasMixComplete<Xs...>::density,
      &IdealGasMixComplete<Xs...>::momentum,
      &IdealGasMixComplete<Xs...>::energy, &IdealGasMixComplete<Xs...>::species,
      &IdealGasMixComplete<Xs...>::pressure,
      &IdealGasMixComplete<Xs...>::speed_of_sound,
      &IdealGasMixComplete<Xs...>::temperature,
      &IdealGasMixComplete<Xs...>::c_p, &IdealGasMixComplete<Xs...>::gamma);

  template <int Rank> using Depths = IdealGasMixCompleteShape<Rank>;
};

template <int N> class IdealGasMix {
public:
  using ConservativeDepths = IdealGasConservativeShape<N>;
  using CompleteDepths = IdealGasMixCompleteShape<N>;

  using Conservative = ::fub::Conservative<IdealGasMix<N>>;
  using ConservativeBase = ::fub::ConservativeBase<IdealGasMix<N>>;
  using Complete = ::fub::Complete<IdealGasMix<N>>;

  using ConservativeArray = ::fub::ConservativeArray<IdealGasMix<N>>;
  using ConservativeArrayBase = ::fub::ConservativeArrayBase<IdealGasMix<N>>;
  using CompleteArray = ::fub::CompleteArray<IdealGasMix<N>>;

  explicit IdealGasMix(const FlameMasterReactor& reactor)
      : reactor_(std::move(reactor)) {}

  static constexpr int Rank() noexcept { return N; }

  void Flux(Conservative& flux, const Complete& state,
            Direction dir = Direction::X) const noexcept;

  void Flux(ConservativeArray& flux, const CompleteArray& state,
            Direction dir = Direction::X) const noexcept;

  void Flux(ConservativeArray& flux, const CompleteArray& state, MaskArray mask,
            Direction dir) const noexcept;

  void CompleteFromCons(Complete& state, const ConservativeBase& cons);

  void CompleteFromCons(CompleteArray& state,
                        const ConservativeArrayBase& cons);

  void CompleteFromCons(CompleteArray& state, const ConservativeArrayBase& cons,
                        MaskArray mask);

  FlameMasterReactor& GetReactor() noexcept { return reactor_; }
  const FlameMasterReactor& GetReactor() const noexcept { return reactor_; }

  void SetReactorStateFromComplete(const Complete& state);

  static Array<double, N, 1> Velocity(const ConservativeBase& cons) noexcept;
  static Array<double, N> Velocity(const ConservativeArrayBase& cons) noexcept;
  static Array<double, N> Velocity(const ConservativeArrayBase& cons,
                                   MaskArray mask) noexcept;

  void CompleteFromReactor(Complete& state,
                           const Eigen::Array<double, N, 1>& velocity =
                               Eigen::Array<double, N, 1>::Zero()) const;

  void CompleteFromReactor(
      CompleteArray& state,
      const Array<double, N>& velocity = Array<double, N>::Zero()) const;

  void CompleteFromReactor(CompleteArray& state,
                           const Array<double, N>& velocit,
                           MaskArray mask) const;

private:
  FlameMasterReactor reactor_;
  Array<double, Eigen::Dynamic, 1> species_buffer_{reactor_.GetNSpecies()};

  template <typename State>
  friend constexpr auto tag_invoke(tag_t<Depths>, const IdealGasMix& equation,
                                   Type<State>) noexcept {
    ToConcreteDepths<typename State::Traits::template Depths<N>> depths{};
    depths.species = equation.reactor_.GetNSpecies();
    return depths;
  }
};

// We define this class only for dimensions 1 to 3.
// The definitions will be found in its source file IdealGasMix.cpp
extern template class IdealGasMix<1>;
extern template class IdealGasMix<2>;
extern template class IdealGasMix<3>;

template <int Dim>
double KineticEnergy(double density,
                     const Eigen::Array<double, Dim, 1>& momentum) noexcept {
  return 0.5 * momentum.matrix().squaredNorm() / density;
}

template <int Dim>
Array1d KineticEnergy(Array1d density,
                      const Eigen::Array<double, Dim, kDefaultChunkSize,
                                         Eigen::RowMajor>& momentum) noexcept {
  Array1d square = Array1d::Zero();
  for (int i = 0; i < Dim; ++i) {
    square += momentum.row(i) * momentum.row(i);
  }
  return Array1d::Constant(0.5) * square / density;
}

template <int Dim>
Array1d KineticEnergy(Array1d density,
                      const Eigen::Array<double, Dim, kDefaultChunkSize,
                                         Eigen::RowMajor>& momentum,
                      MaskArray mask) noexcept {
  Array1d square = Array1d::Zero();
  for (int i = 0; i < Dim; ++i) {
    Array1d rhou = mask.select(momentum.row(i), 0.0);
    square += rhou * rhou;
  }
  Array1d density_s = mask.select(density, 1.0);
  return Array1d::Constant(0.5) * square / density_s;
}

/// @{
/// \brief Defines how to rotate a given state of the euler equations.
///
/// This function is needed when computing the reference state in the boundary
/// flux of the cut-cell stabilizations.
void Rotate(Conservative<IdealGasMix<2>>& rotated,
            const Conservative<IdealGasMix<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation, const IdealGasMix<2>&);

void Rotate(Complete<IdealGasMix<2>>& rotated,
            const Complete<IdealGasMix<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation, const IdealGasMix<2>&);

void Rotate(Conservative<IdealGasMix<3>>& rotated,
            const Conservative<IdealGasMix<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation, const IdealGasMix<3>&);

void Rotate(Complete<IdealGasMix<3>>& rotated,
            const Complete<IdealGasMix<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation, const IdealGasMix<3>&);
/// @}

void Reflect(Complete<IdealGasMix<1>>& reflected,
             const Complete<IdealGasMix<1>>& state,
             const Eigen::Matrix<double, 1, 1>& normal,
             const IdealGasMix<1>& gas);

void Reflect(Complete<IdealGasMix<2>>& reflected,
             const Complete<IdealGasMix<2>>& state,
             const Eigen::Vector2d& normal, const IdealGasMix<2>& gas);

void Reflect(Complete<IdealGasMix<3>>& reflected,
             const Complete<IdealGasMix<3>>& state,
             const Eigen::Vector3d& normal, const IdealGasMix<3>& gas);

} // namespace fub

#endif
