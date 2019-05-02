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
#include "fub/EinfeldtSignalVelocities.hpp"
#include "fub/Equation.hpp"
#include "fub/ExactRiemannSolver.hpp"
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
};

template <typename Density, typename Momentum, typename Energy,
          typename Species, typename Pressure, typename SpeedOfSound>
struct IdealGasMixComplete
    : IdealGasMixConservative<Density, Momentum, Energy, Species> {
  Pressure pressure;
  SpeedOfSound speed_of_sound;
};

// We "register" the complete state with our framework.
// This enables us to name and iterate over all member variables in a given
// compete state.
template <typename... Xs> struct StateTraits<IdealGasMixComplete<Xs...>> {
  static constexpr auto names = std::make_tuple(
      "Density", "Momentum", "Energy", "Species", "Pressure", "SpeedOfSound");
  static constexpr auto pointers_to_member = std::make_tuple(
      &IdealGasMixComplete<Xs...>::density,
      &IdealGasMixComplete<Xs...>::momentum,
      &IdealGasMixComplete<Xs...>::energy, &IdealGasMixComplete<Xs...>::species,
      &IdealGasMixComplete<Xs...>::pressure,
      &IdealGasMixComplete<Xs...>::speed_of_sound);
};

template <int Rank>
using IdealGasConservativeShape =
    IdealGasMixConservative<ScalarDepth, VectorDepth<Rank>, ScalarDepth,
                            VectorDepth<-1>>;

template <int Rank>
using IdealGasMixCompleteShape =
    IdealGasMixComplete<ScalarDepth, VectorDepth<Rank>, ScalarDepth,
                        VectorDepth<-1>, ScalarDepth, ScalarDepth>;

template <int N> class IdealGasMix {
public:
  using ConservativeDepths = IdealGasConservativeShape<N>;
  using CompleteDepths = IdealGasMixCompleteShape<N>;

  using Conservative = ::fub::Conservative<IdealGasMix<N>>;
  using Complete = ::fub::Complete<IdealGasMix<N>>;

  explicit IdealGasMix(const FlameMasterReactor& reactor)
      : reactor_(std::move(reactor)) {}

  static constexpr int Rank() noexcept { return N; }

  void Flux(Conservative& flux, const Complete& state,
            [[maybe_unused]] Direction dir = Direction::X) const noexcept;

  FlameMasterReactor& GetReactor() noexcept { return reactor_; }
  const FlameMasterReactor& GetReactor() const noexcept { return reactor_; }

  void CompleteFromReactor(Complete& state,
                           const Eigen::Array<double, N, 1>& velocity =
                               Eigen::Array<double, N, 1>::Zero()) const;

private:
  FlameMasterReactor reactor_;
};

template <typename Depth> struct ToConcreteDepthImpl { using type = Depth; };

template <> struct ToConcreteDepthImpl<VectorDepth<-1>> { using type = int; };

template <typename Depth>
using ToConcreteDepth = typename ToConcreteDepthImpl<Depth>::type;

template <typename Depths>
using ToConcreteDepths = boost::mp11::mp_transform<ToConcreteDepth, Depths>;

template <int Dim>
struct DepthsImpl<Complete<IdealGasMix<Dim>>, IdealGasMix<Dim>> {
  constexpr ToConcreteDepths<typename IdealGasMix<Dim>::CompleteDepths>
  operator()(const IdealGasMix<Dim>& equation) const noexcept {
    ToConcreteDepths<typename IdealGasMix<Dim>::CompleteDepths> depths;
    depths.species = equation.GetReactor().GetNSpecies();
    return depths;
  }
};

template <int Dim>
struct DepthsImpl<Conservative<IdealGasMix<Dim>>, IdealGasMix<Dim>> {
  constexpr ToConcreteDepths<typename IdealGasMix<Dim>::ConservativeDepths>
  operator()(const IdealGasMix<Dim>& equation) const noexcept {
    ToConcreteDepths<typename IdealGasMix<Dim>::ConservativeDepths> depths;
    depths.species = equation.GetReactor().GetNSpecies();
    return depths;
  }
};

// We define this class only for dimensions 1 to 3.
// The definitions will be found in its source file IdealGasMix.cpp
extern template class IdealGasMix<1>;
extern template class IdealGasMix<2>;
extern template class IdealGasMix<3>;

template <int Dim>
void InitializeState(const IdealGasMix<Dim>& eq,
                     Conservative<IdealGasMix<Dim>>& cons) {
  cons.species.resize(eq.GetReactor().GetNSpecies(), 1);
}

template <int Dim>
void InitializeState(const IdealGasMix<Dim>& eq,
                     Complete<IdealGasMix<Dim>>& state) {
  state.species.resize(eq.GetReactor().GetNSpecies(), 1);
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

template <int Dim> struct CompleteFromConsImpl<IdealGasMix<Dim>> {
  static void apply(IdealGasMix<Dim>& equation,
                    Complete<IdealGasMix<Dim>>& complete,
                    const Conservative<IdealGasMix<Dim>>& cons);
};

extern template struct CompleteFromConsImpl<IdealGasMix<1>>;
extern template struct CompleteFromConsImpl<IdealGasMix<2>>;
extern template struct CompleteFromConsImpl<IdealGasMix<3>>;

template <int Dim> struct EinfeldtSignalVelocitiesImpl<IdealGasMix<Dim>> {
  using Complete = typename IdealGasMix<Dim>::Complete;

  static std::array<double, 2> apply(const IdealGasMix<Dim>& equation,
                                     const Complete& left,
                                     const Complete& right,
                                     [[maybe_unused]] Direction dir) noexcept;
};

extern template struct EinfeldtSignalVelocitiesImpl<IdealGasMix<1>>;
extern template struct EinfeldtSignalVelocitiesImpl<IdealGasMix<2>>;
extern template struct EinfeldtSignalVelocitiesImpl<IdealGasMix<3>>;

} // namespace fub

#endif
