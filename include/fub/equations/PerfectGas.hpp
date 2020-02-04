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

#include "fub/State.hpp"
#include "fub/StateArray.hpp"

#include "fub/CompleteFromCons.hpp"
#include "fub/Equation.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/ext/Eigen.hpp"

#include "fub/EinfeldtSignalVelocities.hpp"
#include "fub/flux_method/HllMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"

#include <array>

namespace fub {

/// This is a template class for constructing conservative states for the
/// perfect gas equations.
template <typename Density, typename Momentum, typename Energy>
struct PerfectGasConservative {
  Density density;
  Momentum momentum;
  Energy energy;
};

// We "register" the conservative state with our framework.
// This enables us to name and iterate over all member variables in a given
// conservative state.
template <typename... Xs> struct StateTraits<PerfectGasConservative<Xs...>> {
  static constexpr auto names =
      std::make_tuple("Density", "Momentum", "Energy");

  static constexpr auto pointers_to_member =
      std::make_tuple(&PerfectGasConservative<Xs...>::density,
                      &PerfectGasConservative<Xs...>::momentum,
                      &PerfectGasConservative<Xs...>::energy);
};

template <typename Density, typename Momentum, typename Energy,
          typename Pressure, typename SpeedOfSound>
struct PerfectGasComplete : PerfectGasConservative<Density, Momentum, Energy> {
  Pressure pressure;
  SpeedOfSound speed_of_sound;
};

// We "register" the complete state with our framework.
// This enables us to name and iterate over all member variables in a given
// compete state.
template <typename... Xs> struct StateTraits<PerfectGasComplete<Xs...>> {
  static constexpr auto names = std::make_tuple("Density", "Momentum", "Energy",
                                                "Pressure", "SpeedOfSound");
  static constexpr auto pointers_to_member = std::make_tuple(
      &PerfectGasComplete<Xs...>::density, &PerfectGasComplete<Xs...>::momentum,
      &PerfectGasComplete<Xs...>::energy, &PerfectGasComplete<Xs...>::pressure,
      &PerfectGasComplete<Xs...>::speed_of_sound);
};

template <int Rank>
using PerfectGasConsShape =
    PerfectGasConservative<ScalarDepth, VectorDepth<Rank>, ScalarDepth>;

template <int Rank>
using PerfectGasCompleteShape =
    PerfectGasComplete<ScalarDepth, VectorDepth<Rank>, ScalarDepth, ScalarDepth,
                       ScalarDepth>;

template <int Rank> struct PerfectGas;

template <int N> struct PerfectGas {
  using ConservativeDepths = PerfectGasConsShape<N>;
  using CompleteDepths = PerfectGasCompleteShape<N>;

  using Conservative = ::fub::Conservative<PerfectGas<N>>;
  using Complete = ::fub::Complete<PerfectGas<N>>;
  using ConservativeArray = ::fub::ConservativeArray<PerfectGas<N>>;
  using CompleteArray = ::fub::CompleteArray<PerfectGas<N>>;

  static constexpr int Rank() noexcept { return N; }

  void Flux(Conservative& flux, const Complete& state,
            [[maybe_unused]] Direction dir = Direction::X) const noexcept;

  void Flux(ConservativeArray& flux, const CompleteArray& state,
            [[maybe_unused]] Direction dir) const noexcept;

  void Flux(ConservativeArray& flux, const CompleteArray& state, MaskArray mask,
            [[maybe_unused]] Direction dir) const noexcept;

  void CompleteFromCons(Complete& complete,
                        const ConservativeBase<PerfectGas>& cons) const
      noexcept;

  void CompleteFromCons(CompleteArray& complete,
                        const ConservativeArrayBase<PerfectGas>& cons) const
      noexcept;

  void CompleteFromCons(CompleteArray& complete,
                        const ConservativeArrayBase<PerfectGas>& cons,
                        MaskArray mask) const noexcept;

  Complete CompleteFromPrim(double density, const Array<double, N, 1>& u,
                            double pressure) const noexcept;
  Array<double, N, 1> Velocity(const Complete& q) const noexcept;
  double Machnumber(const Complete& q) const noexcept;

  double Rspec{287.058};

  double gamma{1.4};
  double gamma_minus_1_inv{1.0 / (gamma - 1.0)};

  Array1d gamma_array_{Array1d::Constant(gamma)};
  Array1d gamma_minus_1_inv_array_{Array1d::Constant(gamma_minus_1_inv)};
};

// We define this class only for dimensions 1 to 3.
// The definitions will be found in its source file PerfetGas.cpp
extern template struct PerfectGas<1>;
extern template struct PerfectGas<2>;
extern template struct PerfectGas<3>;

/// @{
/// \brief Defines how to rotate a given state of the euler equations.
///
/// This function is needed when computing the reference state in the boundary
/// flux of the cut-cell stabilizations.
void Rotate(Conservative<PerfectGas<2>>& rotated,
            const Conservative<PerfectGas<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation, const PerfectGas<2>&);

void Rotate(Complete<PerfectGas<2>>& rotated,
            const Complete<PerfectGas<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation, const PerfectGas<2>&);

void Rotate(Conservative<PerfectGas<3>>& rotated,
            const Conservative<PerfectGas<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation, const PerfectGas<3>&);

void Rotate(Complete<PerfectGas<3>>& rotated,
            const Complete<PerfectGas<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation, const PerfectGas<3>&);
/// @}

void Reflect(Complete<PerfectGas<1>>& reflected,
             const Complete<PerfectGas<1>>& state,
             const Eigen::Matrix<double, 1, 1>& normal,
             const PerfectGas<1>& gas);

void Reflect(Complete<PerfectGas<2>>& reflected,
             const Complete<PerfectGas<2>>& state,
             const Eigen::Vector2d& normal, const PerfectGas<2>& gas);

void Reflect(Complete<PerfectGas<3>>& reflected,
             const Complete<PerfectGas<3>>& state,
             const Eigen::Vector3d& normal, const PerfectGas<3>& gas);

template <int Dim> class ExactRiemannSolver<PerfectGas<Dim>> {
public:
  using Complete = typename PerfectGas<Dim>::Complete;
  using CompleteArray = typename PerfectGas<Dim>::CompleteArray;

  explicit ExactRiemannSolver(const PerfectGas<Dim>& equation)
      : equation_{equation} {}

  /// Returns either left or right, depending on the upwind velocity.
  void SolveRiemannProblem(Complete& state, const Complete& left,
                           const Complete& right, Direction dir);

  void SolveRiemannProblem(CompleteArray& state, const CompleteArray& left,
                           const CompleteArray& right, Direction dir);

  void SolveRiemannProblem(CompleteArray& state, const CompleteArray& left,
                           const CompleteArray& right, MaskArray mask,
                           Direction dir);

  /// Returns the upwind velocity in the specified direction.
  std::array<double, 2> ComputeSignals(const Complete&, const Complete&,
                                       Direction dir);

  std::array<Array1d, 2> ComputeSignals(const CompleteArray&,
                                        const CompleteArray&, Direction dir);

  std::array<double, 2> ComputeMiddleState(const Complete& left,
                                           const Complete& right,
                                           Direction dir);

private:
  PerfectGas<Dim> equation_;
};

extern template class ExactRiemannSolver<PerfectGas<1>>;
extern template class ExactRiemannSolver<PerfectGas<2>>;
extern template class ExactRiemannSolver<PerfectGas<3>>;

template <int Dim> struct EinfeldtSignalVelocities<PerfectGas<Dim>> {
  using Complete = ::fub::Complete<PerfectGas<Dim>>;
  using CompleteArray = ::fub::CompleteArray<PerfectGas<Dim>>;

  std::array<double, 2> operator()(const PerfectGas<Dim>& equation,
                                   const Complete& left, const Complete& right,
                                   Direction dir) const noexcept;

  std::array<Array1d, 2> operator()(const PerfectGas<Dim>& equation,
                                    const CompleteArray& left,
                                    const CompleteArray& right,
                                    Direction dir) const noexcept;

  std::array<Array1d, 2> operator()(const PerfectGas<Dim>& equation,
                                    const CompleteArray& left,
                                    const CompleteArray& right,
                                    const MaskArray& mask,
                                    Direction dir) const noexcept;
};

extern template struct EinfeldtSignalVelocities<PerfectGas<1>>;
extern template struct EinfeldtSignalVelocities<PerfectGas<2>>;
extern template struct EinfeldtSignalVelocities<PerfectGas<3>>;

// extern template class FluxMethod<Godunov<PerfectGas<1>>>;
// extern template class FluxMethod<MusclHancock<PerfectGas<1>>>;
// extern template class FluxMethod<Godunov<PerfectGas<2>>>;
// extern template class FluxMethod<MusclHancock<PerfectGas<2>>>;
// extern template class FluxMethod<Godunov<PerfectGas<3>>>;
// extern template class FluxMethod<MusclHancock<PerfectGas<3>>>;

} // namespace fub

#endif
