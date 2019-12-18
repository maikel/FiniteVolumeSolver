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

#ifndef FUB_EQUATIONS_COMPRESSIBLE_ADVECTION_HPP
#define FUB_EQUATIONS_COMPRESSIBLE_ADVECTION_HPP

#include "fub/State.hpp"
#include "fub/StateArray.hpp"

#include "fub/CompleteFromCons.hpp"
#include "fub/Equation.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/ext/Eigen.hpp"

#include "fub/AMReX/IntegratorContext.hpp"

#include <array>

namespace fub {

/// This is a template class for constructing conservative states for the
/// perfect gas equations.
template <typename Density, typename Momentum, typename PTDensity>
struct CompressibleAdvectionConservative {
  Density density;
  Momentum momentum;
  PTDensity PTdensity;
};

// We "register" the conservative state with our framework.
// This enables us to name and iterate over all member variables in a given
// conservative state.
template <typename... Xs>
struct StateTraits<CompressibleAdvectionConservative<Xs...>> {
  static constexpr auto names =
      std::make_tuple("Density", "Momentum", "PTdensity");

  static constexpr auto pointers_to_member =
      std::make_tuple(&CompressibleAdvectionConservative<Xs...>::density,
                      &CompressibleAdvectionConservative<Xs...>::momentum,
                      &CompressibleAdvectionConservative<Xs...>::PTdensity);
};

template <typename Density, typename Momentum, typename PTDensity,
          typename Velocity, typename PTInverse>
struct CompressibleAdvectionComplete
    : CompressibleAdvectionConservative<Density, Momentum, PTDensity> {
  Velocity velocity;
  PTInverse PTinverse;
};

// We "register" the complete state with our framework.
// This enables us to name and iterate over all member variables in a given
// compete state.
template <typename... Xs>
struct StateTraits<CompressibleAdvectionComplete<Xs...>> {
  static constexpr auto names = std::make_tuple(
      "Density", "Momentum", "PTdensity", "Velocity", "PTinverse");
  static constexpr auto pointers_to_member =
      std::make_tuple(&CompressibleAdvectionComplete<Xs...>::density,
                      &CompressibleAdvectionComplete<Xs...>::momentum,
                      &CompressibleAdvectionComplete<Xs...>::PTdensity,
                      &CompressibleAdvectionComplete<Xs...>::velocity,
                      &CompressibleAdvectionComplete<Xs...>::PTinverse);
};

template <int Rank>
using CompressibleAdvectionConsShape =
    CompressibleAdvectionConservative<ScalarDepth, VectorDepth<Rank>,
                                      ScalarDepth>;

template <int Rank>
using CompressibleAdvectionCompleteShape =
    CompressibleAdvectionComplete<ScalarDepth, VectorDepth<Rank>, ScalarDepth,
                                  VectorDepth<Rank>, ScalarDepth>;

template <int N> struct CompressibleAdvection {
  using ConservativeDepths = CompressibleAdvectionConsShape<N>;
  using CompleteDepths = CompressibleAdvectionCompleteShape<N>;

  using Conservative = ::fub::Conservative<CompressibleAdvection<N>>;
  using Complete = ::fub::Complete<CompressibleAdvection<N>>;

  static constexpr int Rank() noexcept { return N; }

  void CompleteFromCons(Complete& state, const Conservative& cons) {
    state.density = cons.density;
    state.momentum = cons.momentum;
    state.PTdensity = cons.PTdensity;
    state.velocity = cons.momentum / cons.density;
    state.PTinverse = cons.density / cons.PTdensity;
  }

  /// Universal gas constant
  double R{8.314}; // [J / (mol K)]

  /// Heat capacity ratio at constant pressure
  double c_p{2.5 * R};

  /// Gravitational acceleration
  double g{10.0};  //  [m / s^2]

  /// Coriolis parameter in beta plane
  double f{0.0};
};

template <int SpaceDimension> struct CompressibleAdvectionFluxMethod {
  using Conservative =
      typename CompressibleAdvection<SpaceDimension>::Conservative;
  using Complete = typename CompressibleAdvection<SpaceDimension>::Complete;

  constexpr static int GetStencilWidth() { return 2; }

  CompressibleAdvection<2> GetEquation() const noexcept { return {}; }

  Duration ComputeStableDt(amrex::IntegratorContext& context, int level,
                           Direction dir);

  static Duration
  ComputeStableDt(const View<const Complete>& states,
                  const StridedDataView<const double, SpaceDimension> Pv,
                  double dx, Direction dir);

  static Conservative
  ComputeNumericFluxes(const std::array<Complete, 4>& stencil,
                       const std::array<double, 5> Pvs, Duration dt, double dx,
                       Direction dir);

  /// Pv is face centered
  static void
  ComputeNumericFluxes(const View<Conservative>& fluxes,
                       const View<const Complete>& states,
                       const StridedDataView<const double, SpaceDimension>& Pv,
                       Duration dt, double dx, Direction dir);

  static void ComputeNumericFluxes(const View<Conservative>& fluxes,
                                   const View<const Complete>& states,
                                   Duration dt, double dx, Direction dir);

  static void ComputeNumericFluxes(execution::SequentialTag,
                                   const View<Conservative>& fluxes,
                                   const View<const Complete>& states,
                                   Duration dt, double dx, Direction dir) {
    ComputeNumericFluxes(fluxes, states, dt, dx, dir);
  }

  std::function<double(std::array<double, SpaceDimension>, Duration, Direction)>
      Pv_function_;
};

// We define this class only for dimensions 1 to 3.
// The definitions will be found in its source file PerfetGas.cpp
// extern template struct CompressibleAdvection<2>;
extern template struct CompressibleAdvectionFluxMethod<2>;

} // namespace fub

#endif
