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

#ifndef FUB_SOLVER_EULER_HLLE_RIEMANN_SOLVER_HPP
#define FUB_SOLVER_EULER_HLLE_RIEMANN_SOLVER_HPP

#include "fub/solver/DimensionalSplitFluxMethod.hpp"
#include "fub/solver/euler/IdealGas.hpp"

namespace fub::euler {

/// This class is a first order flux method for the IdealGas problem which uses
/// the approximative HLL riemann problem solver with Einfeldts signals
/// velocivties.
struct HlleRiemannSolver
    : fub::DimensionalSplitFluxMethod<IdealGas::FluxState,
                                      IdealGas::CompleteState> {
  using FluxState = IdealGas::FluxState;
  using CompleteState = IdealGas::CompleteState;

  /// Return a stable time step size for a forward euler time integration. If
  /// you use this method as part of a more complicated scheme you might want to
  /// scale this further with some lower CFL condition.
  double
  computeStableDtOnPatch(const CompleteState& state,
                         const SAMRAI::hier::Patch& patch) const override;

  /// This class method writes the flux quantities into fluxes.
  /// The algorithm uses the Einfeldt signal velocities without and the HLL
  /// method without any further correction term for discontiuity wave.
  void computeFluxesOnPatch(const FluxState& flux, const CompleteState& state,
                            const SAMRAI::hier::Patch& patch, double dt,
                            Direction dir) const override;

  /// Returns one in each direction.
  SAMRAI::hier::IntVector
  getStencilWidth(const SAMRAI::tbox::Dimension& dim) const override;
};

} // namespace fub::euler

#endif