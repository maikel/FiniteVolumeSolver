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

#ifndef FUB_SOLVER_HLL_FLUX_METHOD_HPP
#define FUB_SOLVER_HLL_FLUX_METHOD_HPP

#include "fub/solver/euler/HlleRiemannSolver.hpp"
#include "fub/SAMRAI/utility.hpp"

namespace fub {
namespace euler {

namespace {
template <typename T> struct HllSignals {
  T left;
  T right;
};

template <typename T> struct HllState {
  T density;
  T momentum;
  T energy;
  T pressure;
};

template <typename T> struct HllFlux {
  T density;
  T momentum;
  T energy;
};

template <typename T> struct HlleSignalState {
  T density;
  T momentum;
  T speed_of_sound;
};

double computeHllFlux_(double sL, double sR, double uL, double uR, double fL,
                       double fR) {
  if (0 <= sL) {
    return fL;
  }
  if (0 >= sR) {
    return fR;
  }
  FUB_ASSERT(sL < sR);
  return (sR * fL - sL * fR + sL * sR * (uR - uL)) / (sR - sL);
}

HllFlux<double> computeHllFlux_(const HllSignals<double>& signals,
                                const HllState<double>& left,
                                const HllState<double>& right) noexcept {
  const double sL = signals.left;
  const double sR = signals.right;
  const double sLsR = sL * sR;
  const double ds = sR - sL;
  int maskL = (0.0 < sL);
  int maskR = (0.0 > sR);
  auto hllFlux_ = [&](double uL, double uR, double fL, double fR) {
    double hll = (sR * fL - sL * fR + sLsR * (uR - uL)) / ds;
    hll = maskL * fL + (1 - maskL) * hll;
    hll = maskR * fR + (1 - maskR) * hll;
    return hll;
  };
  const double rhoL = left.density;
  const double rhoR = right.density;
  const double rhouL = left.momentum;
  const double rhouR = right.momentum;
  double hllRho = hllFlux_(rhoL, rhoR, rhouL, rhouR);

  const double pL = left.pressure;
  const double pR = right.pressure;
  auto f_rhouL = rhouL * rhouL / rhoL + pL;
  auto f_rhouR = rhouR * rhouR / rhoR + pR;
  double hllRhou = hllFlux_(rhouL, rhouR, f_rhouL, f_rhouR);

  const double rhoeL = left.pressure;
  const double rhoeR = right.pressure;
  auto f_rhoeL = rhouL / rhoL * (rhoeL + pL);
  auto f_rhoeR = rhouR / rhoR * (rhoeR + pR);
  double hllRhoe = hllFlux_(rhoeL, rhoeR, f_rhoeL, f_rhoeR);

  return {hllRho, hllRhou, hllRhoe};
}

HllSignals<double>
computeHlleSignalVelocities_(const HlleSignalState<double>& left,
                             const HlleSignalState<double>& right) noexcept {
  const double rhoL = left.density;
  const double rhoR = right.density;
  const double rhoUL = left.momentum;
  const double rhoUR = right.momentum;
  const double aL = left.speed_of_sound;
  const double aR = right.speed_of_sound;
  const double sqRhoL = std::sqrt(rhoL);
  const double sqRhoR = std::sqrt(rhoR);
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  const double roeU = (sqRhoL * uL + sqRhoR * uR) / (sqRhoL + sqRhoR);
  const double roeA = std::sqrt(
      (sqRhoL * aL * aL + sqRhoR * aR * aR) / (sqRhoL + sqRhoR) +
      0.5f * (sqRhoL * sqRhoR) / ((sqRhoL + sqRhoR) * (sqRhoL + sqRhoR)) *
          (uR - uL) * (uR - uL));
  const double sL1 = uL - 0.5 * aL;
  const double sL2 = roeU - roeA;
  const double sR1 = roeU + roeA;
  const double sR2 = uR + 0.5 * aR;
  return {std::min(sL1, sL2), std::max(sR1, sR2)};
}
} // namespace

/// Return a stable time step size for a forward euler time integration. If
/// you use this method as part of a more complicated scheme you might want to
/// scale this further with some lower CFL condition.
double HlleRiemannSolver::computeStableDtOnPatch(
    const CompleteState& state, const SAMRAI::hier::Patch& patch) const {

  double velocity = 0.0;
  const SAMRAI::hier::Box& patch_box = patch.getBox();
  auto first = SAMRAI::pdat::FaceGeometry::begin(patch_box, 0);
  auto last = SAMRAI::pdat::FaceGeometry::end(patch_box, 0);
  while (first != last) {
    using FaceIndex = SAMRAI::pdat::FaceIndex;
    FaceIndex face = *first++;
    SAMRAI::pdat::CellIndex left(face.toCell(FaceIndex::Lower));
    SAMRAI::pdat::CellIndex right(face.toCell(FaceIndex::Upper));
    const double rhoL = state.density(left);
    const double rhoR = state.density(right);
    const double rhouL = state.momentum(left, 0);
    const double rhouR = state.momentum(right, 0);
    const double aL = state.speed_of_sound(left);
    const double aR = state.speed_of_sound(right);

#ifdef __cpp_structured_bindings
    auto [sL, sR] =
        computeHlleSignalVelocities_({rhoL, rhouL, aL}, {rhoR, rhouR, aR});
#else
    HllSignals<double> signals =
        computeHlleSignalVelocities_({rhoL, rhouL, aL}, {rhoR, rhouR, aR});
    const double& sL = signals.left;
    const double& sR = signals.right;
#endif
    velocity = std::max({velocity, std::abs(sL), std::abs(sR)});
  }
  const double* dx = getCartesianPatchGeometry(patch)->getDx();
  return 0.5 * dx[0] / velocity;
}

/// This class method writes the flux quantities into fluxes.
/// The algorithm uses the Einfeldt signal velocities without and the HLL
/// method without any further correction term for discontinuity wave.
void HlleRiemannSolver::computeFluxesOnPatch(const FluxState& flux,
                                             const CompleteState& state,
                                             const SAMRAI::hier::Patch& patch,
                                             double dt, Direction dir) const {
  const SAMRAI::hier::Box& patch_box = patch.getBox();
  const int d = static_cast<int>(dir);
  auto first = SAMRAI::pdat::FaceGeometry::begin(patch_box, d);
  auto last = SAMRAI::pdat::FaceGeometry::end(patch_box, d);
  while (first != last) {
    using FaceIndex = SAMRAI::pdat::FaceIndex;
    const FaceIndex face = *first++;
    const SAMRAI::pdat::CellIndex left(face.toCell(FaceIndex::Lower));
    const SAMRAI::pdat::CellIndex right(face.toCell(FaceIndex::Upper));

    const double rhoL = state.density(left);
    const double rhoR = state.density(right);
    const double rhouL = state.momentum(left, d);
    const double rhouR = state.momentum(right, d);
    const double aL = state.speed_of_sound(left);
    const double aR = state.speed_of_sound(right);

#ifdef __cpp_structured_bindings
    auto [sL, sR] =
        computeHlleSignalVelocities_({rhoL, rhouL, aL}, {rhoR, rhouR, aR});
#else
    HllSignals<double> signals =
        computeHlleSignalVelocities_({rhoL, rhouL, aL}, {rhoR, rhouR, aR});
    const double& sL = signals.left;
    const double& sR = signals.right;
#endif

    const double rhoeL = state.energy(left);
    const double rhoeR = state.energy(right);
    const double pL = state.pressure(left);
    const double pR = state.pressure(right);

#ifdef __cpp_structured_bindings
    const auto [fRho, fRhoU, fRhoE] = computeHllFlux_(
        {sL, sR}, {rhoL, rhouL, rhoeL, pL}, {rhoR, rhouR, rhoeR, pR});
#else 
    const HllFlux<double> flux_ = computeHllFlux_(
        {sL, sR}, {rhoL, rhouL, rhoeL, pL}, {rhoR, rhouR, rhoeR, pR});
    const double& fRho = flux_.density;
    const double& fRhoU = flux_.momentum;
    const double& fRhoE = flux_.energy;
#endif

    flux.density(face) = fRho;

    const double uL = rhouL / rhoL;
    const double uR = rhouR / rhoR;

    const int dim = patch.getDim().getValue();
    for (int d_ = 0; d_ < dim; ++d_) {
      const double fRhoUL = uL * state.momentum(left, d_);
      const double fRhoUR = uR * state.momentum(right, d_);
      flux.momentum(face, d_) =
          (d_ == d)
              ? fRhoU
              : computeHllFlux_(sL, sR, state.momentum(left, d_),
                                state.momentum(right, d_), fRhoUL, fRhoUR);
    }
    flux.energy(face) = fRhoE;
    const int species_size = state.species.getDepth();
    for (int s = 0; s < species_size; ++s) {
      const double speciesL = state.species(left, s);
      const double speciesR = state.species(right, s);
      const double fSpeciesL = uL * speciesL;
      const double fSpeciesR = uR * speciesR;
      flux.species(face, s) =
          computeHllFlux_(sL, sR, speciesL, speciesR, fSpeciesL, fSpeciesR);
    }
  }
}

SAMRAI::hier::IntVector
HlleRiemannSolver::getStencilWidth(const SAMRAI::tbox::Dimension& dim) const {
  return SAMRAI::hier::IntVector::getOne(dim);
}

} // namespace euler
} // namespace fub

#endif