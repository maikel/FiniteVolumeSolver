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

#ifndef FUB_SAMRAI_HLL_FLUX_METHOD_HPP
#define FUB_SAMRAI_HLL_FLUX_METHOD_HPP

#include "fub/SAMRAI/ideal_gas/HlleRiemannSolver.hpp"

#include "fub/SAMRAI/ideal_gas/IdealGasEquation.hpp"
#include "fub/SAMRAI/utility.hpp"

#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/mechanism/Burke2012.hpp"

namespace fub {
namespace ideal_gas {
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
  span<const T> species;
};

template <typename T> struct HllCons {
  T density;
  T momentum;
  T energy;
  std::vector<T> species;
};

template <typename T> struct HllFlux {
  T density;
  T momentum;
  T energy;
  std::vector<T> species;
};

template <typename T> struct HllResult {
  HllCons<T> state;
  HllFlux<T> flux;
};

template <typename T> struct HlleSignalState {
  T density;
  T momentum;
  T speed_of_sound;
};

double computeHllFlux_(double sL, double sR, double uL, double uR, double fL,
                       double fR) {
  return (sR * fL - sL * fR + sL * sR * (uR - uL)) / (sR - sL);
}

HllResult<double> computeHllFlux_(const HllSignals<double>& signals,
                                  const HllState<double>& left,
                                  const HllState<double>& right) noexcept {
  const double sL = signals.left;
  const double sR = signals.right;
  const double sLsR = sL * sR;
  const double ds = sR - sL;
  auto hllStar_ = [&](double uL, double uR, double fL, double fR) {
    return (fL - fR + sR * uR - sL * uL) / ds;
  };
  auto hllFlux_ = [&](double uL, double uR, double fL, double fR) {
    return (sR * fL - sL * fR + sLsR * (uR - uL)) / ds;
  };
  const double rhoL = left.density;
  const double rhoR = right.density;
  const double rhouL = left.momentum;
  const double rhouR = right.momentum;
  double hllRho = hllStar_(rhoL, rhoR, rhouL, rhouR);
  double hllfRho = hllFlux_(rhoL, rhoR, rhouL, rhouR);

  const double uL = rhouL / rhoL;
  const double uR = rhouR / rhoR;
  const double pL = left.pressure;
  const double pR = right.pressure;
  auto f_rhouL = rhouL * uL + pL;
  auto f_rhouR = rhouR * uR + pR;
  double hllRhou = hllStar_(rhouL, rhouR, f_rhouL, f_rhouR);
  double hllfRhou = hllFlux_(rhouL, rhouR, f_rhouL, f_rhouR);

  const double rhoeL = left.energy;
  const double rhoeR = right.energy;
  auto f_rhoeL = uL * (rhoeL + pL);
  auto f_rhoeR = uR * (rhoeR + pR);
  double hllRhoe = hllStar_(rhoeL, rhoeR, f_rhoeL, f_rhoeR);
  double hllfRhoe = hllFlux_(rhoeL, rhoeR, f_rhoeL, f_rhoeR);

  const int n_species = left.species.size();
  std::vector<double> hllY(n_species);
  std::vector<double> hllfY(n_species);

  span<const double> yL = left.species;
  span<const double> yR = right.species;
  for (int s = 0; s < n_species; ++s) {
    hllY[s] = hllStar_(yL[s], yR[s], uL * yL[s], uR * yR[s]);
    hllfY[s] = hllFlux_(yL[s], yR[s], uL * yL[s], uR * yR[s]);
  }

  return {{hllRho, hllRhou, hllRhoe, std::move(hllY)},
          {hllfRho, hllfRhou, hllfRhoe, std::move(hllfY)}};
}

HllSignals<double> ComputeEinfeldtSignalVelocities_(
    const HlleSignalState<double>& left,
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
      0.5 * (sqRhoL * sqRhoR) / ((sqRhoL + sqRhoR) * (sqRhoL + sqRhoR)) *
          (uR - uL) * (uR - uL));
  const double sL1 = uL - aL;
  const double sL2 = roeU - 0.5 * roeA;
  const double sR1 = roeU + 0.5 * roeA;
  const double sR2 = uR + aR;
  return {std::min(sL1, sL2), std::max(sR1, sR2)};
}
} // namespace

/// Return a stable time step size for a forward euler time integration. If
/// you use this method as part of a more complicated scheme you might want to
/// scale this further with some lower CFL condition.
double
HlleRiemannSolver::ComputeStableDtOnPatch(const CompletePatchData& state,
                                          const SAMRAI::hier::Patch& patch,
                                          Direction dir) const {
  double velocity = 0.0;
  const SAMRAI::hier::Box& patch_box = patch.getBox();
  const int dir_value = static_cast<int>(dir);
  auto first = SAMRAI::pdat::FaceGeometry::begin(patch_box, dir_value);
  auto last = SAMRAI::pdat::FaceGeometry::end(patch_box, dir_value);
  while (first != last) {
    using FaceIndex = SAMRAI::pdat::FaceIndex;
    FaceIndex face = *first++;
    SAMRAI::pdat::CellIndex left(face.toCell(FaceIndex::Lower));
    SAMRAI::pdat::CellIndex right(face.toCell(FaceIndex::Upper));
    const double rhoL = state.density(left);
    const double rhoR = state.density(right);
    const double rhouL = state.momentum(left, dir_value);
    const double rhouR = state.momentum(right, dir_value);
    const double aL = state.speed_of_sound(left);
    const double aR = state.speed_of_sound(right);
    HllSignals<double> signals =
        ComputeEinfeldtSignalVelocities_({rhoL, rhouL, aL}, {rhoR, rhouR, aR});
    const double& bL = signals.left;
    const double& bR = signals.right;
    const double sL = std::min(0.0, bL);
    const double sR = std::max(0.0, bR);

    velocity = std::max({velocity, -sL, sR});
    FUB_ASSERT(velocity >= 0.0);
  }
  const double* dx = GetCartesianPatchGeometry(patch)->getDx();
  const double max_time_step_size = 0.5 * dx[dir_value] / velocity;
  return max_time_step_size;
}

void HlleRiemannSolver::ComputeFluxesOnPatch(const FluxPatchData& flux,
                                             const CompletePatchData& state,
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

    HllSignals<double> signals =
        ComputeEinfeldtSignalVelocities_({rhoL, rhouL, aL}, {rhoR, rhouR, aR});
    const double bL = signals.left;
    const double bR = signals.right;

    const double sL = std::min(0.0, bL);
    const double sR = std::max(0.0, bR);
    const double rhoeL = state.energy(left);
    const double rhoeR = state.energy(right);
    const double pL = state.pressure(left);
    const double pR = state.pressure(right);
    std::vector<double> rhoyL(state.species.getDepth());
    std::vector<double> rhoyR(state.species.getDepth());
    CopyMassFractions(rhoyL, state.species, left);
    CopyMassFractions(rhoyR, state.species, right);

    const HllResult<double> hll =
        computeHllFlux_({sL, sR}, {rhoL, rhouL, rhoeL, pL, rhoyL},
                        {rhoR, rhouR, rhoeR, pR, rhoyR});
    const double fRho = hll.flux.density;
    const double fRhoU = hll.flux.momentum;
    const double fRhoE = hll.flux.energy;

    flux.density(face) = fRho;

    const double uL = rhouL / rhoL;
    const double uR = rhouR / rhoR;

    const int dim = patch.getDim().getValue();
    for (int d_ = 0; d_ < dim; ++d_) {
      const double rhoUL_ = state.momentum(left, d_);
      const double rhoUR_ = state.momentum(right, d_);
      const double fRhoUL = uL * rhoUL_;
      const double fRhoUR = uR * rhoUR_;
      flux.momentum(face, d_) =
          (d_ == d) ? fRhoU
                    : computeHllFlux_(sL, sR, rhoUL_, rhoUR_, fRhoUL, fRhoUR);
    }
    flux.energy(face) = fRhoE;
    const int species_size = state.species.getDepth();
    for (int s = 0; s < species_size; ++s) {
      flux.species(face, s) = hll.flux.species[s];
    }
  }
}

SAMRAI::hier::IntVector
HlleRiemannSolver::GetStencilWidth(const SAMRAI::tbox::Dimension& dim) const {
  return SAMRAI::hier::IntVector::getOne(dim);
}

} // namespace ideal_gas
} // namespace fub

#endif