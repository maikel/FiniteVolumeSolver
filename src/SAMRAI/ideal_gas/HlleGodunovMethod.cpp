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

#include "fub/SAMRAI/ideal_gas/HlleGodunovMethod.hpp"

#include "fub/SAMRAI/ideal_gas/IdealGasEquation.hpp"
#include "fub/SAMRAI/utility.hpp"

#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/mechanism/Burke2012.hpp"

namespace fub {
namespace samrai {
namespace ideal_gas {
namespace {
struct HllSignals {
  double left;
  double right;
};

struct HlleSignalState {
  double density;
  double momentum;
  double speed_of_sound;
};

HllSignals
ComputeEinfeldtSignalVelocities_(const HlleSignalState& left,
                                 const HlleSignalState& right) noexcept {
  const double aL2 = left.speed_of_sound * left.speed_of_sound;
  const double aR2 = right.speed_of_sound * right.speed_of_sound;
  const double sqRhoL = std::sqrt(left.density);
  const double sqRhoR = std::sqrt(right.density);
  const double uL = left.momentum / left.density;
  const double uR = right.momentum / right.density;
  const double roeU = (sqRhoL * uL + sqRhoR * uR) / (sqRhoL + sqRhoR);
  const double roeA = std::sqrt(
      (sqRhoL * aL2 + sqRhoR * aR2) / (sqRhoL + sqRhoR) +
      0.5 * (sqRhoL * sqRhoR) / ((sqRhoL + sqRhoR) * (sqRhoL + sqRhoR)) *
          (uR - uL) * (uR - uL));
  return {std::min(uL - aL, roeU - 0.5 * roeA),
          std::max(uR + aR, roeU + 0.5 * roeA)};
}
} // namespace

/// Return a stable time step size for a forward euler time integration. If
/// you use this method as part of a more complicated scheme you might want to
/// scale this further with some lower CFL condition.
double
HlleGodunovMethod::ComputeStableDtOnPatch(const CompletePatchData& state,
                                          const SAMRAI::hier::Patch& patch,
                                          Direction dir) const {
  double velocity = 0.0;
  const SAMRAI::hier::Box& patch_box = patch.getBox();
  const int dir_v = static_cast<int>(dir);
  auto first = SAMRAI::pdat::FaceGeometry::begin(patch_box, dir_v);
  auto last = SAMRAI::pdat::FaceGeometry::end(patch_box, dir_v);
  while (first != last) {
    using FaceIndex = SAMRAI::pdat::FaceIndex;
    FaceIndex face = *first++;
    SAMRAI::pdat::CellIndex left(face.toCell(FaceIndex::Lower));
    SAMRAI::pdat::CellIndex right(face.toCell(FaceIndex::Upper));

    HllSignals signals = ComputeEinfeldtSignalVelocities_(
        {state.density(left), state.momentum(left, dir_v),
         state.speed_of_sound(left)},
        {state.density(right), state.momentum(right, dir_v),
         state.speed_of_sound(right)});

    velocity =
        std::max({velocity, std::abs(signals.left), std::abs(signals.right)});
    FUB_ASSERT(velocity >= 0.0);
  }
  const double* dx = GetCartesianPatchGeometry(patch)->getDx();
  const double max_time_step_size = 0.5 * dx[dir_v] / velocity;
  return max_time_step_size;
}

void HlleGodunovMethod::ComputeFluxesOnPatch(const FluxPatchData& flux,
                                             const CompletePatchData& state,
                                             const SAMRAI::hier::Patch& patch,
                                             double /* dt */,
                                             Direction dir) const {
  const SAMRAI::hier::Box& patch_box = patch.getBox();
  const int dir_v = static_cast<int>(dir);
  auto first = SAMRAI::pdat::FaceGeometry::begin(patch_box, dir_v);
  auto last = SAMRAI::pdat::FaceGeometry::end(patch_box, dir_v);
  while (first != last) {
    using FaceIndex = SAMRAI::pdat::FaceIndex;
    const FaceIndex face = *first++;
    const SAMRAI::pdat::CellIndex left(face.toCell(FaceIndex::Lower));
    const SAMRAI::pdat::CellIndex right(face.toCell(FaceIndex::Upper));

    const double rhoL = state.density(left);
    const double rhoR = state.density(right);
    const double rhouL = state.momentum(left, dir_v);
    const double rhouR = state.momentum(right, dir_v);
    const double aL = state.speed_of_sound(left);
    const double aR = state.speed_of_sound(right);

    HllSignals signals =
        ComputeEinfeldtSignalVelocities_({rhoL, rhouL, aL}, {rhoR, rhouR, aR});

    const double sL = std::min(0.0, signals.left);
    const double sR = std::max(0.0, signals.right);

    const auto hlleFlux = [&](const auto& qL, const auto& qR, const auto& fL,
                              const auto& fR) {
      FUB_ASSERT(sL < sR);
      return (sL * fR - sR * fL + sL * sR * (qL - qR)) / (sR - sL);
    };

    const double rhoeL = state.energy(left);
    const double rhoeR = state.energy(right);
    const double pL = state.pressure(left);
    const double pR = state.pressure(right);

    flux.density(face) = hlleFlux(rhoL, rhoR, rhouL, rhouR);

    const double uL = rhouL / rhoL;
    const double uR = rhouR / rhoR;

    const int dim = patch.getDim().getValue();
    for (int d = 0; d < dim; ++d) {
      if (d == dir_v) {
        flux.momentum(face, d) =
            hlleFlux(rhouL, rhouR, uL * rhouL + pL, uR * rhouR + pR);
      } else {
        flux.momentum(face, d) =
            hlleFlux(rhouL, rhouR, uL * state.momentum(left, d),
                     state.momentum(right, d));
      }
    }
    flux.energy(face) =
        hlleFlux(rhoeL, rhoeR, uL * (rhoeL + pL), uR * (rhoeR + pR));

    const int n_species = state.species.getDepth();
    for (int s = 0; s < n_species; ++s) {
      flux.species(face, s) =
          hlleFlux(state.species(left, s), state.species(right, s),
                   uL * state.species(left, s), uR * state.species(right, s));
    }
  }
}

SAMRAI::hier::IntVector
HlleGodunovMethod::GetStencilWidth(const SAMRAI::tbox::Dimension& dim) const {
  return SAMRAI::hier::IntVector::getOne(dim);
}

} // namespace ideal_gas
} // namespace samrai
} // namespace fub