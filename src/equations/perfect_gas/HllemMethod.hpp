// Copyright (c) 2020 Maikel Nadolski
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

#include "fub/equations/perfect_gas/HllemMethod.hpp"

namespace fub::perfect_gas {

template <typename EulerEquation, bool Larrouturou>
void Hllem<EulerEquation, Larrouturou>::SolveRiemannProblem(
    Complete& solution, const Complete& left, const Complete& right,
    Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  int d = static_cast<int>(dir);
  left_ = left;
  right_ = right;
  const double Msq = euler::MaSq(equation_);
  const double Minv = 1.0 / std::sqrt(Msq);
  const double rhoL = fub::euler::Density(equation_, left_);
  const double rhoR = fub::euler::Density(equation_, right_);
  [[maybe_unused]] Array<double, 1, Dim> vL_ =
      fub::euler::Velocity(equation_, left_);
  [[maybe_unused]] Array<double, 1, Dim> vR_ =
      fub::euler::Velocity(equation_, right_);
  [[maybe_unused]] double eKinRL = 0.0;
  [[maybe_unused]] double eKinRR = 0.0;
  // If we do the Larrouturou version with velocity we need to set the radial
  // velocity to zero.
  if constexpr (Larrouturou) {
    std::array<int, Dim> dim_mask{};
    for (int i = 0; i < Dim; ++i) {
      dim_mask[i] = i == d;
      eKinRL += (1 - dim_mask[i]) * vL_[i] * vL_[i];
      eKinRR += (1 - dim_mask[i]) * vR_[i] * vR_[i];
      left_.momentum[i] *= dim_mask[i];
      right_.momentum[i] *= dim_mask[i];
    }
    eKinRL *= 0.5 * Msq;
    eKinRR *= 0.5 * Msq;
    left_.energy -= rhoL * eKinRL;
    right_.energy -= rhoR * eKinRR;
  }
  fub::Flux(equation_, fluxL_, left_, dir);
  fub::Flux(equation_, fluxR_, right_, dir);

  // Compute Einfeldt signals velocities

  const double rhoUL = fub::euler::Momentum(equation_, left_, d);
  const double rhoUR = fub::euler::Momentum(equation_, right_, d);
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  const double aL = fub::euler::SpeedOfSound(equation_, left_);
  const double aR = fub::euler::SpeedOfSound(equation_, right_);
  const double rhoEL = fub::euler::Energy(equation_, left_);
  const double rhoER = fub::euler::Energy(equation_, right_);
  const double pL = fub::euler::Pressure(equation_, left_);
  const double pR = fub::euler::Pressure(equation_, right_);
  const double hL = (rhoEL + pL) / rhoL;
  const double hR = (rhoER + pR) / rhoR;
  const double sqRhoL = std::sqrt(rhoL);
  const double sqRhoR = std::sqrt(rhoR);
  const double sqRhoSum = sqRhoL + sqRhoR;
  Array<double, 1, Dim> vL = fub::euler::Velocity(equation_, left_);
  Array<double, 1, Dim> vR = fub::euler::Velocity(equation_, right_);
  FUB_ASSERT(sqRhoSum > 0.0);
  Array<double, 1, Dim> roeU = (sqRhoL * vL + sqRhoR * vR) / sqRhoSum;
  const double roeU0 = roeU[d];
  const double roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;

  const double gammaL = fub::euler::Gamma(equation_, left_);
  const double gammaR = fub::euler::Gamma(equation_, right_);
  const double roeGamma = (sqRhoL * gammaL + sqRhoR * gammaR) / sqRhoSum;
  const double gm1 = roeGamma - 1.0;
  // const double beta = gm1 / (2 * roeGamma);

  const double roeA2 = gm1 * (roeH - 0.5 * Msq * roeU.matrix().squaredNorm());
  const double roeA = std::sqrt(roeA2);
  const double maxA = std::max(aL, aR);
  const double sL1 = uL - maxA * Minv;
  const double sL2 = roeU0 - roeA * Minv;
  const double sR1 = roeU0 + roeA * Minv;
  const double sR2 = uR + maxA * Minv;

  const double sL = std::min(sL1, sL2);
  const double sR = std::max(sR1, sR2);

  const double bL = std::min(sL, 0.0);
  const double bR = std::max(sR, 0.0);
  const double db = bR - bL;
  const double db_positive = int(db <= 0) + int(db > 0) * db;

  ForEachComponent(
      [&](double& w, double fL, double fR, double qL, double qR) {
        w = (bR * qR - bL * qL + fL - fR) / db_positive;
      },
      w_hlle_, fluxL_, fluxR_, AsCons(left_), AsCons(right_));

  const double squaredNormRoeU = Msq * roeU.matrix().squaredNorm();
  const double squaredNormRoeU_half = 0.5 * squaredNormRoeU;
  const double u_bar = 0.5 * (sR + sL);
  const double u_bar_abs = std::abs(u_bar);
  const double delta = roeA / (roeA + u_bar_abs);

  const double deltaRho = rhoR - rhoL;
  const Array<double, 1, Dim> deltaRhoU = right_.momentum - left_.momentum;
  const double deltaRhoE = rhoER - rhoEL;

  if constexpr (Dim == 1) {
    const double l21 = gm1 / roeA2 * (roeH - Msq * roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * Msq * roeU[0];
    const double l23 = gm1 / roeA2 * (-1.0);
    const double alpha_2 =
        l21 * deltaRho + l22 * deltaRhoU[0] + l23 * deltaRhoE;
    const double u_bar_delta_alpha_2 = u_bar * delta * alpha_2;
    w_hllem_.density = w_hlle_.density - u_bar_delta_alpha_2 * 1.0;
    w_hllem_.momentum = w_hlle_.momentum - u_bar_delta_alpha_2 * roeU[0];
    w_hllem_.energy =
        w_hlle_.energy - u_bar_delta_alpha_2 * squaredNormRoeU_half;
    if constexpr (fub::euler::state_with_species<
                                                 Conservative>()) {
      const int n_species = w_hllem_.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left_, i);
        const double rhoYR = fub::euler::Species(equation_, right_, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const int upwind = (u_bar > 0.0);
          w_hllem_.species[i] =
              w_hllem_.density * (upwind * YL + (1 - upwind) * YR);
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          w_hllem_.species[i] = w_hlle_.species[i] -
                                u_bar_delta_alpha_2 * roeY -
                                u_bar * delta * alpha_i;
        }
      }
    }
    if constexpr (fub::euler::state_with_passive_scalars<Conservative>()) {
      const int n_passive_scalars = w_hllem_.passive_scalars.size();
      for (int i = 0; i < n_passive_scalars; ++i) {
        const double rhoYL = left_.passive_scalars[i];
        const double rhoYR = right_.passive_scalars[i];
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const int upwind = (u_bar > 0.0);
          w_hllem_.passive_scalars[i] =
              w_hllem_.density * (upwind * YL + (1 - upwind) * YR);
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          w_hllem_.passive_scalars[i] = w_hlle_.passive_scalars[i] -
                                        u_bar_delta_alpha_2 * roeY -
                                        u_bar * delta * alpha_i;
        }
      }
    }
  } else if constexpr (Dim == 2) {
    const int ix = int(dir);
    const int iy = (ix + 1) % 2;
    const double l21 = gm1 / roeA2 * (roeH - Msq * roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * Msq * roeU[0];
    const double l23 = gm1 / roeA2 * roeU[1];
    const double l24 = gm1 / roeA2 * (-1.0);
    const double alpha_2 = l21 * deltaRho + l22 * deltaRhoU[0] +
                           l23 * deltaRhoU[1] + l24 * deltaRhoE;
    const double l31 = -roeU[iy];
    // const double l32 = 0;
    const double l33 = 1;
    // const double l34 = 0;
    const double alpha_3 = l31 * deltaRho + l33 * deltaRhoU[iy];
    const double u_bar_delta_alpha_2 = u_bar * delta * alpha_2;
    const double u_bar_delta_alpha_3 = u_bar * delta * alpha_3;
    w_hllem_.density = w_hlle_.density - u_bar_delta_alpha_2 * 1.0;
    w_hllem_.momentum[ix] =
        w_hlle_.momentum[ix] - u_bar_delta_alpha_2 * roeU[ix];
    if constexpr (Larrouturou) {
      const int upwind = (u_bar > 0.0);
      w_hllem_.momentum[iy] =
          w_hllem_.density * (upwind * vL_[iy] + (1 - upwind) * vR_[iy]);
      w_hllem_.energy =
          w_hlle_.energy - u_bar_delta_alpha_2 * squaredNormRoeU_half +
          w_hllem_.density * (upwind * eKinRL + (1 - upwind) * eKinRR);
    } else {
      w_hllem_.momentum[iy] = w_hlle_.momentum[iy] -
                              u_bar_delta_alpha_2 * roeU[iy] -
                              u_bar_delta_alpha_3 * 1.0;
      w_hllem_.energy = w_hlle_.energy -
                        u_bar_delta_alpha_2 * squaredNormRoeU_half -
                        u_bar_delta_alpha_3 * roeU[iy];
    }
    if constexpr (fub::euler::state_with_species<
                                                 Conservative>()) {
      const int n_species = w_hllem_.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left_, i);
        const double rhoYR = fub::euler::Species(equation_, right_, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const int upwind = (u_bar > 0.0);
          w_hllem_.species[i] =
              w_hllem_.density * (upwind * YL + (1 - upwind) * YR);
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          w_hllem_.species[i] = w_hlle_.species[i] -
                                u_bar_delta_alpha_2 * roeY -
                                u_bar * delta * alpha_i;
        }
      }
    }
    if constexpr (fub::euler::state_with_passive_scalars<Conservative>()) {
      const int n_passive_scalars = w_hllem_.passive_scalars.size();
      for (int i = 0; i < n_passive_scalars; ++i) {
        const double rhoYL = left_.passive_scalars[i];
        const double rhoYR = right_.passive_scalars[i];
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const int upwind = (u_bar > 0.0);
          w_hllem_.passive_scalars[i] =
              w_hllem_.density * (upwind * YL + (1 - upwind) * YR);
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          w_hllem_.passive_scalars[i] = w_hlle_.passive_scalars[i] -
                                        u_bar_delta_alpha_2 * roeY -
                                        u_bar * delta * alpha_i;
        }
      }
    }
  } else {
    static_assert(Dim == 3);
    const int ix = int(dir);
    const int iy = (ix + 1) % 3;
    const int iz = (iy + 1) % 3;
    const double l21 = gm1 / roeA2 * (roeH - Msq * roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * Msq * roeU[0];
    const double l23 = gm1 / roeA2 * roeU[1];
    const double l24 = gm1 / roeA2 * roeU[2];
    const double l25 = gm1 / roeA2 * (-1.0);
    const double l31 = -roeU[iy];
    // const double l32 = 0;
    const double l33 = 1;
    // const double l34 = 0;
    // const double l34 = 0;
    const double l41 = -roeU[iz];
    // const double l32 = 0;
    // const double l33 = 0;
    const double l44 = 1;
    // const double l34 = 0;
    const double alpha_2 = l21 * deltaRho + l22 * deltaRhoU[0] +
                           l23 * deltaRhoU[1] + l24 * deltaRhoU[2] +
                           l25 * deltaRhoE;
    const double alpha_3 = l31 * deltaRho + l33 * deltaRhoU[iy];
    const double alpha_4 = l41 * deltaRho + l44 * deltaRhoU[iz];
    const double u_bar_delta = u_bar * delta;
    w_hllem_.density = w_hlle_.density - u_bar_delta * alpha_2;
    w_hllem_.momentum[ix] =
        w_hlle_.momentum[ix] - u_bar_delta * alpha_2 * roeU[ix];
    if constexpr (Larrouturou) {
      const int upwind = (u_bar > 0.0);
      w_hllem_.momentum[iy] =
          w_hllem_.density * (upwind * vL_[iy] + (1 - upwind) * vR_[iy]);
      w_hllem_.momentum[iz] =
          w_hllem_.density * (upwind * vL_[iz] + (1 - upwind) * vR_[iz]);
      w_hllem_.energy =
          w_hlle_.energy - u_bar_delta * alpha_2 * squaredNormRoeU_half +
          w_hllem_.density * (upwind * eKinRL + (1 - upwind) * eKinRR);
    } else {
      w_hllem_.momentum[iy] =
          w_hlle_.momentum[iy] -
          u_bar_delta * (alpha_2 * roeU[iy] + alpha_3 * 1.0);
      w_hllem_.momentum[iz] =
          w_hlle_.momentum[iz] -
          u_bar_delta * (alpha_2 * roeU[iz] + alpha_4 * 1.0);
      w_hllem_.energy = w_hlle_.energy -
                        u_bar_delta * (alpha_2 * squaredNormRoeU_half +
                                       alpha_3 * roeU[iy] + alpha_4 * roeU[iz]);
    }
    if constexpr (fub::euler::state_with_species<
                                                 Conservative>()) {
      const int n_species = w_hllem_.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left_, i);
        const double rhoYR = fub::euler::Species(equation_, right_, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const int upwind = (u_bar > 0.0);
          w_hllem_.species[i] =
              w_hllem_.density * (upwind * YL + (1 - upwind) * YR);
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          w_hllem_.species[i] = w_hlle_.species[i] -
                                u_bar_delta * alpha_2 * roeY -
                                u_bar * delta * alpha_i;
        }
      }
    }
    if constexpr (fub::euler::state_with_passive_scalars<Conservative>()) {
      const int n_passive_scalars = w_hllem_.passive_scalars.size();
      for (int i = 0; i < n_passive_scalars; ++i) {
        const double rhoYL = left_.passive_scalars[i];
        const double rhoYR = right_.passive_scalars[i];
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const int upwind = (u_bar > 0.0);
          w_hllem_.passive_scalars[i] =
              w_hllem_.density * (upwind * YL + (1 - upwind) * YR);
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          w_hllem_.passive_scalars[i] = w_hlle_.passive_scalars[i] -
                                u_bar_delta * alpha_2 * roeY -
                                u_bar * delta * alpha_i;
        }
      }
    }
  }

  if (0.0 < sL) {
    solution = left_;
  } else if (sR < 0.0) {
    solution = right_;
  } else {
    AsCons(solution) = w_hllem_;
    CompleteFromCons(equation_, solution, solution);
  }
}

template <typename EulerEquation, bool Larrouturou>
double Hllem<EulerEquation, Larrouturou>::ComputeNumericFlux(
    Conservative& flux, span<const Complete, 2> states, Duration /* dt */,
    double /* dx */, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  int d = static_cast<int>(dir);
  left_ = states[0];
  right_ = states[1];
  const double Msq = euler::MaSq(equation_);
  const double Minv = 1.0 / std::sqrt(Msq);
  const double rhoL = fub::euler::Density(equation_, left_);
  const double rhoR = fub::euler::Density(equation_, right_);
  [[maybe_unused]] Array<double, 1, Dim> vL_ =
      fub::euler::Velocity(equation_, left_);
  [[maybe_unused]] Array<double, 1, Dim> vR_ =
      fub::euler::Velocity(equation_, right_);
  [[maybe_unused]] double eKinRL = 0.0;
  [[maybe_unused]] double eKinRR = 0.0;
  // If we do the Larrouturou version with velocity we need to set the radial
  // velocity to zero.
  if constexpr (Larrouturou) {
    std::array<int, Dim> dim_mask{};
    for (int i = 0; i < Dim; ++i) {
      dim_mask[i] = i == d;
      eKinRL += (1 - dim_mask[i]) * vL_[i] * vL_[i];
      eKinRR += (1 - dim_mask[i]) * vR_[i] * vR_[i];
      left_.momentum[i] *= dim_mask[i];
      right_.momentum[i] *= dim_mask[i];
    }
    eKinRL *= 0.5 * Msq;
    eKinRR *= 0.5 * Msq;
    left_.energy -= rhoL * eKinRL;
    right_.energy -= rhoR * eKinRR;
  }

  fub::Flux(equation_, fluxL_, left_, dir);
  fub::Flux(equation_, fluxR_, right_, dir);

  // Compute Einfeldt signals velocities

  const double rhoUL = fub::euler::Momentum(equation_, left_, d);
  const double rhoUR = fub::euler::Momentum(equation_, right_, d);
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  const double aL = fub::euler::SpeedOfSound(equation_, left_);
  const double aR = fub::euler::SpeedOfSound(equation_, right_);
  const double rhoEL = fub::euler::Energy(equation_, left_);
  const double rhoER = fub::euler::Energy(equation_, right_);
  const double pL = fub::euler::Pressure(equation_, left_);
  const double pR = fub::euler::Pressure(equation_, right_);
  const double hL = (rhoEL + pL) / rhoL;
  const double hR = (rhoER + pR) / rhoR;
  const double sqRhoL = std::sqrt(rhoL);
  const double sqRhoR = std::sqrt(rhoR);
  const double sqRhoSum = sqRhoL + sqRhoR;
  Array<double, 1, Dim> vL = fub::euler::Velocity(equation_, left_);
  Array<double, 1, Dim> vR = fub::euler::Velocity(equation_, right_);
  FUB_ASSERT(sqRhoSum > 0.0);
  Array<double, 1, Dim> roeU = (sqRhoL * vL + sqRhoR * vR) / sqRhoSum;
  const double roeU0 = roeU[d];
  const double roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;

  const double gammaL = fub::euler::Gamma(equation_, left_);
  const double gammaR = fub::euler::Gamma(equation_, right_);
  const double roeGamma = (sqRhoL * gammaL + sqRhoR * gammaR) / sqRhoSum;
  const double gm1 = roeGamma - 1.0;

  const double roeA2 = gm1 * (roeH - 0.5 * Msq * roeU.matrix().squaredNorm());
  const double roeA = std::sqrt(roeA2);
  const double maxA = std::max(aL, aR);
  const double sL1 = uL - maxA * Minv;
  const double sL2 = roeU0 - roeA * Minv;
  const double sR1 = roeU0 + roeA * Minv;
  const double sR2 = uR + maxA * Minv;

  const double sL = std::min(sL1, sL2);
  const double sR = std::max(sR1, sR2);

  const double bL = std::min(sL, 0.0);
  const double bR = std::max(sR, 0.0);
  const double bLbR = bL * bR;
  const double db = bR - bL;
  const double db_positive = int(db <= 0) + int(db > 0) * db;

  ForEachComponent(
      [&](double& nf, double fL, double fR, double qL, double qR) {
        nf = (bR * fL - bL * fR + bLbR * (qR - qL)) / db_positive;
      },
      flux_hlle_, fluxL_, fluxR_, AsCons(left_), AsCons(right_));

  const double squaredNormRoeU = Msq * roeU.matrix().squaredNorm();
  const double squaredNormRoeU_half = 0.5 * squaredNormRoeU;
  const double u_bar = 0.5 * (sR + sL);
  const double u_bar_abs = std::abs(u_bar);
  const double delta = roeA / (roeA + u_bar_abs);
  const double b = bLbR / db_positive;

  const double deltaRho = rhoR - rhoL;
  const Array<double, 1, Dim> deltaRhoU = right_.momentum - left_.momentum;
  const double deltaRhoE = rhoER - rhoEL;

  if constexpr (Dim == 1) {
    const double l21 = gm1 / roeA2 * (roeH - Msq * roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * Msq * roeU[0];
    const double l23 = gm1 / roeA2 * (-1.0);
    const double alpha_2 =
        l21 * deltaRho + l22 * deltaRhoU[0] + l23 * deltaRhoE;
    const double b_delta_alpha_2 = b * delta * alpha_2;
    flux.density = flux_hlle_.density - b_delta_alpha_2 * 1.0;
    flux.momentum = flux_hlle_.momentum - b_delta_alpha_2 * roeU[0];
    flux.energy = flux_hlle_.energy - b_delta_alpha_2 * squaredNormRoeU_half;
    if constexpr (fub::euler::state_with_species<
                                                 Conservative>()) {
      const int n_species = flux.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left_, i);
        const double rhoYR = fub::euler::Species(equation_, right_, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const double f_rho = flux.density;
          const int upwind = (f_rho > 0.0);
          const double f_hllem = f_rho * (upwind * YL - (1 - upwind) * YR);
          flux.species[i] = f_hllem;
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          const double b_delta_alpha_i = b * delta * alpha_i;
          flux.species[i] =
              flux_hlle_.species[i] - b_delta_alpha_2 * roeY - b_delta_alpha_i;
        }
      }
    }
    if constexpr (fub::euler::state_with_passive_scalars<Conservative>()) {
      const int n_passive_scalars = flux.passive_scalars.size();
      for (int i = 0; i < n_passive_scalars; ++i) {
        const double rhoYL = left_.passive_scalars[i];
        const double rhoYR = right_.passive_scalars[i];
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const double f_rho = flux.density;
          const int upwind = (f_rho > 0.0);
          const double f_hllem = f_rho * (upwind * YL - (1 - upwind) * YR);
          flux.passive_scalars[i] = f_hllem;
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          const double b_delta_alpha_i = b * delta * alpha_i;
          flux.passive_scalars[i] = flux_hlle_.passive_scalars[i] -
                                    b_delta_alpha_2 * roeY - b_delta_alpha_i;
        }
      }
    }
  } else if constexpr (Dim == 2) {
    const int ix = int(dir);
    const int iy = (ix + 1) % 2;
    const double l21 = gm1 / roeA2 * (roeH - Msq * roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * Msq * roeU[0];
    const double l23 = gm1 / roeA2 * roeU[1];
    const double l24 = gm1 / roeA2 * (-1.0);
    const double alpha_2 = l21 * deltaRho + l22 * deltaRhoU[0] +
                           l23 * deltaRhoU[1] + l24 * deltaRhoE;
    const double b_delta = b * delta;
    const double b_delta_alpha_2 = b_delta * alpha_2;
    flux.density = flux_hlle_.density - b_delta_alpha_2;
    flux.momentum[ix] = flux_hlle_.momentum[ix] - b_delta_alpha_2 * roeU[ix];
    if constexpr (Larrouturou) {
      const double f_rho = flux.density;
      const int upwind = (f_rho > 0.0);
      const double f_hllem_rhov =
          f_rho * (upwind * vL_[iy] - (1 - upwind) * vR_[iy]);
      const double f_hllem_rhoE =
          f_rho * (upwind * eKinRL - (1 - upwind) * eKinRR);
      flux.momentum[iy] = f_hllem_rhov;
      flux.energy = flux_hlle_.energy - b_delta_alpha_2 * squaredNormRoeU_half +
                    f_hllem_rhoE;
    } else {
      const double l31 = -roeU[iy];
      const double l33 = 1;
      const double alpha_3 = l31 * deltaRho + l33 * deltaRhoU[iy];
      const double b_delta_alpha_3 = b_delta * alpha_3;
      flux.momentum[iy] = flux_hlle_.momentum[iy] - b_delta_alpha_2 * roeU[iy] -
                          b_delta_alpha_3 * 1.0;
      flux.energy = flux_hlle_.energy - b_delta_alpha_2 * squaredNormRoeU_half -
                    b_delta_alpha_3 * roeU[iy];
    }
    if constexpr (fub::euler::state_with_species<
                                                 Conservative>()) {
      const int n_species = flux.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left_, i);
        const double rhoYR = fub::euler::Species(equation_, right_, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const double f_rho = flux.density;
          const int upwind = (f_rho > 0.0);
          const double f_hllem = f_rho * (upwind * YL - (1 - upwind) * YR);
          flux.species[i] = f_hllem;
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          const double b_delta_alpha_i = b * delta * alpha_i;
          flux.species[i] =
              flux_hlle_.species[i] - b_delta_alpha_2 * roeY - b_delta_alpha_i;
        }
      }
    }
    if constexpr (fub::euler::state_with_passive_scalars<Conservative>()) {
      const int n_passive_scalars = flux.passive_scalars.size();
      for (int i = 0; i < n_passive_scalars; ++i) {
        const double rhoYL = left_.passive_scalars[i];
        const double rhoYR = right_.passive_scalars[i];
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const double f_rho = flux.density;
          const int upwind = (f_rho > 0.0);
          const double f_hllem = f_rho * (upwind * YL - (1 - upwind) * YR);
          flux.passive_scalars[i] = f_hllem;
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          const double b_delta_alpha_i = b * delta * alpha_i;
          flux.passive_scalars[i] = flux_hlle_.passive_scalars[i] -
                                    b_delta_alpha_2 * roeY - b_delta_alpha_i;
        }
      }
    }
  } else {
    static_assert(Dim == 3);
    const int ix = int(dir);
    const int iy = (ix + 1) % 3;
    const int iz = (iy + 1) % 3;
    const double l21 = gm1 / roeA2 * (roeH - Msq * roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * Msq * roeU[0];
    const double l23 = gm1 / roeA2 * roeU[1];
    const double l24 = gm1 / roeA2 * roeU[2];
    const double l25 = gm1 / roeA2 * (-1.0);
    const double alpha_2 = l21 * deltaRho + l22 * deltaRhoU[0] +
                           l23 * deltaRhoU[1] + l24 * deltaRhoU[2] +
                           l25 * deltaRhoE;
    const double b_delta = b * delta;
    const double b_delta_alpha_2 = b_delta * alpha_2;
    if constexpr (Larrouturou) {
      const double f_rho = flux.density;
      const int upwind = (f_rho > 0.0);
      const double f_hllem_rhov =
          f_rho * (upwind * vL_[iy] - (1 - upwind) * vR_[iy]);
      const double f_hllem_rhow =
          f_rho * (upwind * vL_[iz] - (1 - upwind) * vR_[iz]);
      const double f_hllem_rhoE =
          f_rho * (upwind * eKinRL - (1 - upwind) * eKinRR);
      flux.momentum[iy] = f_hllem_rhov;
      flux.momentum[iz] = f_hllem_rhow;
      flux.energy = flux_hlle_.energy - b_delta_alpha_2 * squaredNormRoeU_half +
                    f_hllem_rhoE;
    } else {
      const double l31 = -roeU[iy];
      const double l33 = 1;
      const double l41 = -roeU[iz];
      const double l44 = 1;
      const double alpha_3 = l31 * deltaRho + l33 * deltaRhoU[iy];
      const double alpha_4 = l41 * deltaRho + l44 * deltaRhoU[iz];
      const double b_delta_alpha_3 = b_delta * alpha_3;
      const double b_delta_alpha_4 = b_delta * alpha_4;
      flux.momentum[iy] = flux_hlle_.momentum[iy] - b_delta_alpha_2 * roeU[iy] -
                          b_delta_alpha_3;
      flux.momentum[iz] = flux_hlle_.momentum[iz] - b_delta_alpha_2 * roeU[iz] -
                          b_delta_alpha_4;
      flux.energy = flux_hlle_.energy - b_delta_alpha_2 * squaredNormRoeU_half -
                    b_delta_alpha_3 * roeU[iy] - b_delta_alpha_4 * roeU[iz];
    }
    if constexpr (fub::euler::state_with_species<
                                                 Conservative>()) {
      const int n_species = flux.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left_, i);
        const double rhoYR = fub::euler::Species(equation_, right_, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const double f_rho = flux.density;
          const int upwind = (f_rho > 0.0);
          const double f_hllem = f_rho * (upwind * YL - (1 - upwind) * YR);
          flux.species[i] = f_hllem;
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          const double b_delta_alpha_i = b * delta * alpha_i;
          flux.species[i] =
              flux_hlle_.species[i] - b_delta_alpha_2 * roeY - b_delta_alpha_i;
        }
      }
    }
    if constexpr (fub::euler::state_with_passive_scalars<Conservative>()) {
      const int n_passive_scalars = flux.passive_scalars.size();
      for (int i = 0; i < n_passive_scalars; ++i) {
        const double rhoYL = left_.passive_scalars[i];
        const double rhoYR = right_.passive_scalars[i];
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const double f_rho = flux.density;
          const int upwind = (f_rho > 0.0);
          const double f_hllem = f_rho * (upwind * YL - (1 - upwind) * YR);
          flux.passive_scalars[i] = f_hllem;
        } else {
          const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const double deltaRhoY = rhoYR - rhoYL;
          const double li1 = -roeY;
          const double li4 = 1.0;
          const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
          const double b_delta_alpha_i = b * delta * alpha_i;
          flux.passive_scalars[i] = flux_hlle_.passive_scalars[i] -
                                    b_delta_alpha_2 * roeY - b_delta_alpha_i;
        }
      }
    }
  }
  const double p_star = (bR * pL - bL * pR) / db_positive;
  return p_star;
}

template <typename EulerEquation, bool Larrouturou>
Array1d Hllem<EulerEquation, Larrouturou>::ComputeNumericFlux(
    ConservativeArray& flux, span<const CompleteArray, 2> states,
    Duration /* dt */, double /* dx */, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();

  left_array_ = states[0];
  right_array_ = states[1];
  const Array1d Msq = Array1d::Constant(euler::MaSq(equation_));
  const Array1d Minv = 1.0 / Msq.sqrt();
  const Array1d rhoL = fub::euler::Density(equation_, left_array_);
  const Array1d rhoR = fub::euler::Density(equation_, right_array_);
  [[maybe_unused]] Array<double, Dim> vL_;
  [[maybe_unused]] Array<double, Dim> vR_;
  [[maybe_unused]] Array1d eKinRL = Array1d::Zero();
  [[maybe_unused]] Array1d eKinRR = Array1d::Zero();
  // If we do the Larrouturou version with velocity we need to set the radial
  // velocity to zero.
  if constexpr (Larrouturou) {
    std::array<int, Dim> dim_mask{};
    int d = static_cast<int>(dir);
    for (int i = 0; i < Dim; ++i) {
      dim_mask[i] = i == d;
      vL_.row(i) = left_array_.momentum.row(i) / rhoL;
      vR_.row(i) = right_array_.momentum.row(i) / rhoR;
      eKinRL += (1 - dim_mask[i]) * vL_.row(i) * vL_.row(i);
      eKinRR += (1 - dim_mask[i]) * vR_.row(i) * vR_.row(i);
      left_array_.momentum.row(i) *= dim_mask[i];
      right_array_.momentum.row(i) *= dim_mask[i];
    }
    eKinRL *= 0.5 * Msq;
    eKinRR *= 0.5 * Msq;
    left_array_.energy -= rhoL * eKinRL;
    right_array_.energy -= rhoR * eKinRR;
  }

  fub::Flux(equation_, fluxL_array_, left_array_, dir);
  fub::Flux(equation_, fluxR_array_, right_array_, dir);

  // Compute Einfeldt signals velocitie

  const Array1d aL = euler::SpeedOfSound(equation_, left_array_);
  const Array1d aR = euler::SpeedOfSound(equation_, right_array_);
  const Array1d rhoEL = euler::Energy(equation_, left_array_);
  const Array1d rhoER = euler::Energy(equation_, right_array_);
  const Array1d pL = euler::Pressure(equation_, left_array_);
  const Array1d pR = euler::Pressure(equation_, right_array_);
  const Array1d hL = (rhoEL + pL) / rhoL;
  const Array1d hR = (rhoER + pR) / rhoR;
  const Array1d sqRhoL = rhoL.sqrt();
  const Array1d sqRhoR = rhoR.sqrt();
  const Array1d sqRhoSum = sqRhoL + sqRhoR;
  const Array1d sqRhoL_over_Sum = sqRhoL / sqRhoSum;
  const Array1d sqRhoR_over_Sum = sqRhoR / sqRhoSum;

  Array<double, Dim> uL;
  Array<double, Dim> uR;
  Array<double, Dim> roeU;
  for (int i = 0; i < Dim; ++i) {
    uL.row(i) = left_array_.momentum.row(i) / rhoL;
    uR.row(i) = right_array_.momentum.row(i) / rhoR;
    roeU.row(i) = sqRhoL_over_Sum * uL.row(i) + sqRhoR_over_Sum * uR.row(i);
  }
  const Array1d roeH = sqRhoL_over_Sum * hL + sqRhoR_over_Sum * hR;

  Array1d squaredNormRoeU = Array1d::Zero();
  for (int i = 0; i < Dim; ++i) {
    squaredNormRoeU += roeU.row(i) * roeU.row(i);
  }
  squaredNormRoeU *= Msq;
  const Array1d squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  const Array1d gammaL = fub::euler::Gamma(equation_, left_array_);
  const Array1d gammaR = fub::euler::Gamma(equation_, right_array_);
  const Array1d roeGamma = (sqRhoL * gammaL + sqRhoR * gammaR) / sqRhoSum;
  const Array1d gm1 = roeGamma - 1.0;

  const Array1d roeA2 = gm1 * (roeH - squaredNormRoeU_half);
  const Array1d roeA = roeA2.sqrt();
  const Array1d maxA = aL.max(aR);
  const Array1d sL1 = uL.row(int(dir)) - maxA * Minv;
  const Array1d sL2 = roeU.row(int(dir)) - roeA * Minv;
  const Array1d sR1 = roeU.row(int(dir)) + roeA * Minv;
  const Array1d sR2 = uR.row(int(dir)) + maxA * Minv;

  const Array1d sL = sL1.min(sL2);
  const Array1d sR = sR1.max(sR2);

  const Array1d zero = Array1d::Zero();
  const Array1d bL = sL.min(zero);
  const Array1d bR = sR.max(zero);
  const Array1d bLbR = bL * bR;
  const Array1d db = bR - bL;
  const Array1d db_positive = (db > 0.0).select(db, Array1d::Constant(1.0));

  ForEachComponent(
      [&](auto&& nf, Array1d fL, Array1d fR, Array1d qL, Array1d qR) {
        nf = (bR * fL - bL * fR + bLbR * (qR - qL)) / db_positive;
      },
      flux_hlle_array_, fluxL_array_, fluxR_array_, AsCons(left_array_),
      AsCons(right_array_));

  const Array1d u_bar = 0.5 * (sR + sL);
  const Array1d u_bar_abs = u_bar.abs();
  const Array1d delta = roeA / (roeA + u_bar_abs);
  const Array1d b = bLbR / db_positive;

  const Array1d deltaRho = rhoR - rhoL;
  const Array<double, Dim> deltaRhoU =
      right_array_.momentum - left_array_.momentum;
  const Array1d deltaRhoE = rhoER - rhoEL;

  const Array1d gm1_over_roeA2 = gm1 / roeA2;

  if constexpr (Dim == 1) {
    const Array1d l21 = gm1_over_roeA2 * (roeH - squaredNormRoeU);
    const Array1d l22 = gm1_over_roeA2 * Msq * roeU;
    const Array1d l23 = gm1_over_roeA2 * (-1.0);
    const Array1d alpha_2 = l21 * deltaRho + l22 * deltaRhoU + l23 * deltaRhoE;
    const Array1d b_delta_alpha_2 = b * delta * alpha_2;

    flux.density = flux_hlle_array_.density - b_delta_alpha_2;
    flux.momentum = flux_hlle_array_.momentum - b_delta_alpha_2 * roeU;
    flux.energy =
        flux_hlle_array_.energy - b_delta_alpha_2 * squaredNormRoeU_half;
    if constexpr (euler::state_with_species<Conservative>()) {
      const int n_species = flux.species.rows();
      for (int i = 0; i < n_species; ++i) {
        const Array1d rhoYL = fub::euler::Species(equation_, left_array_, i);
        const Array1d rhoYR = fub::euler::Species(equation_, right_array_, i);
        const Array1d YL = rhoYL / rhoL;
        const Array1d YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const Array1d f_rho = flux.density;
          const MaskArray upwind = (f_rho > 0.0);
          const Array1d f_hllem = f_rho * upwind.select(YL, YR);
          flux.species.row(i) = f_hllem;
        } else {
          const Array1d roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const Array1d deltaRhoY = rhoYR - rhoYL;
          const Array1d li1 = -roeY;
          const Array1d alpha_i = li1 * deltaRho + deltaRhoY;
          const Array1d b_delta_alpha_i = b * delta * alpha_i;
          flux.species.row(i) = flux_hlle_array_.species.row(i) -
                                b_delta_alpha_2 * roeY - b_delta_alpha_i;
        }
      }
    }
    if constexpr (euler::state_with_passive_scalars<Conservative>()) {
      const int n_passive_scalars = flux.passive_scalars.rows();
      for (int i = 0; i < n_passive_scalars; ++i) {
        const Array1d rhoYL = left_array_.passive_scalars.row(i);
        const Array1d rhoYR = right_array_.passive_scalars.row(i);
        const Array1d YL = rhoYL / rhoL;
        const Array1d YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const Array1d f_rho = flux.density;
          const MaskArray upwind = (f_rho > 0.0);
          const Array1d f_hllem = f_rho * upwind.select(YL, YR);
          flux.passive_scalars.row(i) = f_hllem;
        } else {
          const Array1d roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const Array1d deltaRhoY = rhoYR - rhoYL;
          const Array1d li1 = -roeY;
          const Array1d alpha_i = li1 * deltaRho + deltaRhoY;
          const Array1d b_delta_alpha_i = b * delta * alpha_i;
          flux.passive_scalars.row(i) =
              flux_hlle_array_.passive_scalars.row(i) - b_delta_alpha_2 * roeY -
              b_delta_alpha_i;
        }
      }
    }
    FUB_ASSERT(!flux.density.isNaN().any());
    FUB_ASSERT(!flux.momentum.isNaN().any());
    FUB_ASSERT(!flux.energy.isNaN().any());
    if constexpr (euler::state_with_species<Conservative>()) {
      FUB_ASSERT(!flux.species.isNaN().any());
    }
    if constexpr (euler::state_with_species<Conservative>()) {
      FUB_ASSERT(!flux.passive_scalars.isNaN().any());
    }
  } else if constexpr (Dim == 2) {
    const int ix = int(dir);
    const int iy = (ix == 0);
    const Array1d l21 = gm1_over_roeA2 * (roeH - squaredNormRoeU);
    const Array1d l22 = gm1_over_roeA2 * Msq * roeU.row(0);
    const Array1d l23 = gm1_over_roeA2 * roeU.row(1);
    const Array1d l24 = -gm1_over_roeA2;
    const Array1d alpha_2 = l21 * deltaRho + l22 * deltaRhoU.row(0) +
                            l23 * deltaRhoU.row(1) + l24 * deltaRhoE;
    const Array1d b_delta = b * delta;
    const Array1d b_delta_alpha_2 = b_delta * alpha_2;
    flux.density = flux_hlle_array_.density - b_delta_alpha_2;
    flux.momentum.row(ix) =
        flux_hlle_array_.momentum.row(ix) - b_delta_alpha_2 * roeU.row(ix);
    if constexpr (Larrouturou) {
      const MaskArray upwind = (flux.density > 0.0);
      const Array1d f_hllem_rhov =
          flux.density * upwind.select(vL_.row(iy), vR_.row(iy));
      const Array1d f_hllem_rhoE = flux.density * upwind.select(eKinRL, eKinRR);
      flux.momentum.row(iy) = f_hllem_rhov;
      flux.energy = flux_hlle_array_.energy -
                    b_delta_alpha_2 * squaredNormRoeU_half + f_hllem_rhoE;
    } else {
      const Array1d l31 = -roeU.row(iy);
      const Array1d alpha_3 = l31 * deltaRho + deltaRhoU.row(iy);
      const Array1d b_delta_alpha_3 = b_delta * alpha_3;
      flux.momentum.row(iy) = flux_hlle_array_.momentum.row(iy) -
                              b_delta_alpha_2 * roeU.row(iy) - b_delta_alpha_3;
      flux.energy = flux_hlle_array_.energy -
                    b_delta_alpha_2 * squaredNormRoeU_half -
                    b_delta_alpha_3 * roeU.row(iy);
    }
    if constexpr (euler::state_with_species<Conservative>()) {
      const int n_species = flux.species.rows();
      for (int i = 0; i < n_species; ++i) {
        const Array1d rhoYL = fub::euler::Species(equation_, left_array_, i);
        const Array1d rhoYR = fub::euler::Species(equation_, right_array_, i);
        const Array1d YL = rhoYL / rhoL;
        const Array1d YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const Array1d f_rho = flux.density;
          const MaskArray upwind = (f_rho > 0.0);
          const Array1d f_hllem = f_rho * upwind.select(YL, YR);
          flux.species.row(i) = f_hllem;
        } else {
          const Array1d roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const Array1d deltaRhoY = rhoYR - rhoYL;
          const Array1d li1 = -roeY;
          const Array1d alpha_i = li1 * deltaRho + deltaRhoY;
          const Array1d b_delta_alpha_i = b_delta * alpha_i;
          // MaskArray mask = rhoL > 0.0 && rhoR > 0.0;
          flux.species.row(i) = flux_hlle_array_.species.row(i) -
                                b_delta_alpha_2 * roeY - b_delta_alpha_i;
          // MaskArray finite_mask = !flux.species.row(i).isNaN();
          // FUB_ASSERT((!mask || finite_mask).all());
        }
      }
    }
    if constexpr (euler::state_with_passive_scalars<Conservative>()) {
      const int n_passive_scalars = flux.passive_scalars.rows();
      for (int i = 0; i < n_passive_scalars; ++i) {
        const Array1d rhoYL = left_array_.passive_scalars.row(i);
        const Array1d rhoYR = right_array_.passive_scalars.row(i);
        const Array1d YL = rhoYL / rhoL;
        const Array1d YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const Array1d f_rho = flux.density;
          const MaskArray upwind = (f_rho > 0.0);
          const Array1d f_hllem = f_rho * upwind.select(YL, YR);
          flux.passive_scalars.row(i) = f_hllem;
        } else {
          const Array1d roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const Array1d deltaRhoY = rhoYR - rhoYL;
          const Array1d li1 = -roeY;
          const Array1d alpha_i = li1 * deltaRho + deltaRhoY;
          const Array1d b_delta_alpha_i = b * delta * alpha_i;
          flux.passive_scalars.row(i) =
              flux_hlle_array_.passive_scalars.row(i) - b_delta_alpha_2 * roeY -
              b_delta_alpha_i;
        }
      }
    }
  } else {
    static_assert(Dim == 3);
    const int ix = int(dir);
    const int iy = (ix + 1) % 3;
    const int iz = (iy + 1) % 3;
    const Array1d l21 = gm1_over_roeA2 * (roeH - squaredNormRoeU);
    const Array1d l22 = gm1_over_roeA2 * Msq * roeU.row(0);
    const Array1d l23 = gm1_over_roeA2 * roeU.row(1);
    const Array1d l24 = gm1_over_roeA2 * roeU.row(2);
    const Array1d l25 = -gm1_over_roeA2;

    const Array1d alpha_2 = l21 * deltaRho + l22 * deltaRhoU.row(0) +
                            l23 * deltaRhoU.row(1) + l24 * deltaRhoU.row(2) +
                            l25 * deltaRhoE;
    const Array1d b_delta = b * delta;
    const Array1d b_delta_alpha_2 = b_delta * alpha_2;
    flux.density = flux_hlle_array_.density - b_delta_alpha_2;
    flux.momentum.row(ix) =
        flux_hlle_array_.momentum.row(ix) - b_delta_alpha_2 * roeU.row(iy);

    if constexpr (Larrouturou) {
      const MaskArray upwind = (flux.density > 0.0);
      const Array1d f_hllem_rhov =
          flux.density * upwind.select(vL_.row(iy), vR_.row(iy));
      const Array1d f_hllem_rhow =
          flux.density * upwind.select(vL_.row(iz), vR_.row(iz));
      const Array1d f_hllem_rhoE = flux.density * upwind.select(eKinRL, eKinRR);
      flux.momentum.row(iy) = f_hllem_rhov;
      flux.momentum.row(iz) = f_hllem_rhow;
      flux.energy = flux_hlle_array_.energy -
                    b_delta_alpha_2 * squaredNormRoeU_half + f_hllem_rhoE;
    } else {
      const Array1d l31 = -roeU.row(iy);
      const Array1d l41 = -roeU.row(iz);
      const Array1d alpha_3 = l31 * deltaRho + deltaRhoU.row(iy);
      const Array1d alpha_4 = l41 * deltaRho + deltaRhoU.row(iz);
      const Array1d b_delta_alpha_3 = b_delta * alpha_3;
      const Array1d b_delta_alpha_4 = b_delta * alpha_4;
      flux.momentum.row(iy) = flux_hlle_array_.momentum.row(iy) -
                              b_delta_alpha_2 * roeU.row(iy) - b_delta_alpha_3;
      flux.momentum.row(iz) = flux_hlle_array_.momentum.row(iz) -
                              b_delta_alpha_2 * roeU.row(iz) - b_delta_alpha_4;
      flux.energy =
          flux_hlle_array_.energy - b_delta_alpha_2 * squaredNormRoeU_half -
          b_delta_alpha_3 * roeU.row(iy) - b_delta_alpha_4 * roeU.row(iz);
    }

    if constexpr (euler::state_with_species<Conservative>()) {
      const int n_species = flux.species.rows();
      for (int i = 0; i < n_species; ++i) {
        const Array1d rhoYL = fub::euler::Species(equation_, left_array_, i);
        const Array1d rhoYR = fub::euler::Species(equation_, right_array_, i);
        const Array1d YL = rhoYL / rhoL;
        const Array1d YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const Array1d f_rho = flux.density;
          const MaskArray upwind = (f_rho > 0.0);
          const Array1d f_hllem = f_rho * upwind.select(YL, YR);
          flux.species.row(i) = f_hllem;
        } else {
          const Array1d roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const Array1d deltaRhoY = rhoYR - rhoYL;
          const Array1d li1 = -roeY;
          const Array1d alpha_i = li1 * deltaRho + deltaRhoY;
          const Array1d b_delta_alpha_i = b_delta * alpha_i;
          flux.species.row(i) = flux_hlle_array_.species.row(i) -
                                b_delta_alpha_2 * roeY - b_delta_alpha_i;
        }
      }
    }
    if constexpr (euler::state_with_passive_scalars<Conservative>()) {
      const int n_passive_scalars = flux.passive_scalars.rows();
      for (int i = 0; i < n_passive_scalars; ++i) {
        const Array1d rhoYL = left_array_.passive_scalars.row(i);
        const Array1d rhoYR = right_array_.passive_scalars.row(i);
        const Array1d YL = rhoYL / rhoL;
        const Array1d YR = rhoYR / rhoR;
        if constexpr (Larrouturou) {
          const Array1d f_rho = flux.density;
          const MaskArray upwind = (f_rho > 0.0);
          const Array1d f_hllem = f_rho * upwind.select(YL, YR);
          flux.passive_scalars.row(i) = f_hllem;
        } else {
          const Array1d roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
          const Array1d deltaRhoY = rhoYR - rhoYL;
          const Array1d li1 = -roeY;
          const Array1d alpha_i = li1 * deltaRho + deltaRhoY;
          const Array1d b_delta_alpha_i = b * delta * alpha_i;
          flux.passive_scalars.row(i) =
              flux_hlle_array_.passive_scalars.row(i) - b_delta_alpha_2 * roeY -
              b_delta_alpha_i;
        }
      }
    }
  }

  const Array1d p_star = (bR * pL - bL * pR) / db_positive;
  return p_star;
}

template <typename EulerEquation, bool Larrouturou>
void Hllem<EulerEquation, Larrouturou>::ComputeNumericFlux(
    ConservativeArray& flux, Array1d face_fractions,
    span<const CompleteArray, 2> states,
    span<const Array1d, 2> /* volume_fractions */, Duration dt, double dx,
    Direction dir) {
  ComputeNumericFlux(flux, states, dt, dx, dir);
  MaskArray face_mask = face_fractions > 0.0;
  ForEachComponent([&face_mask](auto&& f) { f = face_mask.select(f, 0.0); },
                   flux);
  FUB_ASSERT(!flux.density.isNaN().any());
  FUB_ASSERT(!flux.momentum.isNaN().any());
  FUB_ASSERT(!flux.energy.isNaN().any());
  if constexpr (euler::state_with_species<Conservative>()) {
    FUB_ASSERT(!flux.species.isNaN().any());
  }
  if constexpr (euler::state_with_passive_scalars<Conservative>()) {
    FUB_ASSERT(!flux.passive_scalars.isNaN().any());
  }
}

template <typename EulerEquation, bool Larrouturou>
double Hllem<EulerEquation, Larrouturou>::ComputeStableDt(
    span<const Complete, 2> states, double dx, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  const Complete& left = states[0];
  const Complete& right = states[1];
  const double Msq = euler::MaSq(equation_);
  const double Minv = 1.0 / std::sqrt(Msq);
  const double gm1 = equation_.gamma - 1.0;

  // Compute Einfeldt signals velocities

  int d = static_cast<int>(dir);

  const double rhoL = left.density;
  const double rhoR = right.density;
  const double rhoUL = left.momentum[d];
  const double rhoUR = right.momentum[d];
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  const double aL = left.speed_of_sound;
  const double aR = right.speed_of_sound;
  const double rhoEL = left.energy;
  const double rhoER = right.energy;
  const double pL = left.pressure;
  const double pR = right.pressure;
  const double hL = (rhoEL + pL) / rhoL;
  const double hR = (rhoER + pR) / rhoR;
  const double sqRhoL = std::sqrt(rhoL);
  const double sqRhoR = std::sqrt(rhoR);
  const double sqRhoSum = sqRhoL + sqRhoR;
  Array<double, 1, Dim> vL = left.momentum / rhoL;
  Array<double, 1, Dim> vR = right.momentum / rhoR;
  FUB_ASSERT(sqRhoSum > 0.0);
  Array<double, 1, Dim> roeU = (sqRhoL * vL + sqRhoR * vR) / sqRhoSum;
  const double roeU0 = roeU[d];
  const double roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;
  const double roeA2 = gm1 * (roeH - 0.5 * Msq * roeU.matrix().squaredNorm());
  const double roeA = std::sqrt(roeA2);
  const double maxA = std::max(aL, aR);
  const double sL1 = uL - maxA * Minv;
  const double sL2 = roeU0 - roeA * Minv;
  const double sR1 = roeU0 + roeA * Minv;
  const double sR2 = uR + maxA * Minv;

  const double sL = std::min(sL1, sL2);
  const double sR = std::max(sR1, sR2);

  const double sL_abs = std::abs(sL);
  const double sR_abs = std::abs(sR);

  const double maxS = std::max(sL_abs, sR_abs);
  const double max_dt = dx / maxS;

  return max_dt;
}

template <typename EulerEquation, bool Larrouturou>
Array1d Hllem<EulerEquation, Larrouturou>::ComputeStableDt(
    span<const CompleteArray, 2> states, double dx, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  const Array1d gm1 = (equation_.gamma_array_ - Array1d::Constant(1.0));
  const Array1d Msq = Array1d::Constant(euler::MaSq(equation_));
  const Array1d Minv = 1.0 / Msq.sqrt();
  // Compute Einfeldt signals velocities

  const Array1d rhoL = left.density;
  const Array1d rhoR = right.density;
  const Array1d aL = left.speed_of_sound;
  const Array1d aR = right.speed_of_sound;
  const Array1d rhoEL = left.energy;
  const Array1d rhoER = right.energy;
  const Array1d pL = left.pressure;
  const Array1d pR = right.pressure;
  const Array1d hL = (rhoEL + pL) / rhoL;
  const Array1d hR = (rhoER + pR) / rhoR;
  const Array1d sqRhoL = rhoL.sqrt();
  const Array1d sqRhoR = rhoR.sqrt();
  const Array1d sqRhoSum = sqRhoL + sqRhoR;
  const Array1d sqRhoL_over_Sum = sqRhoL / sqRhoSum;
  const Array1d sqRhoR_over_Sum = sqRhoR / sqRhoSum;

  Array<double, Dim> uL;
  Array<double, Dim> uR;
  Array<double, Dim> roeU;
  for (int i = 0; i < Dim; ++i) {
    uL.row(i) = left.momentum.row(i) / rhoL;
    uR.row(i) = right.momentum.row(i) / rhoR;
    roeU.row(i) = sqRhoL_over_Sum * uL.row(i) + sqRhoR_over_Sum * uR.row(i);
  }
  const Array1d roeH = sqRhoL_over_Sum * hL + sqRhoR_over_Sum * hR;

  Array1d squaredNormRoeU = Array1d::Zero();
  for (int i = 0; i < Dim; ++i) {
    squaredNormRoeU += roeU.row(i) * roeU.row(i);
  }
  squaredNormRoeU *= Msq;
  const Array1d squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  const Array1d roeA2 = gm1 * (roeH - squaredNormRoeU_half);
  const Array1d roeA = roeA2.sqrt();
  const Array1d maxA = aL.max(aR);

  const Array1d sL1 = uL.row(int(dir)) - maxA * Minv;
  const Array1d sL2 = roeU.row(int(dir)) - roeA * Minv;
  const Array1d sR1 = roeU.row(int(dir)) + roeA * Minv;
  const Array1d sR2 = uR.row(int(dir)) + maxA * Minv;

  const Array1d sL = sL1.min(sL2);
  const Array1d sR = sR1.max(sR2);

  const Array1d sL_abs = sL.abs();
  const Array1d sR_abs = sR.abs();

  const Array1d maxS = sL_abs.max(sR_abs);
  const Array1d max_dt = Array1d::Constant(dx) / maxS;

  return max_dt;
}

template <typename EulerEquation, bool Larrouturou>
Array1d Hllem<EulerEquation, Larrouturou>::ComputeStableDt(
    span<const CompleteArray, 2> states, Array1d face_fraction,
    span<const Array1d, 2>, double dx, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  MaskArray face_mask = face_fraction > 0.0;

  const Array1d gm1 = (equation_.gamma_array_ - Array1d::Constant(1.0));
  const Array1d Msq = Array1d::Constant(euler::MaSq(equation_));
  const Array1d Minv = 1.0 / Msq.sqrt();

  // Compute Einfeldt signals velocities

  const Array1d ones = Array1d::Constant(1.0);
  const Array1d zeros = Array1d::Constant(0.0);

  const Array1d rhoL = face_mask.select(left.density, ones);
  const Array1d rhoR = face_mask.select(right.density, ones);
  const Array1d aL = face_mask.select(left.speed_of_sound, zeros);
  const Array1d aR = face_mask.select(right.speed_of_sound, zeros);
  const Array1d rhoEL = face_mask.select(left.energy, zeros);
  const Array1d rhoER = face_mask.select(right.energy, zeros);
  const Array1d pL = face_mask.select(left.pressure, zeros);
  const Array1d pR = face_mask.select(right.pressure, zeros);
  const Array1d hL = (rhoEL + pL) / rhoL;
  const Array1d hR = (rhoER + pR) / rhoR;
  const Array1d sqRhoL = rhoL.sqrt();
  const Array1d sqRhoR = rhoR.sqrt();
  const Array1d sqRhoSum = sqRhoL + sqRhoR;
  const Array1d sqRhoL_over_Sum = sqRhoL / sqRhoSum;
  const Array1d sqRhoR_over_Sum = sqRhoR / sqRhoSum;

  Array<double, Dim> uL;
  Array<double, Dim> uR;
  Array<double, Dim> roeU;
  for (int i = 0; i < Dim; ++i) {
    uL.row(i) = left.momentum.row(i) / rhoL;
    uR.row(i) = right.momentum.row(i) / rhoR;
    roeU.row(i) = sqRhoL_over_Sum * uL.row(i) + sqRhoR_over_Sum * uR.row(i);
  }
  const Array1d roeH = sqRhoL_over_Sum * hL + sqRhoR_over_Sum * hR;

  Array1d squaredNormRoeU = Array1d::Zero();
  for (int i = 0; i < Dim; ++i) {
    squaredNormRoeU += roeU.row(i) * roeU.row(i);
  }
  squaredNormRoeU *= Msq;
  const Array1d squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  const Array1d roeA2 = gm1 * (roeH - squaredNormRoeU_half);
  const Array1d roeA = roeA2.sqrt();
  const Array1d maxA = aL.max(aR);

  const Array1d sL1 = uL.row(int(dir)) - maxA * Minv;
  const Array1d sL2 = roeU.row(int(dir)) - roeA * Minv;
  const Array1d sR1 = roeU.row(int(dir)) + roeA * Minv;
  const Array1d sR2 = uR.row(int(dir)) + maxA * Minv;

  const Array1d sL = sL1.min(sL2);
  const Array1d sR = sR1.max(sR2);

  const Array1d sL_abs = sL.abs();
  const Array1d sR_abs = sR.abs();

  const Array1d maxS = sL_abs.max(sR_abs);
  const Array1d maxS_positive = face_mask.select(maxS, ones);

  const Array1d max_dt = Array1d::Constant(dx) / maxS_positive;
  const Array1d dummy_value =
      Array1d::Constant(std::numeric_limits<double>::max());
  const Array1d max_dt_result = face_mask.select(max_dt, dummy_value);
  return max_dt_result;
}

} // namespace fub::perfect_gas