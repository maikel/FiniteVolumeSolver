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

template <typename EulerEquation>
void Hllem<EulerEquation>::SolveRiemannProblem(Complete& solution,
                                               const Complete& left,
                                               const Complete& right,
                                               Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  fub::Flux(equation_, fluxL_, left, dir);
  fub::Flux(equation_, fluxR_, right, dir);

  // Compute Einfeldt signals velocities

  int d = static_cast<int>(dir);

  const double rhoL = fub::euler::Density(equation_, left);
  const double rhoR = fub::euler::Density(equation_, right);
  const double rhoUL = fub::euler::Momentum(equation_, left, d);
  const double rhoUR = fub::euler::Momentum(equation_, right, d);
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  const double aL = fub::euler::SpeedOfSound(equation_, left);
  const double aR = fub::euler::SpeedOfSound(equation_, right);
  const double rhoEL = fub::euler::Energy(equation_, left);
  const double rhoER = fub::euler::Energy(equation_, right);
  const double pL = fub::euler::Pressure(equation_, left);
  const double pR = fub::euler::Pressure(equation_, right);
  const double hL = (rhoEL + pL) / rhoL;
  const double hR = (rhoER + pR) / rhoR;
  const double sqRhoL = std::sqrt(rhoL);
  const double sqRhoR = std::sqrt(rhoR);
  const double sqRhoSum = sqRhoL + sqRhoR;
  Array<double, 1, Dim> vL = fub::euler::Velocity(equation_, left);
  Array<double, 1, Dim> vR = fub::euler::Velocity(equation_, right);
  FUB_ASSERT(sqRhoSum > 0.0);
  Array<double, 1, Dim> roeU = (sqRhoL * vL + sqRhoR * vR) / sqRhoSum;
  const double roeU0 = roeU[d];
  const double roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;

  const double gammaL = fub::euler::Gamma(equation_, left);
  const double gammaR = fub::euler::Gamma(equation_, right);
  const double roeGamma = (sqRhoL * gammaL + sqRhoR * gammaR) / sqRhoSum;
  const double gm1 = roeGamma - 1.0;
  const double beta = gm1 / (2 * roeGamma);

  const double roeA2 = gm1 * (roeH - 0.5 * roeU.matrix().squaredNorm());
  const double roeA = std::sqrt(roeA2);
  const double sL1 = uL - beta * aL;
  const double sL2 = roeU0 - roeA;
  const double sR1 = roeU0 + roeA;
  const double sR2 = uR + beta * aR;

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
      w_hlle_, fluxL_, fluxR_, AsCons(left), AsCons(right));

  const double squaredNormRoeU = roeU.matrix().squaredNorm();
  const double squaredNormRoeU_half = 0.5 * squaredNormRoeU;
  const double u_bar = 0.5 * (sR + sL);
  const double u_bar_abs = std::abs(u_bar);
  const double delta = roeA / (roeA + u_bar_abs);

  const double deltaRho = rhoR - rhoL;
  const Array<double, 1, Dim> deltaRhoU = right.momentum - left.momentum;
  const double deltaRhoE = rhoER - rhoEL;

  if constexpr (Dim == 1) {
    const double l21 = gm1 / roeA2 * (roeH - roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * roeU[0];
    const double l23 = gm1 / roeA2 * (-1.0);
    const double alpha_2 =
        l21 * deltaRho + l22 * deltaRhoU[0] + l23 * deltaRhoE;
    const double u_bar_delta_alpha_2 = u_bar * delta * alpha_2;
    w_hllem_.density = w_hlle_.density - u_bar_delta_alpha_2 * 1.0;
    w_hllem_.momentum = w_hlle_.momentum - u_bar_delta_alpha_2 * roeU[0];
    w_hllem_.energy =
        w_hlle_.energy - u_bar_delta_alpha_2 * squaredNormRoeU_half;
    if constexpr (fub::euler::state_with_species<EulerEquation,
                                                 Conservative>()) {
      const int n_species = w_hllem_.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left, i);
        const double rhoYR = fub::euler::Species(equation_, right, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
        const double deltaRhoY = rhoYR - rhoYL;
        const double li1 = -roeY;
        const double li4 = 1.0;
        const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
        w_hllem_.species[i] = w_hlle_.species[i] - u_bar_delta_alpha_2 * roeY -
                              u_bar * delta * alpha_i;
      }
    }
  } else if constexpr (Dim == 2) {
    const int ix = int(dir);
    const int iy = (ix + 1) % 2;
    const double l21 = gm1 / roeA2 * (roeH - roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * roeU[0];
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
    w_hllem_.momentum[iy] = w_hlle_.momentum[iy] -
                            u_bar_delta_alpha_2 * roeU[iy] -
                            u_bar_delta_alpha_3 * 1.0;
    w_hllem_.energy = w_hlle_.energy -
                      u_bar_delta_alpha_2 * squaredNormRoeU_half -
                      u_bar_delta_alpha_3 * roeU[iy];
    if constexpr (fub::euler::state_with_species<EulerEquation,
                                                 Conservative>()) {
      const int n_species = w_hllem_.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left, i);
        const double rhoYR = fub::euler::Species(equation_, right, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
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
  } else {
    static_assert(Dim == 3);
    const int ix = int(dir);
    const int iy = (ix + 1) % 3;
    const int iz = (iy + 1) % 3;
    const double l21 = gm1 / roeA2 * (roeH - roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * roeU[0];
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
    w_hllem_.momentum[iy] = w_hlle_.momentum[iy] -
                            u_bar_delta * (alpha_2 * roeU[iy] + alpha_3 * 1.0);
    w_hllem_.momentum[iz] = w_hlle_.momentum[iz] -
                            u_bar_delta * (alpha_2 * roeU[iz] + alpha_4 * 1.0);
    w_hllem_.energy = w_hlle_.energy -
                      u_bar_delta * (alpha_2 * squaredNormRoeU_half +
                                     alpha_3 * roeU[iy] + alpha_4 * roeU[iz]);
    if constexpr (fub::euler::state_with_species<EulerEquation,
                                                 Conservative>()) {
      const int n_species = w_hllem_.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left, i);
        const double rhoYR = fub::euler::Species(equation_, right, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
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

  if (0.0 < sL) {
    solution = left;
  } else if (sR < 0.0) {
    solution = right;
  } else {
    AsCons(solution) = w_hllem_;
    CompleteFromCons(equation_, solution, solution);
  }
}

template <typename EulerEquation>
void Hllem<EulerEquation>::ComputeNumericFlux(Conservative& flux,
                                              span<const Complete, 2> states,
                                              Duration /* dt */,
                                              double /* dx */, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  const Complete& left = states[0];
  const Complete& right = states[1];

  fub::Flux(equation_, fluxL_, left, dir);
  fub::Flux(equation_, fluxR_, right, dir);

  // Compute Einfeldt signals velocities

  int d = static_cast<int>(dir);

  const double rhoL = fub::euler::Density(equation_, left);
  const double rhoR = fub::euler::Density(equation_, right);
  const double rhoUL = fub::euler::Momentum(equation_, left, d);
  const double rhoUR = fub::euler::Momentum(equation_, right, d);
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  const double aL = fub::euler::SpeedOfSound(equation_, left);
  const double aR = fub::euler::SpeedOfSound(equation_, right);
  const double rhoEL = fub::euler::Energy(equation_, left);
  const double rhoER = fub::euler::Energy(equation_, right);
  const double pL = fub::euler::Pressure(equation_, left);
  const double pR = fub::euler::Pressure(equation_, right);
  const double hL = (rhoEL + pL) / rhoL;
  const double hR = (rhoER + pR) / rhoR;
  const double sqRhoL = std::sqrt(rhoL);
  const double sqRhoR = std::sqrt(rhoR);
  const double sqRhoSum = sqRhoL + sqRhoR;
  Array<double, 1, Dim> vL = fub::euler::Velocity(equation_, left);
  Array<double, 1, Dim> vR = fub::euler::Velocity(equation_, right);
  FUB_ASSERT(sqRhoSum > 0.0);
  Array<double, 1, Dim> roeU = (sqRhoL * vL + sqRhoR * vR) / sqRhoSum;
  const double roeU0 = roeU[d];
  const double roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;

  const double gammaL = fub::euler::Gamma(equation_, left);
  const double gammaR = fub::euler::Gamma(equation_, right);
  const double roeGamma = (sqRhoL * gammaL + sqRhoR * gammaR) / sqRhoSum;
  const double gm1 = roeGamma - 1.0;
  const double beta = std::sqrt(gm1 / (2 * roeGamma));

  const double roeA2 = gm1 * (roeH - 0.5 * roeU.matrix().squaredNorm());
  const double roeA = std::sqrt(roeA2);
  const double sL1 = uL - beta * aL;
  const double sL2 = roeU0 - roeA;
  const double sR1 = roeU0 + roeA;
  const double sR2 = uR + beta * aR;

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
      flux_hlle_, fluxL_, fluxR_, AsCons(left), AsCons(right));

  const double squaredNormRoeU = roeU.matrix().squaredNorm();
  const double squaredNormRoeU_half = 0.5 * squaredNormRoeU;
  const double u_bar = 0.5 * (sR + sL);
  const double u_bar_abs = std::abs(u_bar);
  const double delta = roeA / (roeA + u_bar_abs);
  const double b = bLbR / db_positive;

  const double deltaRho = rhoR - rhoL;
  const Array<double, 1, Dim> deltaRhoU = right.momentum - left.momentum;
  const double deltaRhoE = rhoER - rhoEL;

  if constexpr (Dim == 1) {
    const double l21 = gm1 / roeA2 * (roeH - roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * roeU[0];
    const double l23 = gm1 / roeA2 * (-1.0);
    const double alpha_2 =
        l21 * deltaRho + l22 * deltaRhoU[0] + l23 * deltaRhoE;
    const double b_delta_alpha_2 = b * delta * alpha_2;
    flux.density = flux_hlle_.density - b_delta_alpha_2 * 1.0;
    flux.momentum = flux_hlle_.momentum - b_delta_alpha_2 * roeU[0];
    flux.energy = flux_hlle_.energy - b_delta_alpha_2 * squaredNormRoeU_half;
    if constexpr (fub::euler::state_with_species<EulerEquation,
                                                 Conservative>()) {
      const int n_species = flux.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left, i);
        const double rhoYR = fub::euler::Species(equation_, right, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        // const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
        // const double deltaRhoY = rhoYR - rhoYL;
        // const double li1 = -roeY;
        // const double li4 = 1.0;
        // const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
        // const double b_delta_alpha_i = b * delta * alpha_i;
        const double f_rho = flux.density;
        const int upwind = (f_rho > 0.0);
        const double f_hllem = f_rho * (upwind * YL - (1 - upwind) * YR);
        flux.species[i] = f_hllem;
      }
    }
  } else if constexpr (Dim == 2) {
    const int ix = int(dir);
    const int iy = (ix + 1) % 2;
    const double l21 = gm1 / roeA2 * (roeH - roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * roeU[0];
    const double l23 = gm1 / roeA2 * roeU[1];
    const double l24 = gm1 / roeA2 * (-1.0);
    const double alpha_2 = l21 * deltaRho + l22 * deltaRhoU[0] +
                           l23 * deltaRhoU[1] + l24 * deltaRhoE;
    const double l31 = -roeU[iy];
    // const double l32 = 0;
    const double l33 = 1;
    // const double l34 = 0;
    const double alpha_3 = l31 * deltaRho + l33 * deltaRhoU[iy];
    const double b_delta = b * delta;
    const double b_delta_alpha_2 = b_delta * alpha_2;
    const double b_delta_alpha_3 = b_delta * alpha_3;
    flux.density = flux_hlle_.density - b_delta_alpha_2 * 1.0;
    flux.momentum[ix] = flux_hlle_.momentum[ix] - b_delta_alpha_2 * roeU[ix];
    flux.momentum[iy] = flux_hlle_.momentum[iy] - b_delta_alpha_2 * roeU[iy] -
                        b_delta_alpha_3 * 1.0;
    flux.energy = flux_hlle_.energy - b_delta_alpha_2 * squaredNormRoeU_half -
                  b_delta_alpha_3 * roeU[iy];
    if constexpr (fub::euler::state_with_species<EulerEquation,
                                                 Conservative>()) {
      const int n_species = flux.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left, i);
        const double rhoYR = fub::euler::Species(equation_, right, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
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
  } else {
    static_assert(Dim == 3);
    const int ix = int(dir);
    const int iy = (ix + 1) % 3;
    const int iz = (iy + 1) % 3;
    const double l21 = gm1 / roeA2 * (roeH - roeU.matrix().squaredNorm());
    const double l22 = gm1 / roeA2 * roeU[0];
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
    const double b_delta = b * delta;
    flux.density = flux_hlle_.density - b_delta * alpha_2;
    flux.momentum[ix] = flux_hlle_.momentum[ix] - b_delta * alpha_2 * roeU[ix];
    flux.momentum[iy] = flux_hlle_.momentum[iy] -
                        b_delta * (alpha_2 * roeU[iy] + alpha_3 * 1.0);
    flux.momentum[iz] = flux_hlle_.momentum[iz] -
                        b_delta * (alpha_2 * roeU[iz] + alpha_4 * 1.0);
    flux.energy =
        flux_hlle_.energy - b_delta * (alpha_2 * squaredNormRoeU_half +
                                       alpha_3 * roeU[iy] + alpha_4 * roeU[iz]);
    if constexpr (fub::euler::state_with_species<EulerEquation,
                                                 Conservative>()) {
      const int n_species = flux.species.size();
      for (int i = 0; i < n_species; ++i) {
        const double rhoYL = fub::euler::Species(equation_, left, i);
        const double rhoYR = fub::euler::Species(equation_, right, i);
        const double YL = rhoYL / rhoL;
        const double YR = rhoYR / rhoR;
        const double roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
        const double deltaRhoY = rhoYR - rhoYL;
        const double li1 = -roeY;
        const double li4 = 1.0;
        const double alpha_i = li1 * deltaRho + li4 * deltaRhoY;
        const double b_delta_alpha_i = b * delta * alpha_i;
        flux.species[i] =
            flux_hlle_.species[i] - b_delta * alpha_2 * roeY - b_delta_alpha_i;
      }
    }
  }
}

template <typename EulerEquation>
void Hllem<EulerEquation>::ComputeNumericFlux(
    ConservativeArray& flux, span<const CompleteArray, 2> states,
    Duration /* dt */, double /* dx */, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  ConservativeArray fluxL;
  ConservativeArray fluxR;
  fub::Flux(equation_, fluxL_array_, left, dir);
  fub::Flux(equation_, fluxR_array_, right, dir);

  // Compute Einfeldt signals velocitie

  const Array1d rhoL = euler::Density(equation_, left);
  const Array1d rhoR = euler::Density(equation_, right);
  const Array1d aL = euler::SpeedOfSound(equation_, left);
  const Array1d aR = euler::SpeedOfSound(equation_, right);
  const Array1d rhoEL = euler::Energy(equation_, left);
  const Array1d rhoER = euler::Energy(equation_, right);
  const Array1d pL = euler::Pressure(equation_, left);
  const Array1d pR = euler::Pressure(equation_, right);
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
  const Array1d squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  const Array1d gammaL = fub::euler::Gamma(equation_, left);
  const Array1d gammaR = fub::euler::Gamma(equation_, right);
  const Array1d roeGamma = (sqRhoL * gammaL + sqRhoR * gammaR) / sqRhoSum;
  const Array1d gm1 = roeGamma - 1.0;
  const Array1d beta = (gm1 / (2 * roeGamma)).sqrt();

  const Array1d roeA2 = gm1 * (roeH - squaredNormRoeU_half);
  const Array1d roeA = roeA2.sqrt();

  const Array1d sL1 = uL.row(int(dir)) - beta * aL;
  const Array1d sL2 = roeU.row(int(dir)) - roeA;
  const Array1d sR1 = roeU.row(int(dir)) + roeA;
  const Array1d sR2 = uR.row(int(dir)) + beta * aR;

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
      flux_hlle_array_, fluxL_array_, fluxR_array_, AsCons(left),
      AsCons(right));

  const Array1d u_bar = 0.5 * (sR + sL);
  const Array1d u_bar_abs = u_bar.abs();
  const Array1d delta = roeA / (roeA + u_bar_abs);
  const Array1d b = bLbR / db_positive;

  const Array1d deltaRho = rhoR - rhoL;
  const Array<double, Dim> deltaRhoU = right.momentum - left.momentum;
  const Array1d deltaRhoE = rhoER - rhoEL;

  const Array1d gm1_over_roeA2 = gm1 / roeA2;

  if constexpr (Dim == 1) {
    const Array1d l21 = gm1_over_roeA2 * (roeH - squaredNormRoeU);
    const Array1d l22 = gm1_over_roeA2 * roeU;
    const Array1d l23 = gm1_over_roeA2 * (-1.0);
    const Array1d alpha_2 = l21 * deltaRho + l22 * deltaRhoU + l23 * deltaRhoE;
    const Array1d b_delta_alpha_2 = b * delta * alpha_2;

    flux.density = flux_hlle_array_.density - b_delta_alpha_2;
    flux.momentum = flux_hlle_array_.momentum - b_delta_alpha_2 * roeU;
    flux.energy =
        flux_hlle_array_.energy - b_delta_alpha_2 * squaredNormRoeU_half;
    if constexpr (euler::state_with_species<EulerEquation, Conservative>()) {
      const int n_species = flux.species.rows();
      for (int i = 0; i < n_species; ++i) {
        const Array1d rhoYL = fub::euler::Species(equation_, left, i);
        const Array1d rhoYR = fub::euler::Species(equation_, right, i);
        const Array1d YL = rhoYL / rhoL;
        const Array1d YR = rhoYR / rhoR;
        const Array1d roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
        const Array1d deltaRhoY = rhoYR - rhoYL;
        const Array1d li1 = -roeY;
        const Array1d alpha_i = li1 * deltaRho + deltaRhoY;
        const Array1d b_delta_alpha_i = b * delta * alpha_i;
        flux.species.row(i) = flux_hlle_array_.species.row(i) -
                              b_delta_alpha_2 * roeY - b_delta_alpha_i;
      }
    }
  } else if constexpr (Dim == 2) {
    const int ix = int(dir);
    const int iy = (ix == 0);
    const Array1d l21 = gm1_over_roeA2 * (roeH - squaredNormRoeU);
    const Array1d l22 = gm1_over_roeA2 * roeU.row(0);
    const Array1d l23 = gm1_over_roeA2 * roeU.row(1);
    const Array1d l24 = -gm1_over_roeA2;
    const Array1d alpha_2 = l21 * deltaRho + l22 * deltaRhoU.row(0) +
                            l23 * deltaRhoU.row(1) + l24 * deltaRhoE;
    const Array1d l31 = -roeU.row(iy);
    // const Array1d l32 = 0;
    // const Array1d l33 = 1;
    // const Array1d l34 = 0;
    const Array1d alpha_3 = l31 * deltaRho + deltaRhoU.row(iy);
    const Array1d b_delta = b * delta;
    const Array1d b_delta_alpha_2 = b_delta * alpha_2;
    const Array1d b_delta_alpha_3 = b_delta * alpha_3;
    flux.density = flux_hlle_array_.density - b_delta_alpha_2;
    flux.momentum.row(ix) =
        flux_hlle_array_.momentum.row(ix) - b_delta_alpha_2 * roeU.row(ix);
    flux.momentum.row(iy) = flux_hlle_array_.momentum.row(iy) -
                            b_delta_alpha_2 * roeU.row(iy) - b_delta_alpha_3;
    flux.energy = flux_hlle_array_.energy -
                  b_delta_alpha_2 * squaredNormRoeU_half -
                  b_delta_alpha_3 * roeU.row(iy);
    if constexpr (euler::state_with_species<EulerEquation, Conservative>()) {
      const int n_species = flux.species.rows();
      for (int i = 0; i < n_species; ++i) {
        const Array1d rhoYL = fub::euler::Species(equation_, left, i);
        const Array1d rhoYR = fub::euler::Species(equation_, right, i);
        const Array1d YL = rhoYL / rhoL;
        const Array1d YR = rhoYR / rhoR;
        const Array1d roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
        const Array1d deltaRhoY = rhoYR - rhoYL;
        const Array1d li1 = -roeY;
        const Array1d alpha_i = li1 * deltaRho + deltaRhoY;
        const Array1d b_delta_alpha_i = b_delta * alpha_i;
        flux.species.row(i) = flux_hlle_array_.species.row(i) -
                              b_delta_alpha_2 * roeY - b_delta_alpha_i;
      }
    }
  } else {
    static_assert(Dim == 3);
    const int ix = int(dir);
    const int iy = (ix + 1) % 3;
    const int iz = (iy + 1) % 3;
    const Array1d l21 = gm1_over_roeA2 * (roeH - squaredNormRoeU);
    const Array1d l22 = gm1_over_roeA2 * roeU.row(0);
    const Array1d l23 = gm1_over_roeA2 * roeU.row(1);
    const Array1d l24 = gm1_over_roeA2 * roeU.row(2);
    const Array1d l25 = -gm1_over_roeA2;
    const Array1d l31 = -roeU.row(iy);
    // const Array1d l32 = 0;
    // const Array1d l33 = 1.0;
    // const Array1d l34 = 0;
    const Array1d l41 = -roeU.row(iz);
    // const Array1d l32 = 0;
    // const Array1d l33 = 0;
    // const Array1d l44 = 1.0;
    const Array1d alpha_2 = l21 * deltaRho + l22 * deltaRhoU.row(0) +
                            l23 * deltaRhoU.row(1) + l24 * deltaRhoU.row(2) +
                            l25 * deltaRhoE;
    const Array1d alpha_3 = l31 * deltaRho + deltaRhoU.row(iy);
    const Array1d alpha_4 = l41 * deltaRho + deltaRhoU.row(iz);
    const Array1d b_delta = b * delta;
    flux.density = flux_hlle_array_.density - b_delta * alpha_2;
    flux.momentum.row(ix) =
        flux_hlle_array_.momentum.row(ix) - b_delta * alpha_2 * roeU.row(iy);
    flux.momentum.row(iy) = flux_hlle_array_.momentum.row(iy) -
                            b_delta * (alpha_2 * roeU.row(iy) + alpha_3);
    flux.momentum.row(iz) = flux_hlle_array_.momentum.row(iz) -
                            b_delta * (alpha_2 * roeU.row(iz) + alpha_4);
    flux.energy = flux_hlle_array_.energy -
                  b_delta * (alpha_2 * squaredNormRoeU_half +
                             alpha_3 * roeU.row(iy) + alpha_4 * roeU.row(iz));
    if constexpr (euler::state_with_species<EulerEquation, Conservative>()) {
      const int n_species = flux.species.rows();
      for (int i = 0; i < n_species; ++i) {
        const Array1d rhoYL = fub::euler::Species(equation_, left, i);
        const Array1d rhoYR = fub::euler::Species(equation_, right, i);
        const Array1d YL = rhoYL / rhoL;
        const Array1d YR = rhoYR / rhoR;
        const Array1d roeY = (sqRhoL * YL + sqRhoR * YR) / sqRhoSum;
        const Array1d deltaRhoY = rhoYR - rhoYL;
        const Array1d li1 = -roeY;
        const Array1d alpha_i = li1 * deltaRho + deltaRhoY;
        const Array1d b_delta_alpha_i = b_delta * alpha_i;
        flux.species.row(i) = flux_hlle_array_.species.row(i) -
                              b_delta * alpha_2 * roeY - b_delta_alpha_i;
      }
    }
  }
}

template <typename EulerEquation>
void Hllem<EulerEquation>::ComputeNumericFlux(
    ConservativeArray& flux, Array1d face_fractions,
    span<const CompleteArray, 2> states,
    span<const Array1d, 2> /* volume_fractions */, Duration /* dt */,
    double /* dx */, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  MaskArray face_mask = face_fractions > 0.0;

  ConservativeArray fluxL;
  ConservativeArray fluxR;
  equation_.Flux(fluxL, left, face_mask, dir);
  equation_.Flux(fluxR, right, face_mask, dir);

  const Array1d gm1 = (equation_.gamma_array_ - Array1d::Constant(1.0));
  const Array1d beta = (gm1 / (2 * equation_.gamma_array_)).sqrt();

  // Compute Einfeldt signals velocities

  const Array1d ones = Array1d::Constant(1.0);
  const Array1d zeros = Array1d::Constant(0.0);

  const Array1d rhoL = face_mask.select(left.density, ones);
  const Array1d rhoR = face_mask.select(right.density, ones);
  FUB_ASSERT((rhoL > 0.0).all());
  FUB_ASSERT((rhoR > 0.0).all());
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
    uL.row(i) = face_mask.select(left.momentum.row(i) / rhoL, zeros);
    uR.row(i) = face_mask.select(right.momentum.row(i) / rhoR, zeros);
    roeU.row(i) = sqRhoL_over_Sum * uL.row(i) + sqRhoR_over_Sum * uR.row(i);
  }
  const Array1d roeH = sqRhoL_over_Sum * hL + sqRhoR_over_Sum * hR;

  Array1d squaredNormRoeU = Array1d::Zero();
  for (int i = 0; i < Dim; ++i) {
    squaredNormRoeU += roeU.row(i) * roeU.row(i);
  }
  const Array1d squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  const Array1d roeA2 = gm1 * (roeH - squaredNormRoeU_half);
  const Array1d roeA = roeA2.sqrt();

  const Array1d sL1 = uL.row(int(dir)) - beta * aL;
  const Array1d sL2 = roeU.row(int(dir)) - roeA;
  const Array1d sR1 = roeU.row(int(dir)) + roeA;
  const Array1d sR2 = uR.row(int(dir)) + beta * aR;

  const Array1d sL = sL1.min(sL2);
  const Array1d sR = sR1.max(sR2);

  const Array1d bL = sL.min(zeros);
  const Array1d bR = sR.max(zeros);
  const Array1d bLbR = bL * bR;
  const Array1d db = bR - bL;
  const Array1d db_positive = (db > 0.0).select(db, ones);

  ConservativeArray flux_hlle{};
  ForEachComponent(
      [&](auto&& nf, Array1d fL, Array1d fR, Array1d qL, Array1d qR) {
        nf = (bR * fL - bL * fR + bLbR * (qR - qL)) / db_positive;
      },
      flux_hlle, fluxL, fluxR, AsCons(left), AsCons(right));

  const Array1d u_bar = 0.5 * (sR + sL);
  const Array1d u_bar_abs = u_bar.abs();

  const Array1d u_signal = face_mask.select(roeA + u_bar_abs, ones);
  const Array1d delta = roeA / u_signal;
  const Array1d b = bLbR / db_positive;

  const Array1d deltaRho = rhoR - rhoL;
  const Array<double, Dim> deltaRhoU = right.momentum - left.momentum;
  const Array1d deltaRhoE = rhoER - rhoEL;

  const Array1d gm1_over_roeA2 = face_mask.select(gm1 / roeA2, zeros);

  if constexpr (Dim == 1) {
    const Array1d l21 = gm1_over_roeA2 * (roeH - squaredNormRoeU);
    const Array1d l22 = gm1_over_roeA2 * roeU;
    const Array1d l23 = gm1_over_roeA2 * (-1.0);
    const Array1d alpha_2 = l21 * deltaRho + l22 * deltaRhoU + l23 * deltaRhoE;
    const Array1d b_delta_alpha_2 = b * delta * alpha_2;

    flux.density = flux_hlle.density - b_delta_alpha_2;
    flux.momentum = flux_hlle.momentum - b_delta_alpha_2 * roeU;
    flux.energy = flux_hlle.energy - b_delta_alpha_2 * squaredNormRoeU_half;
  } else if constexpr (Dim == 2) {
    const int ix = int(dir);
    const int iy = (ix == 0);
    const Array1d l21 = gm1_over_roeA2 * (roeH - squaredNormRoeU);
    const Array1d l22 = gm1_over_roeA2 * roeU.row(0);
    const Array1d l23 = gm1_over_roeA2 * roeU.row(1);
    const Array1d l24 = -gm1_over_roeA2;
    const Array1d alpha_2 = l21 * deltaRho + l22 * deltaRhoU.row(0) +
                            l23 * deltaRhoU.row(1) + l24 * deltaRhoE;
    const Array1d l31 = -roeU.row(iy);
    // const Array1d l32 = 0;
    // const Array1d l33 = 1;
    // const Array1d l34 = 0;
    const Array1d alpha_3 = l31 * deltaRho + deltaRhoU.row(iy);
    const Array1d b_delta = b * delta;
    const Array1d b_delta_alpha_2 = b_delta * alpha_2;
    const Array1d b_delta_alpha_3 = b_delta * alpha_3;
    flux.density = flux_hlle.density - b_delta_alpha_2;
    flux.momentum.row(ix) =
        flux_hlle.momentum.row(ix) - b_delta_alpha_2 * roeU.row(ix);
    flux.momentum.row(iy) = flux_hlle.momentum.row(iy) -
                            b_delta_alpha_2 * roeU.row(iy) - b_delta_alpha_3;
    flux.energy = flux_hlle.energy - b_delta_alpha_2 * squaredNormRoeU_half -
                  b_delta_alpha_3 * roeU.row(iy);
  } else {
    static_assert(Dim == 3);
    const int ix = int(dir);
    const int iy = (ix + 1) % 3;
    const int iz = (iy + 1) % 3;
    const Array1d l21 = gm1_over_roeA2 * (roeH - squaredNormRoeU);
    const Array1d l22 = gm1_over_roeA2 * roeU.row(0);
    const Array1d l23 = gm1_over_roeA2 * roeU.row(1);
    const Array1d l24 = gm1_over_roeA2 * roeU.row(2);
    const Array1d l25 = -gm1_over_roeA2;
    const Array1d l31 = -roeU.row(iy);
    // const Array1d l32 = 0;
    // const Array1d l33 = 1.0;
    // const Array1d l34 = 0;
    const Array1d l41 = -roeU.row(iz);
    // const Array1d l32 = 0;
    // const Array1d l33 = 0;
    // const Array1d l44 = 1.0;
    const Array1d alpha_2 = l21 * deltaRho + l22 * deltaRhoU.row(0) +
                            l23 * deltaRhoU.row(1) + l24 * deltaRhoU.row(2) +
                            l25 * deltaRhoE;
    const Array1d alpha_3 = l31 * deltaRho + deltaRhoU.row(iy);
    const Array1d alpha_4 = l41 * deltaRho + deltaRhoU.row(iz);
    const Array1d b_delta = b * delta;
    flux.density = flux_hlle.density - b_delta * alpha_2;
    flux.momentum.row(ix) =
        flux_hlle.momentum.row(ix) - b_delta * alpha_2 * roeU.row(iy);
    flux.momentum.row(iy) = flux_hlle.momentum.row(iy) -
                            b_delta * (alpha_2 * roeU.row(iy) + alpha_3);
    flux.momentum.row(iz) = flux_hlle.momentum.row(iz) -
                            b_delta * (alpha_2 * roeU.row(iz) + alpha_4);
    flux.energy = flux_hlle.energy -
                  b_delta * (alpha_2 * squaredNormRoeU_half +
                             alpha_3 * roeU.row(iy) + alpha_4 * roeU.row(iz));
  }
  FUB_ASSERT(!flux.density.isNaN().any());
  FUB_ASSERT(!flux.momentum.isNaN().any());
  FUB_ASSERT(!flux.energy.isNaN().any());
}

template <typename EulerEquation>
double Hllem<EulerEquation>::ComputeStableDt(span<const Complete, 2> states,
                                             double dx, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  const Complete& left = states[0];
  const Complete& right = states[1];

  const double gm1 = equation_.gamma - 1.0;
  const double beta = gm1 / (2 * equation_.gamma);

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
  const double roeA2 = gm1 * (roeH - 0.5 * roeU.matrix().squaredNorm());
  const double roeA = std::sqrt(roeA2);
  const double sL1 = uL - beta * aL;
  const double sL2 = roeU0 - roeA;
  const double sR1 = roeU0 + roeA;
  const double sR2 = uR + beta * aR;

  const double sL = std::min(sL1, sL2);
  const double sR = std::max(sR1, sR2);

  const double sL_abs = std::abs(sL);
  const double sR_abs = std::abs(sR);

  const double maxS = std::max(sL_abs, sR_abs);
  const double max_dt = dx / maxS;

  return max_dt;
}

template <typename EulerEquation>
Array1d
Hllem<EulerEquation>::ComputeStableDt(span<const CompleteArray, 2> states,
                                      double dx, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  const Array1d gm1 = (equation_.gamma_array_ - Array1d::Constant(1.0));
  const Array1d beta = gm1 / (2 * equation_.gamma_array_);

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
  const Array1d squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  const Array1d roeA2 = gm1 * (roeH - squaredNormRoeU_half);
  const Array1d roeA = roeA2.sqrt();

  const Array1d sL1 = uL.row(int(dir)) - beta * aL;
  const Array1d sL2 = roeU.row(int(dir)) - roeA;
  const Array1d sR1 = roeU.row(int(dir)) + roeA;
  const Array1d sR2 = uR.row(int(dir)) + beta * aR;

  const Array1d sL = sL1.min(sL2);
  const Array1d sR = sR1.max(sR2);

  const Array1d sL_abs = sL.abs();
  const Array1d sR_abs = sR.abs();

  const Array1d maxS = sL_abs.max(sR_abs);
  const Array1d max_dt = Array1d::Constant(dx) / maxS;

  return max_dt;
}

template <typename EulerEquation>
Array1d Hllem<EulerEquation>::ComputeStableDt(
    span<const CompleteArray, 2> states, Array1d face_fraction,
    span<const Array1d, 2>, double dx, Direction dir) {
  static constexpr int Dim = EulerEquation::Rank();
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  MaskArray face_mask = face_fraction > 0.0;

  const Array1d gm1 = (equation_.gamma_array_ - Array1d::Constant(1.0));
  const Array1d beta = gm1 / (2 * equation_.gamma_array_);

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
  const Array1d squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  const Array1d roeA2 = gm1 * (roeH - squaredNormRoeU_half);
  const Array1d roeA = roeA2.sqrt();

  const Array1d sL1 = uL.row(int(dir)) - beta * aL;
  const Array1d sL2 = roeU.row(int(dir)) - roeA;
  const Array1d sR1 = roeU.row(int(dir)) + roeA;
  const Array1d sR2 = uR.row(int(dir)) + beta * aR;

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