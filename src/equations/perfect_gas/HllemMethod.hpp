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

template <int Dim>
void Hllem<Dim>::ComputeNumericFlux(
    Conservative& flux, span<const Complete, 2> states,
    Duration /* dt */, double /* dx */, Direction dir) {
  const Complete& left = states[0];
  const Complete& right = states[1];

  Conservative fluxL;
  Conservative fluxR;
  equation_.Flux(fluxL, left, dir);
  equation_.Flux(fluxR, right, dir);

  const double gm1 = equation_.gamma - 1.0;

  // Compute Einfeldt signals velocities
  
  int d = static_cast<int>(dir);

  const double rhoL = left.density;
  const double rhoR = right.density;
  const double rhoUL = left.momentum[d];
  const double rhoUR = right.momentum[d];
  const double aL = left.speed_of_sound;
  const double aR = right.speed_of_sound;
  const double hL = (left.energy + left.pressure) / rhoL;
  const double hR = (right.energy + right.pressure) / rhoR;
  const double sqRhoL = std::sqrt(rhoL);
  const double sqRhoR = std::sqrt(rhoR);
  const double sqRhoSum = sqRhoL + sqRhoR;
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  const double roeU0 = (sqRhoL * uL + sqRhoR * uR) / sqRhoSum;
  const double roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;
  const double roeA2 = gm1 * (roeH  - 0.5 * roeU0 * roeU0);
  const double roeA = std::sqrt(roeA2);
  const double sL1 = uL - aL;
  const double sL2 = roeU0 - 0.5 * roeA;
  const double sR1 = roeU0 + 0.5 * roeA;
  const double sR2 = uR + aR;
  
  const double sL = std::min(sL1, sL2);
  const double sR = std::max(sR1, sR2);

  const double bL = std::min(sL, 0.0);
  const double bR = std::max(sR, 0.0);
  const double bLbR = bL * bR;
  const double db = bR - bL;
  const double db_positive = int(db == 0) + int(db != 0) * db;

  Conservative flux_hlle;
  ForEachComponent(
      [&](double& nf, double fL, double fR, double qL, double qR) {
        nf = (bR * fL - bL * fR + bLbR * (qR - qL)) / db_positive;
      },
      flux_hlle, fluxL, fluxR, AsCons(left), AsCons(right));

//  const double roeRho = sqRhoL * sqRhoR;
  std::array<double, Dim> vL;
  std::array<double, Dim> vR;
  std::array<double, Dim> roeU;
  for (int i = 0; i < Dim; ++i) {
    vL[i] = left.momentum[i] / rhoL;
    vR[i] = right.momentum[i] / rhoR;
    roeU[i] = (sqRhoL * vL[i] + sqRhoR * vR[i]) / sqRhoSum;
  }
  double squaredNormRoeU = 0.0;
  for (int i = 0; i < Dim; ++i) {
    squaredNormRoeU += roeU[i] * roeU[i];
  }
  const double u_bar = 0.5 * (sL + sR);
  const double u_bar_abs = std::abs(u_bar);
  const double delta = roeA / (roeA + u_bar_abs);
  
  const double b = bLbR / db_positive;
  const double b_delta = b * delta * roeU0;
  const double squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  if constexpr (Dim == 1) {
    flux.density  = flux_hlle.density  - b_delta * 1.0;
    flux.momentum = flux_hlle.momentum - b_delta * roeU[0];
    flux.energy   = flux_hlle.energy   - b_delta * squaredNormRoeU_half;
  } else if constexpr (Dim == 2) {
    const int i = int(dir);
    const int j = int(i == 0);
    flux.density         = flux_hlle.density     - b_delta;
    flux.momentum[i] = flux_hlle.momentum[i] - b_delta *  roeU[i];
    flux.momentum[j] = flux_hlle.momentum[j] - b_delta * (roeU[j]              + 1.0);
    flux.energy          = flux_hlle.energy      - b_delta * (squaredNormRoeU_half + roeU[j]);
  } else {
    static_assert(Dim == 3);
    const int i = int(dir);
    const int j = (i + 1) % 3;
    const int k = (j + 1) % 3;
    flux.density     = flux_hlle.density     - b_delta;
    flux.momentum[i] = flux_hlle.momentum[i] - b_delta *  roeU[i];
    flux.momentum[j] = flux_hlle.momentum[j] - b_delta * (roeU[j]              + 1.0);
    flux.momentum[k] = flux_hlle.momentum[k] - b_delta * (roeU[k]                        + 1.0);
    flux.energy      = flux_hlle.energy      - b_delta * (squaredNormRoeU_half + roeU[j] + roeU[k]);
  }
}

template <int Dim>
void Hllem<Dim>::ComputeNumericFlux(
    ConservativeArray& flux, span<const CompleteArray, 2> states,
    Duration /* dt */, double /* dx */, Direction dir) {
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  ConservativeArray fluxL;
  ConservativeArray fluxR;
  equation_.Flux(fluxL, left, dir);
  equation_.Flux(fluxR, right, dir);

  const Array1d gm1 = (equation_.gamma_array_ - Array1d::Constant(1.0));

  // Compute Einfeldt signals velocities
    
  const Array1d rhoL = left.density;
  const Array1d rhoR = right.density;
  const Array1d rhoUL = left.momentum.row(int(dir));
  const Array1d rhoUR = right.momentum.row(int(dir));
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
  const Array1d uL = rhoUL / rhoL;
  const Array1d uR = rhoUR / rhoR;
  const Array1d roeU0 = (sqRhoL * uL + sqRhoR * uR) / sqRhoSum;
  const Array1d roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;
  const Array1d roeA2 = gm1 * (roeH  - 0.5 * roeU0 * roeU0);
  const Array1d roeA = roeA2.sqrt();
  const Array1d sL1 = uL - aL;
  const Array1d sL2 = roeU0 - 0.5 * roeA;
  const Array1d sR1 = roeU0 + 0.5 * roeA;
  const Array1d sR2 = uR + aR;
  
  const Array1d sL = sL1.min(sL2);
  const Array1d sR = sR1.max(sR2);

  const Array1d zero = Array1d::Zero();
  const Array1d bL = sL.min(zero);
  const Array1d bR = sR.max(zero);
  const Array1d bLbR = bL * bR;
  const Array1d db = bR - bL;
  const Array1d db_positive = (db > 0.0).select(db, Array1d::Constant(1.0));

  ConservativeArray flux_hlle{};
  ForEachComponent(
      [&](auto&& nf, Array1d fL, Array1d fR, Array1d qL, Array1d qR) {
        nf = (bR * fL - bL * fR + bLbR * (qR - qL)) / db_positive;
      },
      flux_hlle, fluxL, fluxR, AsCons(left), AsCons(right));

  Array<double, Dim> vL;
  Array<double, Dim> vR;
  Array<double, Dim> roeU;
  for (int i = 0; i < Dim; ++i) {
    vL.row(i) = left.momentum.row(i) / rhoL;
    vR.row(i) = right.momentum.row(i) / rhoR;
    roeU.row(i) = (sqRhoL * vL.row(i) + sqRhoR * vR.row(i)) / sqRhoSum;
  }
  Array1d squaredNormRoeU = Array1d::Zero();
  for (int i = 0; i < Dim; ++i) {
    squaredNormRoeU += roeU.row(i) * roeU.row(i);
  }
  const Array1d u_bar = 0.5 * (sR + sL);
  const Array1d c_bar = 0.5 * (sR - sL);
  const Array1d u_bar_abs = u_bar.abs();
  const Array1d delta = c_bar / (c_bar + u_bar_abs);
  const Array1d b = bLbR / db_positive;

  const Array1d deltaRho = rhoR - rhoL;
  const Array1d deltaRhoU = rhoUR - rhoUL;
  const Array1d deltaRhoE = rhoER - rhoEL;
  [[maybe_unused]] const Array<double, Dim> deltaRhoV = right.momentum - left.momentum;

  const Array1d C1 = (deltaRhoU - deltaRho * roeU0) / roeA;

  const Array1d squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  if constexpr (Dim == 1) {
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - C1 * roeU * roeA) / (roeH - squaredNormRoeU_half);
    const Array1d alpha_2 = deltaRho - C2;
    const Array1d b_delta_alpha = b * delta * alpha_2;

    flux.density  = flux_hlle.density  - b_delta_alpha * 1.0;
    flux.momentum = flux_hlle.momentum - b_delta_alpha * roeU;
    flux.energy   = flux_hlle.energy   - b_delta_alpha * squaredNormRoeU_half;
  } else if constexpr (Dim == 2) {
    const int i = int(dir);
    const int j = int(i == 0);
    const Array1d alpha_3 = deltaRhoV.row(j) - deltaRho * roeU.row(j);
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - alpha_3 * roeU.row(j) - C1 * roeU0 * roeA) / (roeH - squaredNormRoeU_half);
    const Array1d alpha_2 = deltaRho - C2;
    const Array1d b_delta = b * delta;
    flux.density         = flux_hlle.density         - b_delta * alpha_2 * 1.0;
    flux.momentum.row(i) = flux_hlle.momentum.row(i) - b_delta * alpha_2 * roeU.row(i);
    flux.momentum.row(j) = flux_hlle.momentum.row(j) - b_delta * (alpha_2 * roeU.row(j) + alpha_3 * 1.0);
    flux.energy          = flux_hlle.energy          - b_delta * (alpha_2 * squaredNormRoeU_half + alpha_3 * roeU.row(j));
  } else {
    static_assert(Dim == 3);
    const int i = int(dir);
    const int j = (i + 1) % 3;
    const int k = (j + 1) % 3;
    const Array1d alpha_3 = deltaRhoV.row(j) - deltaRho * roeU.row(j);
    const Array1d alpha_4 = deltaRhoV.row(k) - deltaRho * roeU.row(k);
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - alpha_3 * roeU.row(j) - alpha_4 * roeU.row(k) - C1 * roeU0 * roeA) / (roeH - squaredNormRoeU_half);
    const Array1d alpha_2 = deltaRho - C2;
    const Array1d b_delta = b * delta;
    flux.density         = flux_hlle.density         - b_delta * alpha_2;
    flux.momentum.row(i) = flux_hlle.momentum.row(i) - b_delta * alpha_2 * roeU.row(i);
    flux.momentum.row(j) = flux_hlle.momentum.row(j) - b_delta * (alpha_2 * roeU.row(j) + alpha_3 * 1.0);
    flux.momentum.row(k) = flux_hlle.momentum.row(k) - b_delta * (alpha_2 * roeU.row(k) + alpha_4 * 1.0);
    flux.energy          = flux_hlle.energy          - b_delta * (alpha_2 * squaredNormRoeU_half + alpha_3 * roeU.row(j) + alpha_4 * roeU.row(k));
  }
}

template <int Dim>
void Hllem<Dim>::ComputeNumericFlux(
    ConservativeArray& flux, Array1d face_fractions,
    span<const CompleteArray, 2> states, span<const Array1d, 2> /* volume_fractions */,
    Duration /* dt */, double /* dx */, Direction dir) {
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  MaskArray face_mask = face_fractions > 0.0;

  ConservativeArray fluxL;
  ConservativeArray fluxR;
  equation_.Flux(fluxL, left, face_mask, dir);
  equation_.Flux(fluxR, right, face_mask, dir);

  const Array1d gm1 = (equation_.gamma_array_ - Array1d::Constant(1.0));

  // Compute Einfeldt signals velocities
    
  const Array1d ones = Array1d::Constant(1.0);
  const Array1d zeros = Array1d::Constant(0.0);

  const Array1d rhoL = face_mask.select(left.density, ones);
  const Array1d rhoR = face_mask.select(right.density, ones);
  const Array1d rhoUL = face_mask.select(left.momentum.row(int(dir)), zeros);
  const Array1d rhoUR = face_mask.select(right.momentum.row(int(dir)), zeros);
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
  const Array1d uL = rhoUL / rhoL;
  const Array1d uR = rhoUR / rhoR;
  const Array1d roeU0 = (sqRhoL * uL + sqRhoR * uR) / sqRhoSum;
  const Array1d roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;
  const Array1d roeA2 = gm1 * (roeH  - 0.5 * roeU0 * roeU0);
  const Array1d roeA = roeA2.sqrt();
  const Array1d sL1 = uL - aL;
  const Array1d sL2 = roeU0 - 0.5 * roeA;
  const Array1d sR1 = roeU0 + 0.5 * roeA;
  const Array1d sR2 = uR + aR;
  
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

//  const Array1d roeRho = sqRhoL * sqRhoR;
  Array<double, Dim> vL;
  Array<double, Dim> vR;
  Array<double, Dim> roeU;
  for (int i = 0; i < Dim; ++i) {
    vL.row(i) = face_mask.select(left.momentum.row(i) / rhoL, zeros);
    vR.row(i) = face_mask.select(right.momentum.row(i) / rhoR, zeros);
    roeU.row(i) = (sqRhoL * vL.row(i) + sqRhoR * vR.row(i)) / sqRhoSum;
  }
  Array1d squaredNormRoeU = Array1d::Zero();
  for (int i = 0; i < Dim; ++i) {
    squaredNormRoeU += roeU.row(i) * roeU.row(i);
  }
  const Array1d u_bar = 0.5 * (sL + sR);
  const Array1d u_bar_abs = u_bar.abs();
  const Array1d u_signal = face_mask.select(roeA + u_bar_abs, ones);
  const Array1d delta = roeA / u_signal;

  const Array1d deltaRho = rhoR - rhoL;
  const Array<double, Dim> deltaRhoU = right.momentum - left.momentum;
  const Array1d deltaRhoE = rhoER - rhoEL;

  const Array1d roeA_positive = face_mask.select(roeA, ones);
  const Array1d C1 = (deltaRhoU.row(int(dir)) - deltaRho * roeU0) / roeA_positive;

  const Array1d squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  const Array1d b = bLbR / db_positive;
  const Array1d b_delta = b * delta;

  const Array1d roeH_minus_kinetic = face_mask.select(roeH - squaredNormRoeU_half, ones);

  if constexpr (Dim == 1) {
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - C1 * roeU * roeA) / roeH_minus_kinetic;
    const Array1d alpha_2 = deltaRho - C2;
    const Array1d b_delta_alpha = b_delta * alpha_2;

    flux.density  = flux_hlle.density  - b_delta_alpha * 1.0;
    flux.momentum = flux_hlle.momentum - b_delta_alpha * roeU;
    flux.energy   = flux_hlle.energy   - b_delta_alpha * squaredNormRoeU_half;
  } else if constexpr (Dim == 2) {
    const int i = int(dir);
    const int j = int(i == 0);
    const Array1d alpha_3 = deltaRhoU.row(j) - deltaRho * roeU.row(j);
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - alpha_3 * roeU.row(j) - C1 * roeU0 * roeA) / roeH_minus_kinetic;
    const Array1d alpha_2 = deltaRho - C2;
    flux.density         = flux_hlle.density         - b_delta * alpha_2 * 1.0;
    flux.momentum.row(i) = flux_hlle.momentum.row(i) - b_delta * alpha_2 * roeU.row(i);
    flux.momentum.row(j) = flux_hlle.momentum.row(j) - b_delta * (alpha_2 * roeU.row(j) + alpha_3 * 1.0);
    flux.energy          = flux_hlle.energy          - b_delta * (alpha_2 * squaredNormRoeU_half + alpha_3 * roeU.row(j));
  } else {
    static_assert(Dim == 3);
    const int i = int(dir);
    const int j = (i + 1) % 3;
    const int k = (j + 1) % 3;
    const Array1d alpha_3 = deltaRhoU.row(j) - deltaRho * roeU.row(j);
    const Array1d alpha_4 = deltaRhoU.row(k) - deltaRho * roeU.row(k);
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - alpha_3 * roeU.row(j) - alpha_4 * roeU.row(k) - C1 * roeU0 * roeA) / roeH_minus_kinetic;
    const Array1d alpha_2 = deltaRho - C2;
    flux.density         = flux_hlle.density         - b_delta * alpha_2;
    flux.momentum.row(i) = flux_hlle.momentum.row(i) - b_delta * alpha_2 * roeU.row(i);
    flux.momentum.row(j) = flux_hlle.momentum.row(j) - b_delta * (alpha_2 * roeU.row(j) + alpha_3 * 1.0);
    flux.momentum.row(k) = flux_hlle.momentum.row(k) - b_delta * (alpha_2 * roeU.row(k) + alpha_4 * 1.0);
    flux.energy          = flux_hlle.energy          - b_delta * (alpha_2 * squaredNormRoeU_half + alpha_3 * roeU.row(j) + alpha_4 * roeU.row(k));
  }
}

template <int Dim>
double Hllem<Dim>::ComputeStableDt(
    span<const Complete, 2> states, double dx, Direction dir)
{
  const Complete& left = states[0];
  const Complete& right = states[1];

  const double gm1 = equation_.gamma - 1.0;

  // Compute Einfeldt signals velocities
  
  int d = static_cast<int>(dir);

  const double rhoL = left.density;
  const double rhoR = right.density;
  const double rhoUL = left.momentum[d];
  const double rhoUR = right.momentum[d];
  const double aL = left.speed_of_sound;
  const double aR = right.speed_of_sound;
  const double hL = (left.energy + left.pressure) / rhoL;
  const double hR = (right.energy + right.pressure) / rhoR;
  const double sqRhoL = std::sqrt(rhoL);
  const double sqRhoR = std::sqrt(rhoR);
  const double sqRhoSum = sqRhoL + sqRhoR;
  const double uL = rhoUL / rhoL;
  const double uR = rhoUR / rhoR;
  const double roeU0 = (sqRhoL * uL + sqRhoR * uR) / sqRhoSum;
  const double roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;
  const double roeA2 = gm1 * (roeH  - 0.5 * roeU0 * roeU0);
  const double roeA = std::sqrt(roeA2);
  const double sL1 = uL - aL;
  const double sL2 = roeU0 - 0.5 * roeA;
  const double sR1 = roeU0 + 0.5 * roeA;
  const double sR2 = uR + aR;
  
  const double sL = std::min(sL1, sL2);
  const double sR = std::max(sR1, sR2);

  const double sL_abs = std::abs(sL);
  const double sR_abs = std::abs(sR);

  const double maxS = std::max(sL_abs, sR_abs);
  const double max_dt = dx / maxS;

  return max_dt;
}


template <int Dim>
Array1d Hllem<Dim>::ComputeStableDt(
    span<const CompleteArray, 2> states, double dx, Direction dir)
{
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  const Array1d gm1 = (equation_.gamma_array_ - Array1d::Constant(1.0));

  // Compute Einfeldt signals velocities
    
  const Array1d rhoL = left.density;
  const Array1d rhoR = right.density;
  const Array1d rhoUL = left.momentum.row(int(dir));
  const Array1d rhoUR = right.momentum.row(int(dir));
  const Array1d aL = left.speed_of_sound;
  const Array1d aR = right.speed_of_sound;
  const Array1d hL = (left.energy + left.pressure) / rhoL;
  const Array1d hR = (right.energy + right.pressure) / rhoR;
  const Array1d sqRhoL = rhoL.sqrt();
  const Array1d sqRhoR = rhoR.sqrt();
  const Array1d sqRhoSum = sqRhoL + sqRhoR;
  const Array1d uL = rhoUL / rhoL;
  const Array1d uR = rhoUR / rhoR;
  const Array1d roeU0 = (sqRhoL * uL + sqRhoR * uR) / sqRhoSum;
  const Array1d roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;
  const Array1d roeA2 = gm1 * (roeH  - 0.5 * roeU0 * roeU0);
  const Array1d roeA = roeA2.sqrt();
  const Array1d sL1 = uL - aL;
  const Array1d sL2 = roeU0 - 0.5 * roeA;
  const Array1d sR1 = roeU0 + 0.5 * roeA;
  const Array1d sR2 = uR + aR;
  
  const Array1d sL = sL1.min(sL2);
  const Array1d sR = sR1.max(sR2);

  const Array1d sL_abs = sL.abs();
  const Array1d sR_abs = sR.abs();

  const Array1d maxS = sL_abs.max(sR_abs);
  const Array1d max_dt = Array1d::Constant(dx) / maxS;

  return max_dt;
}

template <int Dim>
Array1d Hllem<Dim>::ComputeStableDt(
    span<const CompleteArray, 2> states, Array1d face_fraction,
    span<const Array1d, 2>, double dx, Direction dir)
{
  const CompleteArray& left = states[0];
  const CompleteArray& right = states[1];

  MaskArray mask = face_fraction > 0.0;

  const Array1d gm1 = (equation_.gamma_array_ - Array1d::Constant(1.0));

  // Compute Einfeldt signals velocities
    
  const Array1d ones = Array1d::Constant(1.0);
  const Array1d zeros = Array1d::Constant(0.0);

  const Array1d rhoL = mask.select(left.density, ones);
  const Array1d rhoR = mask.select(right.density, ones);
  const Array1d rhoUL = mask.select(left.momentum.row(int(dir)), zeros);
  const Array1d rhoUR = mask.select(right.momentum.row(int(dir)), zeros);
  const Array1d aL = mask.select(left.speed_of_sound, zeros);
  const Array1d aR = mask.select(right.speed_of_sound, zeros);
  const Array1d eL = mask.select(left.energy, zeros);
  const Array1d eR = mask.select(right.energy, zeros);
  const Array1d pL = mask.select(left.pressure, zeros);
  const Array1d pR = mask.select(right.pressure, zeros);
  const Array1d hL = (eL + pL) / rhoL;
  const Array1d hR = (eR + pR) / rhoR;
  const Array1d sqRhoL = rhoL.sqrt();
  const Array1d sqRhoR = rhoR.sqrt();
  const Array1d sqRhoSum = sqRhoL + sqRhoR;
  const Array1d uL = rhoUL / rhoL;
  const Array1d uR = rhoUR / rhoR;
  const Array1d roeU0 = (sqRhoL * uL + sqRhoR * uR) / sqRhoSum;
  const Array1d roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;
  const Array1d roeA2 = (gm1 * (roeH  - 0.5 * roeU0 * roeU0), zeros);
  const Array1d roeA = roeA2.sqrt();
  const Array1d sL1 = uL - aL;
  const Array1d sL2 = roeU0 - 0.5 * roeA;
  const Array1d sR1 = roeU0 + 0.5 * roeA;
  const Array1d sR2 = uR + aR;
  
  const Array1d sL = sL1.min(sL2);
  const Array1d sR = sR1.max(sR2);

  const Array1d sL_abs = sL.abs();
  const Array1d sR_abs = sR.abs();

  const Array1d maxS = sL_abs.max(sR_abs);
  const Array1d maxS_positive = mask.select(maxS, ones);

  const Array1d max_dt = Array1d::Constant(dx) / maxS_positive;
  const Array1d dummy_value = Array1d::Constant(std::numeric_limits<double>::max());
  const Array1d max_dt_result = mask.select(max_dt, dummy_value);
  return max_dt_result;
}

}