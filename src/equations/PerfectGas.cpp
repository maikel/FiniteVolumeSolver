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

#include "fub/equations/PerfectGas.hpp"
#include "fub/NewtonIteration.hpp"

namespace fub {
namespace {
template <int Dim>
double KineticEnergy(double density,
                     const Eigen::Array<double, Dim, 1>& momentum) noexcept {
  return 0.5 * momentum.matrix().squaredNorm() / density;
}

template <int Dim, int N, int O, int MR, int MC>
Array1d KineticEnergy(
    const Array1d& density,
    const Eigen::Array<double, Dim, N, O, MR, MC>& momentum) noexcept {
  Array1d squaredMomentum = momentum.matrix().colwise().squaredNorm();
  return 0.5 * squaredMomentum / density;
}

template <int Dim, int N, int O, int MR, int MC>
Array1d KineticEnergy(const Array1d& density,
                      const Eigen::Array<double, Dim, N, O, MR, MC>& momentum,
                      MaskArray mask) noexcept {
  mask = mask && (density > 0.0);
  Array1d squaredMomentum = momentum.matrix().colwise().squaredNorm();
  Array1d safe_density = density;
  safe_density = mask.select(density, 1.0);
  FUB_ASSERT((safe_density > 0.0).all());
  return 0.5 * squaredMomentum / safe_density;
}
} // namespace

template <int Dim>
void PerfectGas<Dim>::Flux(Conservative& flux, const Complete& state,
                           Direction dir) const noexcept {
  FUB_ASSERT(state.density > 0);
  const int d0 = static_cast<int>(dir);
  const double velocity = state.momentum[d0] / state.density;
  flux.density = state.momentum[d0];
  for (int d = 0; d < Dim; ++d) {
    flux.momentum[d] = velocity * state.momentum[d];
  }
  flux.momentum[d0] += state.pressure;
  flux.energy = velocity * (state.energy + state.pressure);
}

template <int Dim>
void PerfectGas<Dim>::Flux(ConservativeArray& flux, const CompleteArray& state,
                           Direction dir) const noexcept {
  const int d0 = static_cast<int>(dir);
  const Array1d velocity = state.momentum.row(d0) / state.density;
  flux.density = state.momentum.row(d0);
  for (int d = 0; d < Dim; ++d) {
    flux.momentum.row(d) = velocity * state.momentum.row(d);
  }
  flux.momentum.row(d0) += state.pressure;
  flux.energy = velocity * (state.energy + state.pressure);
}

template <int Dim>
void PerfectGas<Dim>::Flux(ConservativeArray& flux, const CompleteArray& state, MaskArray mask,
                           Direction dir) const noexcept {
  const int d0 = static_cast<int>(dir);
  mask = mask && (state.density > 0.0);
  const Array1d density = mask.select(state.density, 1.0);
  FUB_ASSERT((density > 0.0).all());
  const Array1d velocity = state.momentum.row(d0) / density;
  flux.density = state.momentum.row(d0);
  for (int d = 0; d < Dim; ++d) {
    flux.momentum.row(d) = velocity * state.momentum.row(d);
  }
  flux.momentum.row(d0) += state.pressure;
  flux.energy = velocity * (state.energy + state.pressure);
}


template <int Dim>
void PerfectGas<Dim>::CompleteFromCons(
    Complete& complete, const ConservativeBase<PerfectGas>& cons) const
    noexcept {
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  const double e_kin = KineticEnergy(cons.density, cons.momentum);
  FUB_ASSERT(e_kin < cons.energy);
  const double e_int = cons.energy - e_kin;
  complete.pressure = e_int / gamma_minus_1_inv;
  FUB_ASSERT(complete.pressure > 0.0 && complete.density > 0.0);
  complete.speed_of_sound =
      std::sqrt(gamma * complete.pressure / complete.density);
}

template <int Dim>
void PerfectGas<Dim>::CompleteFromCons(
    CompleteArray& complete,
    const ConservativeArrayBase<PerfectGas>& cons) const noexcept {
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  const Array1d e_kin = KineticEnergy(cons.density, cons.momentum);
  const Array1d e_int = cons.energy - e_kin;
  complete.pressure = e_int / gamma_minus_1_inv;
  complete.speed_of_sound =
      (gamma_array_ * complete.pressure / complete.density).sqrt();
}

template <int Dim>
void PerfectGas<Dim>::CompleteFromCons(
    CompleteArray& complete, const ConservativeArrayBase<PerfectGas>& cons,
    MaskArray mask) const noexcept {
  Array1d zero = Array1d::Zero();
  mask = mask && (cons.density > 0.0);
  complete.density = mask.select(cons.density, zero);
  for (int d = 0; d < Dim; ++d) {
    complete.momentum.row(d) = mask.select(cons.momentum.row(d), zero);
  }
  complete.energy = mask.select(cons.energy, zero);
  const Array1d e_kin = KineticEnergy(cons.density, cons.momentum, mask);
  const Array1d e_int = cons.energy - e_kin;
  complete.pressure = mask.select(e_int / gamma_minus_1_inv, zero);
  Array1d safe_density = mask.select(complete.density, 1.0);
  FUB_ASSERT((safe_density > 0.0).all());
  complete.speed_of_sound = mask.select(
      (gamma_array_ * complete.pressure / safe_density).sqrt(), zero);
}

template <int Dim>
Complete<PerfectGas<Dim>>
PerfectGas<Dim>::CompleteFromPrim(double rho, const Array<double, Dim, 1>& v,
                                  double p) const noexcept {
  Complete q{};
  q.density = rho;
  q.momentum = rho * v;
  q.pressure = p;
  const double e_kin = KineticEnergy(q.density, q.momentum);
  q.energy = e_kin + p * gamma_minus_1_inv;
  q.speed_of_sound = std::sqrt(gamma * q.pressure / q.density);
  return q;
}

template <int Dim>
Array<double, Dim, 1> PerfectGas<Dim>::Velocity(const Complete& complete) const
    noexcept {
  Array<double, Dim, 1> u = complete.momentum / complete.density;
  return u;
}

template <int Dim>
double PerfectGas<Dim>::Machnumber(const Complete& complete) const noexcept {
  double u = Velocity(complete).matrix().norm();
  return u / complete.speed_of_sound;
}

template struct PerfectGas<1>;
template struct PerfectGas<2>;
template struct PerfectGas<3>;

void Rotate(Conservative<PerfectGas<2>>& rotated,
            const Conservative<PerfectGas<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation, const PerfectGas<2>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Rotate(Complete<PerfectGas<2>>& rotated,
            const Complete<PerfectGas<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation, const PerfectGas<2>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.pressure = state.pressure;
  rotated.speed_of_sound = state.speed_of_sound;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Rotate(Conservative<PerfectGas<3>>& rotated,
            const Conservative<PerfectGas<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation, const PerfectGas<3>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Rotate(Complete<PerfectGas<3>>& rotated,
            const Complete<PerfectGas<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation, const PerfectGas<3>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.pressure = state.pressure;
  rotated.speed_of_sound = state.speed_of_sound;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Reflect(Complete<PerfectGas<1>>& reflected,
             const Complete<PerfectGas<1>>& state,
             const Eigen::Matrix<double, 1, 1>& normal, const PerfectGas<1>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
}

void Reflect(Complete<PerfectGas<2>>& reflected,
             const Complete<PerfectGas<2>>& state,
             const Eigen::Vector2d& normal, const PerfectGas<2>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
}

void Reflect(Complete<PerfectGas<3>>& reflected,
             const Complete<PerfectGas<3>>& state,
             const Eigen::Vector3d& normal, const PerfectGas<3>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
}

template <int Dim>
std::array<double, 2> EinfeldtSignalVelocities<PerfectGas<Dim>>::
operator()(const PerfectGas<Dim>&, const Complete& left, const Complete& right,
           Direction dir) const noexcept {
  FUB_ASSERT(left.density > 0.0);
  FUB_ASSERT(right.density > 0.0);
  const double rhoL = left.density;
  const double rhoR = right.density;
  double rhoUL;
  double rhoUR;
  rhoUL = left.momentum[int(dir)];
  rhoUR = right.momentum[int(dir)];
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

template <int Dim>
std::array<Array1d, 2> EinfeldtSignalVelocities<PerfectGas<Dim>>::
operator()(const PerfectGas<Dim>&, const CompleteArray& left,
           const CompleteArray& right, Direction dir) const noexcept {
  const Array1d rhoL = left.density;
  const Array1d rhoR = right.density;
  const Array1d rhoUL = left.momentum.row(int(dir));
  const Array1d rhoUR = right.momentum.row(int(dir));
  const Array1d aL = left.speed_of_sound;
  const Array1d aR = right.speed_of_sound;
  const Array1d sqRhoL = rhoL.sqrt();
  const Array1d sqRhoR = rhoR.sqrt();
  const Array1d uL = rhoUL / rhoL;
  const Array1d uR = rhoUR / rhoR;
  const Array1d roeU = (sqRhoL * uL + sqRhoR * uR) / (sqRhoL + sqRhoR);
  const Array1d roeA =
      ((sqRhoL * aL * aL + sqRhoR * aR * aR) / (sqRhoL + sqRhoR) +
       0.5 * (sqRhoL * sqRhoR) / ((sqRhoL + sqRhoR) * (sqRhoL + sqRhoR)) *
           (uR - uL) * (uR - uL))
          .sqrt();
  const Array1d sL1 = uL - aL;
  const Array1d sL2 = roeU - 0.5 * roeA;
  const Array1d sR1 = roeU + 0.5 * roeA;
  const Array1d sR2 = uR + aR;
  return {sL1.min(sL2), sR1.max(sR2)};
}

template <int Dim>
std::array<Array1d, 2> EinfeldtSignalVelocities<PerfectGas<Dim>>::
operator()(const PerfectGas<Dim>&, const CompleteArray& left,
           const CompleteArray& right, const MaskArray& mask,
           Direction dir) const noexcept {
  const Array1d rhoL = left.density;
  const Array1d rhoR = right.density;
  const Array1d zero = Array1d::Zero();
  const Array1d one = Array1d::Constant(1.0);
  const Array1d rhoUL = mask.select(left.momentum.row(int(dir)), zero);
  const Array1d rhoUR = mask.select(right.momentum.row(int(dir)), zero);
  const Array1d aL = mask.select(left.speed_of_sound, zero);
  const Array1d aR = mask.select(right.speed_of_sound, zero);
  const Array1d rhoLs = mask.select(rhoL, one);
  const Array1d rhoRs = mask.select(rhoR, one);
  const Array1d sqRhoL = rhoLs.sqrt();
  const Array1d sqRhoR = rhoRs.sqrt();
  const Array1d sqRho = sqRhoL + sqRhoR;
  FUB_ASSERT((rhoLs > 0.0).all());
  FUB_ASSERT((rhoRs > 0.0).all());
  FUB_ASSERT((sqRho > 0.0).all());
  const Array1d uL = rhoUL / rhoLs;
  const Array1d uR = rhoUR / rhoRs;
  const Array1d sqRhoL_over_sqRho = sqRhoL / sqRho;
  const Array1d sqRhoR_over_sqRho = sqRhoR / sqRho;
  const Array1d roeU = sqRhoL_over_sqRho * uL + sqRhoR_over_sqRho * uR;
  const Array1d roeA =
      (sqRhoL_over_sqRho * aL * aL + sqRhoR_over_sqRho * aR * aR  +
       Array1d::Constant(0.5) * sqRhoL_over_sqRho * sqRhoR_over_sqRho * (uR - uL) * (uR - uL))
          .sqrt();
  const Array1d sL1 = uL - aL;
  const Array1d sL2 = roeU - Array1d::Constant(0.5) * roeA;
  const Array1d sR1 = roeU + Array1d::Constant(0.5) * roeA;
  const Array1d sR2 = uR + aR;
  return {sL1.min(sL2), sR1.max(sR2)};
}

template struct EinfeldtSignalVelocities<PerfectGas<1>>;
template struct EinfeldtSignalVelocities<PerfectGas<2>>;
template struct EinfeldtSignalVelocities<PerfectGas<3>>;

template class FluxMethod<
    Hll<PerfectGas<3>, EinfeldtSignalVelocities<PerfectGas<3>>>>;

template <int Dim>
std::array<double, 2> ExactRiemannSolver<PerfectGas<Dim>>::ComputeMiddleState(
    const Complete& left, const Complete& right, Direction dir) {
  const double g = equation_.gamma;
  const double gp1 = g + 1.0;
  const double gm1 = g - 1.0;
  const double rhoL = left.density;
  const double pL = left.pressure;
  const double aL = left.speed_of_sound;
  const double A_l = 2.0 / (gp1 * rhoL);
  const double B_l = gm1 / gp1 * pL;
  const double exponent = 0.5 * gm1 / g;
  const double d_exponent = -0.5 * gp1 / g;
  auto f_left = [=](double p) {
    if (p > pL) {
      return (p - pL) * std::sqrt(A_l / (p + B_l));
    } else {
      return 2.0 * aL / gm1 * (std::pow(p / pL, exponent) - 1.0);
    }
  };
  auto df_left = [=](double p) {
    if (p > pL) {
      return std::sqrt(A_l / (B_l + p)) * (1.0 - 0.5 * (p - pL) / (B_l + p));
    } else {
      return std::pow(p / pL, d_exponent) / (rhoL * aL);
    }
  };

  const double rhoR = right.density;
  const double pR = right.pressure;
  const double aR = right.speed_of_sound;
  const double A_r = 2.0 / (gp1 * rhoR);
  const double B_r = gm1 / gp1 * pR;
  auto f_right = [=](double p) {
    if (p > pR) {
      return (p - pR) * std::sqrt(A_r / (p + B_r));
    } else {
      return 2 * aR / gm1 * (std::pow(p / pR, exponent) - 1.0);
    }
  };
  auto df_right = [=](double p) {
    if (p > pR) {
      return std::sqrt(A_r / (B_r + p)) * (1.0 - 0.5 * (p - pR) / (B_r + p));
    } else {
      return std::pow(p / pR, d_exponent) / (rhoR * aR);
    }
  };

  auto velocity = [dir](const Complete& state) {
    return state.momentum[static_cast<int>(dir)] / state.density;
  };
  const double uL = velocity(left);
  const double uR = velocity(right);
  const double du = uR - uL;
  auto f = [=](double h) { return f_right(h) + f_left(h) + du; };
  auto df = [=](double h) { return df_right(h) + df_left(h); };

  const double pM = NewtonIteration(f, df, 0.5 * (pL + pR));
  const double fL = f_left(pM);
  const double fR = f_right(pM);
  const double uM = 0.5 * (uL + uR) + 0.5 * (fR - fL);

  return {pM, uM};
}

template <int Dim>
std::array<double, 2> ExactRiemannSolver<PerfectGas<Dim>>::ComputeSignals(
    const Complete& left, const Complete& right, Direction dir) {
  auto velocity = [dir](const Complete& state) {
    return state.momentum[static_cast<int>(dir)] / state.density;
  };
  FUB_ASSERT(0.0 < left.density);
  FUB_ASSERT(0.0 < left.pressure);
  FUB_ASSERT(0.0 < right.density);
  FUB_ASSERT(0.0 < right.pressure);
  std::array<double, 2> states = ComputeMiddleState(left, right, dir);
  const double pM = states[0];
  FUB_ASSERT(0.0 < pM);
  const double g = equation_.gamma;
  const double gp1 = g + 1.0;
  const double gm1 = g - 1.0;

  std::array<double, 2> signals{};

  // LEFT SIDE

  const double uL = velocity(left);
  const double aL = left.speed_of_sound;
  const double pL = left.pressure;
  const double pL_ratio = pM / pL;
  if (pM > pL) {
    signals[0] = uL - aL * std::sqrt((gp1 * pL_ratio + gm1) / (2.0 * g));
  } else {
    signals[0] = uL - aL;
  }

  // RIGHT SIDE
  const double uR = velocity(right);
  const double aR = right.speed_of_sound;
  const double pR = right.pressure;
  const double pR_ratio = pM / pR;
  if (pM > pR) {
    signals[1] = uR + aR * std::sqrt((gp1 * pR_ratio + gm1) / (2.0 * g));
  } else {
    signals[1] = uR + aR;
  }

  return signals;
}

template <int Dim>
void ExactRiemannSolver<PerfectGas<Dim>>::SolveRiemannProblem(
    Complete& state, const Complete& left, const Complete& right,
    Direction dir) {
  auto velocity = [dir](const Complete& state) {
    return state.momentum[static_cast<int>(dir)] / state.density;
  };
  auto set_momentum = [dir](Complete& state, double momentum,
                            const Complete& side) {
    if constexpr (Dim == 1) {
      (void)dir;
      state.momentum = momentum;
    } else {
      state.momentum = state.density * (side.momentum / side.density);
      state.momentum[static_cast<int>(dir)] = momentum;
    }
  };
  auto from_prim = [eq = equation_](Complete& state) {
    const double rhoE_kin = KineticEnergy(state.density, state.momentum);
    const double rhoE_internal = state.pressure * eq.gamma_minus_1_inv;
    state.energy = rhoE_internal + rhoE_kin;
    state.speed_of_sound = std::sqrt(eq.gamma * state.pressure / state.density);
  };
  if (0.0 < left.density && 0.0 < right.density) {

  FUB_ASSERT(0.0 < left.density);
  FUB_ASSERT(0.0 < left.pressure);
  FUB_ASSERT(0.0 < right.density);
  FUB_ASSERT(0.0 < right.pressure);
  std::array<double, 2> states = ComputeMiddleState(left, right, dir);
  const double pM = states[0];
  FUB_ASSERT(0.0 < pM);
  const double uM = states[1];
  const double g = equation_.gamma;
  const double gp1 = g + 1.0;
  const double gm1 = g - 1.0;

  // (II) Determine Structure of Solution

  // Left side is interesting
  if (0.0 <= uM) {
    const double uL = velocity(left);
    const double aL = left.speed_of_sound;
    const double pL = left.pressure;
    const double p_ratio = pM / pL;
    const double g_ratio = gm1 / gp1;

    // SHOCK CASE
    if (pM > pL) {
      const double shock_speed =
          uL - aL * std::sqrt((gp1 * p_ratio + gm1) / (2.0 * g));
      if (shock_speed < 0.0) {
        state.density =
            left.density * (p_ratio + g_ratio) / (g_ratio * p_ratio + 1.0);
        set_momentum(state, state.density * uM, left);
        state.pressure = pM;
        from_prim(state);
      } else {
        state.density = left.density;
        state.momentum = left.momentum;
        state.energy = left.energy;
        state.pressure = left.pressure;
        state.speed_of_sound = left.speed_of_sound;
      }
    }

    // RAREFACTION CASE
    else {
      const double aM = aL * std::pow(p_ratio, gm1 / (2.0 * g));
      const double shock_speed_TL = uM - aM;
      const double shock_speed_HL = uL - aL;
      if (0.0 < shock_speed_HL) {
        state.density = left.density;
        state.momentum = left.momentum;
        state.energy = left.energy;
        state.pressure = left.pressure;
        state.speed_of_sound = left.speed_of_sound;
      } else if (shock_speed_TL < 0.0) {
        state.density = left.density * std::pow(p_ratio, 1 / g);
        set_momentum(state, state.density * uM, left);
        state.pressure = pM;
        from_prim(state);
      } else {
        state.density =
            left.density * std::pow(2.0 / gp1 + g_ratio * uL / aL, 2.0 / gm1);
        const double u = 2.0 / gp1 * (aL + 0.5 * gm1 * uL);
        set_momentum(state, state.density * u, left);
        state.pressure =
            pL * std::pow(2.0 / gp1 + g_ratio * uL / aL, 2.0 * g / gm1);
        from_prim(state);
      }
    }
  }

  // Right Hand Side is interesting
  else {
    const double uR = velocity(right);
    const double aR = right.speed_of_sound;
    const double pR = right.pressure;
    const double p_ratio = pM / pR;
    const double g_ratio = gm1 / gp1;

    // SHOCK CASE
    if (pM > pR) {
      const double shock_speed =
          uR + aR * std::sqrt((gp1 * p_ratio + gm1) / (2.0 * g));
      if (0.0 < shock_speed) {
        state.density =
            right.density * (p_ratio + g_ratio) / (g_ratio * p_ratio + 1.0);
        set_momentum(state, state.density * uM, right);
        state.pressure = pM;
        from_prim(state);
      } else {
        state.density = right.density;
        state.momentum = right.momentum;
        state.energy = right.energy;
        state.pressure = right.pressure;
        state.speed_of_sound = right.speed_of_sound;
      }
    }

    // RAREFACTION CASE
    else {
      const double aM = aR * std::pow(p_ratio, gm1 / (2.0 * g));
      const double shock_speed_TR = uM + aM;
      const double shock_speed_HR = uR + aR;
      if (shock_speed_HR <= 0.0) {
        state.density = right.density;
        state.momentum = right.momentum;
        state.energy = right.energy;
        state.pressure = right.pressure;
        state.speed_of_sound = right.speed_of_sound;
      } else if (shock_speed_TR >= 0.0) {
        state.density = right.density * std::pow(p_ratio, 1 / g);
        set_momentum(state, state.density * uM, right);
        state.pressure = pM;
        from_prim(state);
      } else {
        state.density =
            right.density * std::pow(2.0 / gp1 - g_ratio * uR / aR, 2.0 / gm1);
        const double u = 2.0 / gp1 * (-aR + 0.5 * gm1 * uR);
        set_momentum(state, state.density * u, right);
        state.pressure =
            pR * std::pow(2.0 / gp1 - g_ratio * uR / aR, 2.0 * g / gm1);
        from_prim(state);
      }
    }
  }

  }
}

template <int Dim>
void ExactRiemannSolver<PerfectGas<Dim>>::SolveRiemannProblem(
    CompleteArray& state, const CompleteArray& left, const CompleteArray& right,
    MaskArray mask, Direction dir) {
  for (int i = 0; i < mask.size(); ++i) {
    if (mask(i)) {
      typename PerfectGas<Dim>::Complete qL;
      qL.density = left.density(i);
      for (int d = 0; d < Dim; ++d) {
        qL.momentum(d) = left.momentum(d, i);
      }
      qL.energy = left.energy(i);
      qL.pressure = left.pressure(i);
      qL.speed_of_sound = left.speed_of_sound(i);

      typename PerfectGas<Dim>::Complete qR;
      qR.density = right.density(i);
      for (int d = 0; d < Dim; ++d) {
        qR.momentum(d) = right.momentum(d, i);
      }
      qR.energy = right.energy(i);
      qR.pressure = right.pressure(i);
      qR.speed_of_sound = right.speed_of_sound(i);

      typename PerfectGas<Dim>::Complete q;
      SolveRiemannProblem(q, qL, qR, dir);
      state.density(i) = q.density;
      for (int d = 0; d < Dim; ++d) {
        state.momentum(d, i) = q.momentum(d);
      }
      state.energy(i) = q.energy;
      state.pressure(i) = q.pressure;
      state.speed_of_sound(i) = q.speed_of_sound;
    }
  }
}

template class ExactRiemannSolver<PerfectGas<1>>;
template class ExactRiemannSolver<PerfectGas<2>>;
template class ExactRiemannSolver<PerfectGas<3>>;

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
  const double factor = b * delta * roeU0;
  const double squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  if constexpr (Dim == 1) {
    flux.density  = flux_hlle.density  - factor * 1.0;
    flux.momentum = flux_hlle.momentum - factor * roeU[0];
    flux.energy   = flux_hlle.energy   - factor * squaredNormRoeU_half;
  } else if constexpr (Dim == 2) {
    const int i = int(dir);
    const int j = int(i == 0);
    flux.density         = flux_hlle.density     - factor;
    flux.momentum[i] = flux_hlle.momentum[i] - factor *  roeU[i];
    flux.momentum[j] = flux_hlle.momentum[j] - factor * (roeU[j]              + 1.0);
    flux.energy          = flux_hlle.energy      - factor * (squaredNormRoeU_half + roeU[j]);
  } else {
    static_assert(Dim == 3);
    const int i = int(dir);
    const int j = (i + 1) % 3;
    const int k = (j + 1) % 3;
    flux.density     = flux_hlle.density     - factor;
    flux.momentum[i] = flux_hlle.momentum[i] - factor *  roeU[i];
    flux.momentum[j] = flux_hlle.momentum[j] - factor * (roeU[j]              + 1.0);
    flux.momentum[k] = flux_hlle.momentum[k] - factor * (roeU[k]                        + 1.0);
    flux.energy      = flux_hlle.energy      - factor * (squaredNormRoeU_half + roeU[j] + roeU[k]);
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
  const Array<double, Dim> deltaRhoV = right.momentum - left.momentum;

  const Array1d C1 = (deltaRhoU - deltaRho * roeU0) / roeA;

  const Array1d squaredNormRoeU_half = 0.5 * squaredNormRoeU;

  if constexpr (Dim == 1) {
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - C1 * roeU * roeA) / (roeH - squaredNormRoeU_half);
    const Array1d alpha_2 = deltaRho - C2;
    const Array1d factor = b * delta * alpha_2;

    flux.density  = flux_hlle.density  - factor * 1.0;
    flux.momentum = flux_hlle.momentum - factor * roeU;
    flux.energy   = flux_hlle.energy   - factor * squaredNormRoeU_half;
  } else if constexpr (Dim == 2) {
    const int i = int(dir);
    const int j = int(i == 0);
    const Array1d alpha_3 = deltaRhoV.row(j) - deltaRho * roeU.row(j);
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - alpha_3 * roeU.row(j) - C1 * roeU0 * roeA) / (roeH - squaredNormRoeU_half);
    const Array1d alpha_2 = deltaRho - C2;
    const Array1d factor = b * delta;
    flux.density         = flux_hlle.density         - factor * alpha_2 * 1.0;
    flux.momentum.row(i) = flux_hlle.momentum.row(i) - factor * alpha_2 * roeU.row(i);
    flux.momentum.row(j) = flux_hlle.momentum.row(j) - factor * (alpha_2 * roeU.row(j) + alpha_3 * 1.0);
    flux.energy          = flux_hlle.energy          - factor * (alpha_2 * squaredNormRoeU_half + alpha_3 * roeU.row(j));
  } else {
    static_assert(Dim == 3);
    const int i = int(dir);
    const int j = (i + 1) % 3;
    const int k = (j + 1) % 3;
    const Array1d alpha_3 = deltaRhoV.row(j) - deltaRho * roeU.row(j);
    const Array1d alpha_4 = deltaRhoV.row(k) - deltaRho * roeU.row(k);
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - alpha_3 * roeU.row(j) - alpha_4 * roeU.row(k) - C1 * roeU0 * roeA) / (roeH - squaredNormRoeU_half);
    const Array1d alpha_2 = deltaRho - C2;
    const Array1d factor = b * delta;
    flux.density         = flux_hlle.density         - factor * alpha_2;
    flux.momentum.row(i) = flux_hlle.momentum.row(i) - factor * alpha_2 * roeU.row(i);
    flux.momentum.row(j) = flux_hlle.momentum.row(j) - factor * (alpha_2 * roeU.row(j) + alpha_3 * 1.0);
    flux.momentum.row(k) = flux_hlle.momentum.row(k) - factor * (alpha_2 * roeU.row(k) + alpha_4 * 1.0);
    flux.energy          = flux_hlle.energy          - factor * (alpha_2 * squaredNormRoeU_half + alpha_3 * roeU.row(j) + alpha_4 * roeU.row(k));
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
  const Array1d factor = b * delta;

  const Array1d roeH_minus_kinetic = face_mask.select(roeH - squaredNormRoeU_half, ones);

  if constexpr (Dim == 1) {
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - C1 * roeU * roeA) / roeH_minus_kinetic;
    const Array1d alpha_2 = deltaRho - C2;
    const Array1d factor = b * delta * alpha_2;

    flux.density  = flux_hlle.density  - factor * 1.0;
    flux.momentum = flux_hlle.momentum - factor * roeU;
    flux.energy   = flux_hlle.energy   - factor * squaredNormRoeU_half;
  } else if constexpr (Dim == 2) {
    const int i = int(dir);
    const int j = int(i == 0);
    const Array1d alpha_3 = deltaRhoU.row(j) - deltaRho * roeU.row(j);
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - alpha_3 * roeU.row(j) - C1 * roeU0 * roeA) / roeH_minus_kinetic;
    const Array1d alpha_2 = deltaRho - C2;
    flux.density         = flux_hlle.density         - factor * alpha_2 * 1.0;
    flux.momentum.row(i) = flux_hlle.momentum.row(i) - factor * alpha_2 * roeU.row(i);
    flux.momentum.row(j) = flux_hlle.momentum.row(j) - factor * (alpha_2 * roeU.row(j) + alpha_3 * 1.0);
    flux.energy          = flux_hlle.energy          - factor * (alpha_2 * squaredNormRoeU_half + alpha_3 * roeU.row(j));
  } else {
    static_assert(Dim == 3);
    const int i = int(dir);
    const int j = (i + 1) % 3;
    const int k = (j + 1) % 3;
    const Array1d alpha_3 = deltaRhoU.row(j) - deltaRho * roeU.row(j);
    const Array1d alpha_4 = deltaRhoU.row(k) - deltaRho * roeU.row(k);
    const Array1d C2 = (deltaRhoE - deltaRho * squaredNormRoeU_half - alpha_3 * roeU.row(j) - alpha_4 * roeU.row(k) - C1 * roeU0 * roeA) / roeH_minus_kinetic;
    const Array1d alpha_2 = deltaRho - C2;
    flux.density         = flux_hlle.density         - factor * alpha_2;
    flux.momentum.row(i) = flux_hlle.momentum.row(i) - factor * alpha_2 * roeU.row(i);
    flux.momentum.row(j) = flux_hlle.momentum.row(j) - factor * (alpha_2 * roeU.row(j) + alpha_3 * 1.0);
    flux.momentum.row(k) = flux_hlle.momentum.row(k) - factor * (alpha_2 * roeU.row(k) + alpha_4 * 1.0);
    flux.energy          = flux_hlle.energy          - factor * (alpha_2 * squaredNormRoeU_half + alpha_3 * roeU.row(j) + alpha_4 * roeU.row(k));
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

// template <int Dim>
// void Hllem<Dim>::SolveRiemannProblem(Complete& solution, const Complete& left,
//                            const Complete& right, Duration dt, Direction dir) {
//   const double gm1 = equation_.gamma_ - 1.0;

//   // Compute Einfeldt signals velocities
    
//   const double rhoL = left.density;
//   const double rhoR = right.density;
//   const double rhoUL = left.momentum.row(int(dir));
//   const double rhoUR = right.momentum.row(int(dir));
//   const double aL = left.speed_of_sound;
//   const double aR = right.speed_of_sound;
//   const double hL = (left.energy + left.pressure) / rhoL;
//   const double hR = (right.energy + right.pressure) / rhoR;
//   const double sqRhoL = std::sqrt(rhoL);
//   const double sqRhoR = std::sqrt(rhoR);
//   const double sqRhoSum = sqRhoL + sqRhoR;
//   const double uL = rhoUL / rhoL;
//   const double uR = rhoUR / rhoR;
//   const double roeU0 = (sqRhoL * uL + sqRhoR * uR) / sqRhoSum;
//   const double roeH = (sqRhoL * hL + sqRhoR * hR) / sqRhoSum;
//   const double roeA2 = gm1 * (roeH  - 0.5 * roeU0 * roeU0);
//   const double roeA = std::sqrt(roeA2);
//   const double sL1 = uL - aL;
//   const double sL2 = roeU0 - 0.5 * roeA;
//   const double sR1 = roeU0 + 0.5 * roeA;
//   const double sR2 = uR + aR;
  
//   const double sL = std::min(sL1, sL2);
//   const double sR = std::max(sR1, sR2);

//   const double roeRho = sqRhoL * sqRhoR;
//   double roeU[Dim];
//   for (int i = 0; i < Dim; ++i) {
//     roeU[i] = (sqRhoL * uL[i] + sqRhoR * uR[i]) / sqRhoSum;
//   }
//   double squaredNormRoeU = 0.0;
//   for (int i = 0; i < Dim; ++i) {
//     squaredNormRoeU += roeU[i] * roeU[i];
//   }
//   const double u_bar = 0.5 * (sL + sR);
//   const double u_bar_abs = std::abs(u_bar);
//   const double delta = roeA / (roeA + u_bar_abs);
  
//   const double b = bLbR / db_positive;
//   const double factor = b * delta * roeU0;
//   const double squaredNormRoeU_half = 0.5 * squaredNormRoeU;

//   if (0 < sL) {
//     solution = left;
//   } else if (sR < 0) {
//     solution = right;
//   } else {
//     solution.density = 
//   }
// }

template struct Hllem<1>;
template struct Hllem<2>;
template struct Hllem<3>;

template class FluxMethod<Hllem<1>>;
template class FluxMethod<Hllem<2>>;
template class FluxMethod<Hllem<3>>;

} // namespace fub
