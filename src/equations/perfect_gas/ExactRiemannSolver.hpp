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

#include "fub/equations/perfect_gas/ExactRiemannSolver.hpp"
#include "fub/NewtonIteration.hpp"

namespace fub {
namespace {
template <int Dim>
double KineticEnergy_(double density,
                     const Eigen::Array<double, Dim, 1>& momentum) noexcept {
  return 0.5 * momentum.matrix().squaredNorm() / density;
}
} // namespace

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
    const double rhoE_kin = KineticEnergy_(state.density, state.momentum);
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

}