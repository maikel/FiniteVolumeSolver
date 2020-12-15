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

#ifndef FUB_FLUX_METHOD_GRADIENT
#define FUB_FLUX_METHOD_GRADIENT

#include "fub/Direction.hpp"
#include "fub/State.hpp"

namespace fub {
struct NoLimiter2 {
  double operator()(double) const noexcept { return 1.0; }
};

struct UpwindLimiter {
  double operator()(double) const noexcept { return 0.0; }
};

struct MinModLimiter {
  double operator()(double R) const noexcept {
    const double two_over_one_plus_R = 2.0 / (1.0 + R);
    return std::min(two_over_one_plus_R, R * two_over_one_plus_R);
  }
};

struct VanLeerLimiter {
  double operator()(double R) const noexcept {
    return 4.0 * R / ((1.0 + R) * (1.0 + R));
  }
};

template <typename Limiter> class CentralDifferenceGradient {
public:
  double ComputeGradient(double sL, double sR);

  template <typename State>
  void ComputeGradient(State& grad, const State& sL, const State& sR);

private:
  Limiter limiter_;
};

template <typename Equation,
          typename GradientMethod = CentralDifferenceGradient<MinModLimiter>>
class ConservativeGradient {
public:
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;

  using Gradient = Conservative;
  using GradientArray = ConservativeArray;

  ConservativeGradient() = default;

  explicit ConservativeGradient(const Equation& eq) : equation_{eq} {}

  explicit ConservativeGradient(const GradientMethod& method)
      : gradient_{method} {}

  ConservativeGradient(const Equation&, const GradientMethod& method)
      : gradient_{method} {}

  void ComputeGradient(Gradient& dudx, span<const Complete, 3> q, double dx,
                       Direction dir);

  void ComputeGradient(GradientArray& dudx, span<const CompleteArray, 3> q,
                       double dx, Direction dir);

private:
  Equation equation_{};
  Conservative dudx_L_{equation_};
  Conservative dudx_R_{equation_};
  GradientMethod gradient_{};
};

template <typename EulerEquation,
          typename GradientMethod = CentralDifferenceGradient<MinModLimiter>>
class PrimitiveGradient {
public:
  using Complete = ::fub::Complete<EulerEquation>;
  using CompleteArray = ::fub::CompleteArray<EulerEquation>;
  using Primitive = ::fub::Primitive<EulerEquation>;
  using PrimitiveArray = ::fub::PrimitiveArray<EulerEquation>;

  using Gradient = Primitive;
  using GradientArray = PrimitiveArray;

  PrimitiveGradient() = default;

  explicit PrimitiveGradient(const EulerEquation& equation)
      : equation_{equation} {}

  explicit PrimitiveGradient(const GradientMethod& method)
      : gradient_{method} {}

  PrimitiveGradient(const EulerEquation& equation, const GradientMethod& method)
      : equation_{equation}, gradient_{method} {}

  void ComputeGradient(Gradient& gradient, span<const Complete, 3> q, double dx,
                       Direction dir);

private:
  EulerEquation equation_{};
  Primitive wL_{equation_};
  Primitive wM_{equation_};
  Primitive wR_{equation_};
  Primitive dwdx_L_{equation_};
  Primitive dwdx_R_{equation_};
  GradientMethod gradient_{};
};

template <typename EulerEquation,
          typename GradientMethod = CentralDifferenceGradient<MinModLimiter>>
class CharacteristicsGradient {
public:
  using Complete = ::fub::Complete<EulerEquation>;
  using CompleteArray = ::fub::CompleteArray<EulerEquation>;
  using Characteristics = ::fub::Characteristics<EulerEquation>;
  using CharacteristicsArray = ::fub::CharacteristicsArray<EulerEquation>;
  using Primitive = ::fub::Primitive<EulerEquation>;
  using PrimitiveArray = ::fub::PrimitiveArray<EulerEquation>;

  using Gradient = Characteristics;
  using GradientArray = CharacteristicsArray;

  CharacteristicsGradient() = default;

  explicit CharacteristicsGradient(const EulerEquation& equation)
      : equation_{equation} {}

  explicit CharacteristicsGradient(const GradientMethod& method)
      : gradient_{method} {}

  CharacteristicsGradient(const EulerEquation& equation,
                          const GradientMethod& method)
      : equation_{equation}, gradient_{method} {}

  void ComputeGradient(Gradient& gradient, span<const Complete, 3> q, double dx,
                       Direction dir);

private:
  EulerEquation equation_{};
  Primitive wL_{equation_};
  Primitive wM_{equation_};
  Primitive wR_{equation_};
  Primitive dwdx_L_{equation_};
  Primitive dwdx_R_{equation_};
  Characteristics dKdx_L_{equation_};
  Characteristics dKdx_R_{equation_};
  GradientMethod gradient_{};
};

template <typename Equation, typename GradientMethod>
void ConservativeGradient<Equation, GradientMethod>::ComputeGradient(
    Gradient& dudx, span<const Complete, 3> q, double dx, Direction) {
  ForEachComponent(
      [dx](double& slope, double uL, double uM) { slope = (uM - uL) / dx; },
      dudx_L_, q[0], q[1]);
  ForEachComponent(
      [dx](double& slope, double uM, double uR) { slope = (uR - uM) / dx; },
      dudx_R_, q[1], q[2]);
  gradient_.ComputeGradient(dudx, dudx_L_, dudx_R_);
}

template <typename Equation, typename GradientMethod>
void PrimitiveGradient<Equation, GradientMethod>::ComputeGradient(
    Gradient& dwdx, span<const Complete, 3> q, double dx, Direction) {
  PrimFromComplete(equation_, wL_, q[0]);
  PrimFromComplete(equation_, wM_, q[1]);
  PrimFromComplete(equation_, wR_, q[2]);
  ForEachComponent(
      [dx](double& slope, double wL, double wR) { slope = (wR - wL) / dx; },
      dwdx_L_, wL_, wM_);
  ForEachComponent(
      [dx](double& slope, double wL, double wR) { slope = (wR - wL) / dx; },
      dwdx_R_, wM_, wR_);
  gradient_.ComputeGradient(dwdx, dwdx_L_, dwdx_R_);
}

template <typename Equation>
void ComputeAmplitudes(Characteristics<Equation>& amplitudes,
                       const Primitive<Equation>& left,
                       const Primitive<Equation>& right, double rhoc,
                       double ooc2, int ix) {
  const double dp = right.pressure - left.pressure;
  const double du = right.velocity[ix] - left.velocity[ix];
  const double drho = right.density - left.density;
  amplitudes.minus = dp - rhoc * du;
  amplitudes.plus = dp + rhoc * du;
  amplitudes.zero[0] = drho - ooc2 * dp;
  constexpr int Rank = Equation::Rank();
  const int n = amplitudes.zero.size();
  for (int i = 1; i < Rank; ++i) {
    const int iy = (ix + i) % Rank;
    amplitudes.zero[i] = right.velocity[iy] - left.velocity[iy];
  }
}

template <typename Equation, typename GradientMethod>
void CharacteristicsGradient<Equation, GradientMethod>::ComputeGradient(
    Gradient& dKdx, span<const Complete, 3> q, double dx, Direction dir) {
  PrimFromComplete(equation_, wL_, q[0]);
  PrimFromComplete(equation_, wM_, q[1]);
  PrimFromComplete(equation_, wR_, q[2]);
  const int ix = static_cast<int>(dir);
  const double rho = q[1].density;
  const double c = q[1].speed_of_sound;
  const double ooc2 = 1.0 / (c * c);
  const double rhoc = rho * c;
  ComputeAmplitudes(dKdx_L_, wL_, wM_, rhoc, ooc2, ix);
  ComputeAmplitudes(dKdx_R_, wM_, wR_, rhoc, ooc2, ix);
  gradient_.ComputeGradient(dKdx, dKdx_L_, dKdx_R_);
}

template <typename Limiter>
double CentralDifferenceGradient<Limiter>::ComputeGradient(double sL,
                                                           double sR) {
  if (sL * sR <= 0.0) {
    return 0.0;
  }
  const double R = sR / sL;
  return 0.5 * limiter_(R) * (sL + sR);
}

template <typename Limiter>
template <typename State>
void CentralDifferenceGradient<Limiter>::ComputeGradient(State& grad,
                                                         const State& sL,
                                                         const State& sR) {
  ForEachComponent([this](double& grad_u, double sL,
                          double sR) { grad_u = ComputeGradient(sL, sR); },
                   grad, sL, sR);
}

} // namespace fub

#endif