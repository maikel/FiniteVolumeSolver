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

#ifndef FUB_FLUX_METHOD_RECONSTRUCT
#define FUB_FLUX_METHOD_RECONSTRUCT

#include "fub/State.hpp"
#include "fub/flux_method/Gradient.hpp"

namespace fub {

template <typename Equation> class ConservativeReconstruction {
public:
  using Conservative = ::fub::Conservative<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using Gradient = Conservative;

  explicit ConservativeReconstruction(const Equation& equation)
      : equation_{equation} {}

  void Reconstruct(Complete& reconstruction, const Complete& q0,
                   const Gradient& du_dx, Duration dt, double dx,
                   Direction dir) noexcept;

private:
  Equation equation_;

  Complete q_left_{equation_};
  Complete q_right_{equation_};
  Conservative flux_left_{equation_};
  Conservative flux_right_{equation_};
};

template <typename EulerEquation> class PrimitiveReconstruction {
public:
  using Primitive = ::fub::Primitive<EulerEquation>;
  using Complete = ::fub::Complete<EulerEquation>;
  using Gradient = Primitive;

  explicit PrimitiveReconstruction(const EulerEquation& equation)
      : equation_{equation} {}

  void Reconstruct(Complete& reconstruction, const Complete& q0,
                   const Gradient& dw_dx, Duration dt, double dx,
                   Direction dir) noexcept;

private:
  EulerEquation equation_;

  Primitive w_rec_{equation_};
  Primitive w_{equation_};
};

template <typename EulerEquation> class CharacteristicsReconstruction {
public:
  using Characteristics = ::fub::Characteristics<EulerEquation>;
  using Primitive = ::fub::Primitive<EulerEquation>;
  using Complete = ::fub::Complete<EulerEquation>;
  using Gradient = Characteristics;

  explicit CharacteristicsReconstruction(const EulerEquation& equation)
      : equation_{equation} {}

  void Reconstruct(Complete& reconstruction, const Complete& q0,
                   const Gradient& dw_dx, Duration dt, double dx,
                   Direction dir) noexcept;

private:
  EulerEquation equation_;

  Primitive w_rec_{equation_};
  Primitive dwdt_{equation_};
  Characteristics dKdt_{equation_};
};

template <typename Equation>
void ConservativeReconstruction<Equation>::Reconstruct(Complete& reconstruction,
                                                       const Complete& q0,
                                                       const Gradient& dq_dx,
                                                       Duration dt, double dx,
                                                       Direction dir) noexcept {
  ForEachComponent(
      [dx](double& uL, double& uR, double u, double du_dx) {
        const double du = 0.5 * dx * du_dx;
        uL = u - du;
        uR = u + du;
      },
      AsCons(q_left_), AsCons(q_right_), q0, dq_dx);
  CompleteFromCons(equation_, q_left_, q_left_);
  CompleteFromCons(equation_, q_right_, q_right_);
  Flux(equation_, flux_left_, q_left_, dir);
  Flux(equation_, flux_right_, q_right_, dir);
  const double lambda_half = 0.5 * dt.count() / dx;
  ForEachComponent(
      [&lambda_half](double& rec, double u, double fL, double fR) {
        const double dF = fL - fR;
        const double dU = lambda_half * dF;
        rec = u + dU;
      },
      AsCons(reconstruction), AsCons(q_left_), flux_left_, flux_right_);
  CompleteFromCons(equation_, reconstruction, reconstruction);
}

template <typename EulerEquation>
void PrimitiveReconstruction<EulerEquation>::Reconstruct(
    Complete& reconstruction, const Complete& q0, const Primitive& dw_dx,
    Duration dt, double dx, Direction dir) noexcept {
  PrimFromComplete(equation_, w_, q0);
  const int ix = static_cast<int>(dir);
  const double lambda = dt.count() * dx;
  const double u = w_.velocity[ix];
  const double rho = q0.density;
  const double a = q0.speed_of_sound;
  const double a2 = a * a;
  const double rho_a2 = rho * a2;
  w_rec_.density =
      w_.density + 0.5 * dx *
                       (dw_dx.density + lambda * (u * dw_dx.density +
                                                  rho * dw_dx.velocity[ix]));
  w_rec_.pressure =
      w_.pressure +
      0.5 * dx *
          (dw_dx.pressure +
           lambda * (rho_a2 * dw_dx.velocity[ix] + u * dw_dx.pressure));
  w_rec_.velocity[ix] =
      w_.velocity[ix] +
      0.5 * dx *
          (dw_dx.velocity[ix] +
           lambda * (u * dw_dx.velocity[ix] + dw_dx.pressure / rho));
  constexpr int Rank = EulerEquation::Rank();
  for (int i = 1; i < Rank; ++i) {
    const int iy = (ix + i) % Rank;
    w_rec_.velocity[iy] = w_.velocity[iy] + 0.5 * dx * u * dw_dx.velocity[iy];
  }
  CompleteFromPrim(equation_, reconstruction, w_rec_);
}

template <typename EulerEquation>
void CharacteristicsReconstruction<EulerEquation>::Reconstruct(
    Complete& reconstruction, const Complete& q0, const Characteristics& dKdx,
    Duration dt, double dx, Direction dir) noexcept {
  PrimFromComplete(equation_, w_rec_, q0);
  const int ix = static_cast<int>(dir);
  const double u = w_rec_.velocity[ix];
  const double c = q0.speed_of_sound;
  const double lambda = dt.count() / dx;
  dKdt_.minus = 0.5 * dx * dKdx.minus * (1.0 + lambda * (u - c));
  dKdt_.plus = 0.5 * dx * dKdx.plus * (1.0 + lambda * (u + c));
  for (int i = 0; i < dKdt_.zero.size(); ++i) {
    dKdt_.zero[i] = 0.5 * dx * dKdx.zero[i] * (1.0 + lambda * u);
  }
  const double rho = w_rec_.density;
  const double rhoc = rho * c;
  const double c2 = c * c;
  dwdt_.pressure = 0.5 * (dKdt_.minus + dKdt_.plus);
  dwdt_.density = dwdt_.pressure / c2 + dKdt_.zero;
  dwdt_.velocity.row(ix) = 0.5 * (dKdt_.plus - dKdt_.minus) / rhoc;
  constexpr int Rank = EulerEquation::Rank();
  for (int i = 1; i < Rank; ++i) {
    const int iy = (ix + i) % Rank;
    dwdt_.velocity.row(iy) = dKdt_.zero[i];
  }
  w_rec_ += dwdt_;
  CompleteFromPrim(equation_, reconstruction, w_rec_);
}

} // namespace fub

#endif