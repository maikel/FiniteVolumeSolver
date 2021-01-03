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

#include "fub/Direction.hpp"
#include "fub/State.hpp"
#include "fub/flux_method/Gradient.hpp"

namespace fub {

template <typename Equation> class ConservativeReconstruction {
public:
  using Conservative = ::fub::Conservative<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using Gradient = Conservative;

  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using GradientArray = ConservativeArray;

  explicit ConservativeReconstruction(const Equation& equation)
      : equation_{equation} {}

  void Reconstruct(Complete& reconstruction, const Complete& q0,
                   const Gradient& du_dx, Duration dt, double dx, Direction dir,
                   Side side) noexcept;

  void Reconstruct(CompleteArray& reconstruction, const CompleteArray& q0,
                   const GradientArray& du_dx, Duration dt, double dx,
                   Direction dir, Side side) noexcept;

private:
  Equation equation_;

  Complete q_left_{equation_};
  Complete q_right_{equation_};
  Conservative flux_left_{equation_};
  Conservative flux_right_{equation_};

  CompleteArray q_left_array_{equation_};
  CompleteArray q_right_array_{equation_};
  ConservativeArray flux_left_array_{equation_};
  ConservativeArray flux_right_array_{equation_};
};

template <typename EulerEquation> class PrimitiveReconstruction {
public:
  using Primitive = ::fub::Primitive<EulerEquation>;
  using Complete = ::fub::Complete<EulerEquation>;
  using Gradient = Primitive;

  using PrimitiveArray = ::fub::PrimitiveArray<EulerEquation>;
  using CompleteArray = ::fub::CompleteArray<EulerEquation>;
  using GradientArray = PrimitiveArray;

  explicit PrimitiveReconstruction(const EulerEquation& equation)
      : equation_{equation} {}

  void Reconstruct(Complete& reconstruction, const Complete& q0,
                   const Gradient& dw_dx, Duration dt, double dx, Direction dir,
                   Side side) noexcept;

  void Reconstruct(CompleteArray& reconstruction, const CompleteArray& q0,
                   const GradientArray& dw_dx, Duration dt, double dx,
                   Direction dir, Side side) noexcept;

private:
  EulerEquation equation_;

  Primitive w_rec_{equation_};
  Primitive w_{equation_};

  PrimitiveArray w_rec_array_{equation_};
  PrimitiveArray w_array_{equation_};
};

template <typename EulerEquation> class CharacteristicsReconstruction {
public:
  using Characteristics = ::fub::Characteristics<EulerEquation>;
  using Primitive = ::fub::Primitive<EulerEquation>;
  using Complete = ::fub::Complete<EulerEquation>;
  using Gradient = Characteristics;

  using CharacteristicsArray = ::fub::CharacteristicsArray<EulerEquation>;
  using PrimitiveArray = ::fub::PrimitiveArray<EulerEquation>;
  using CompleteArray = ::fub::CompleteArray<EulerEquation>;
  using GradientArray = CharacteristicsArray;

  explicit CharacteristicsReconstruction(const EulerEquation& equation)
      : equation_{equation} {}

  void Reconstruct(Complete& reconstruction, const Complete& q0,
                   const Gradient& dw_dx, Duration dt, double dx, Direction dir,
                   Side side) noexcept;

  void Reconstruct(CompleteArray& reconstruction, const CompleteArray& q0,
                   const GradientArray& dw_dx, Duration dt, double dx,
                   Direction dir, Side side) noexcept;

private:
  EulerEquation equation_;

  Primitive w_rec_{equation_};
  Primitive dwdt_{equation_};
  Characteristics dKdt_{equation_};

  PrimitiveArray w_rec_array_{equation_};
  PrimitiveArray dwdt_array_{equation_};
  CharacteristicsArray dKdt_array_{equation_};
};

template <typename Equation>
void ConservativeReconstruction<Equation>::Reconstruct(
    Complete& reconstruction, const Complete& q0, const Gradient& dq_dx,
    Duration dt, double dx, Direction dir, Side side) noexcept {
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
  const Complete& q = side == Side::Lower ? q_left_ : q_right_;
  ForEachComponent(
      [&lambda_half](double& rec, double u, double fL, double fR) {
        const double dF = fL - fR;
        const double dU = lambda_half * dF;
        rec = u + dU;
      },
      AsCons(reconstruction), AsCons(q), flux_left_, flux_right_);
  CompleteFromCons(equation_, reconstruction, reconstruction);
}

template <typename Equation>
void ConservativeReconstruction<Equation>::Reconstruct(
    CompleteArray& reconstruction, const CompleteArray& q0,
    const GradientArray& dq_dx, Duration dt, double dx, Direction dir,
    Side side) noexcept {
  ForEachComponent(
      [dx](auto&& uL, auto&& uR, auto&& u, auto&& du_dx) {
        const Array1d du = 0.5 * dx * du_dx;
        uL = u - du;
        uR = u + du;
      },
      AsCons(q_left_array_), AsCons(q_right_array_), q0, dq_dx);
  CompleteFromCons(equation_, q_left_array_, q_left_array_);
  CompleteFromCons(equation_, q_right_array_, q_right_array_);
  Flux(equation_, flux_left_array_, q_left_array_, dir);
  Flux(equation_, flux_right_array_, q_right_array_, dir);
  const double lambda_half = 0.5 * dt.count() / dx;
  const CompleteArray& q = side == Side::Lower ? q_left_array_ : q_right_array_;
  ForEachComponent(
      [&lambda_half](auto&& rec, auto&& u, auto&& fL, auto&& fR) {
        const Array1d dF = fL - fR;
        const Array1d dU = lambda_half * dF;
        rec = u + dU;
      },
      AsCons(reconstruction), AsCons(q), flux_left_array_, flux_right_array_);
  CompleteFromCons(equation_, reconstruction, reconstruction);
}

template <typename EulerEquation>
void PrimitiveReconstruction<EulerEquation>::Reconstruct(
    Complete& reconstruction, const Complete& q0, const Primitive& dw_dx,
    Duration dt, double dx, Direction dir, Side side) noexcept {
  PrimFromComplete(equation_, w_, q0);
  const int ix = static_cast<int>(dir);
  const int sign = (side == Side::Upper) - (side == Side::Lower);
  const double lambda = sign * dt.count() / dx;
  const double dx_half = sign * 0.5 * dx;
  const double u = w_.velocity[ix];
  const double rho = q0.density;
  const double a = q0.speed_of_sound;
  const double a2 = a * a;
  const double rho_a2 = rho * a2;
  w_rec_.density =
      w_.density +
      dx_half * (dw_dx.density -
                 lambda * (u * dw_dx.density + rho * dw_dx.velocity[ix]));
  w_rec_.pressure =
      w_.pressure +
      dx_half * (dw_dx.pressure -
                 lambda * (rho_a2 * dw_dx.velocity[ix] + u * dw_dx.pressure));
  w_rec_.velocity[ix] =
      w_.velocity[ix] +
      dx_half * (dw_dx.velocity[ix] +
                 lambda * (u * dw_dx.velocity[ix] + dw_dx.pressure / rho));
  constexpr int Rank = EulerEquation::Rank();
  for (int i = 1; i < Rank; ++i) {
    const int iy = (ix + i) % Rank;
    w_rec_.velocity[iy] =
        w_.velocity[iy] +
        dx_half * (dw_dx.velocity[iy] - lambda * u * dw_dx.velocity[iy]);
  }
  if constexpr (fub::euler::state_with_species<EulerEquation, Primitive>()) {
    for (int i = 0; i < w_rec_.species.size(); ++i) {
      w_rec_.species[i] =
          w_.species[i] +
          dx_half * (dw_dx.species[i] - lambda * u * dw_dx.species[i]);
    }
  }
  CompleteFromPrim(equation_, reconstruction, w_rec_);
}

template <typename EulerEquation>
void PrimitiveReconstruction<EulerEquation>::Reconstruct(
    CompleteArray& reconstruction, const CompleteArray& q0,
    const PrimitiveArray& dw_dx, Duration dt, double dx, Direction dir,
    Side side) noexcept {
  PrimFromComplete(equation_, w_array_, q0);
  const int ix = static_cast<int>(dir);
  const int sign = (side == Side::Upper) - (side == Side::Lower);
  const Array1d lambda = Array1d::Constant(sign * dt.count() / dx);
  const Array1d dx_half = Array1d::Constant(sign * 0.5 * dx);
  const Array1d u = w_array_.velocity.row(ix);
  const Array1d rho = q0.density;
  const Array1d a = q0.speed_of_sound;
  const Array1d a2 = a * a;
  const Array1d rho_a2 = rho * a2;
  w_rec_array_.density =
      w_array_.density +
      dx_half * (dw_dx.density -
                 lambda * (u * dw_dx.density + rho * dw_dx.velocity.row(ix)));
  w_rec_array_.pressure =
      w_array_.pressure +
      dx_half * (dw_dx.pressure - lambda * (rho_a2 * dw_dx.velocity.row(ix) +
                                            u * dw_dx.pressure));
  w_rec_array_.velocity.row(ix) =
      w_array_.velocity.row(ix) +
      dx_half * (dw_dx.velocity.row(ix) +
                 lambda * (u * dw_dx.velocity.row(ix) + dw_dx.pressure / rho));
  constexpr int Rank = EulerEquation::Rank();
  for (int i = 1; i < Rank; ++i) {
    const int iy = (ix + i) % Rank;
    w_rec_array_.velocity.row(iy) =
        w_array_.velocity.row(iy) + dx_half * (dw_dx.velocity.row(iy) -
                                         lambda * u * dw_dx.velocity.row(iy));
  }
  if constexpr (fub::euler::state_with_species<EulerEquation, Primitive>()) {
    for (int i = 0; i < w_rec_array_.species.rows(); ++i) {
      w_rec_array_.species.row(i) =
          w_array_.species.row(i) +
          dx_half * (dw_dx.species.row(i) - lambda * u * dw_dx.species.row(i));
    }
  }
  CompleteFromPrim(equation_, reconstruction, w_rec_array_);
}

template <typename EulerEquation>
void CharacteristicsReconstruction<EulerEquation>::Reconstruct(
    Complete& reconstruction, const Complete& q0, const Characteristics& dKdx,
    Duration dt, double dx, Direction dir, Side side) noexcept {
  PrimFromComplete(equation_, w_rec_, q0);
  const int ix = static_cast<int>(dir);
  const double u = w_rec_.velocity[ix];
  const double c = q0.speed_of_sound;
  const int sign = (side == Side::Upper) - (side == Side::Lower);
  const double lambda = sign * dt.count() / dx;
  const double dx_half = sign * 0.5 * dx;
  dKdt_.minus = dx_half * dKdx.minus * (1.0 - lambda * (u - c));
  dKdt_.plus = dx_half * dKdx.plus * (1.0 - lambda * (u + c));
  for (int i = 0; i < dKdt_.zero.size(); ++i) {
    dKdt_.zero[i] = dx_half * dKdx.zero[i] * (1.0 - lambda * u);
  }
  if constexpr (fub::euler::state_with_species<EulerEquation, Primitive>()) {
    for (int i = 0; i < dKdt_.species.size(); ++i) {
      dKdt_.species[i] = dx_half * dKdx.species[i] * (1.0 - lambda * u);
    }
  }
  const double rho = w_rec_.density;
  const double rhoc = rho * c;
  const double c2 = c * c;
  dwdt_.pressure = 0.5 * (dKdt_.minus + dKdt_.plus);
  dwdt_.density = dwdt_.pressure / c2 + dKdt_.zero[0];
  dwdt_.velocity.row(ix) = 0.5 * (dKdt_.plus - dKdt_.minus) / rhoc;
  constexpr int Rank = EulerEquation::Rank();
  for (int i = 1; i < Rank; ++i) {
    const int iy = (ix + i) % Rank;
    dwdt_.velocity.row(iy) = dKdt_.zero[i];
  }
  if constexpr (fub::euler::state_with_species<EulerEquation, Primitive>()) {
    for (int i = 0; i < dKdt_.species.size(); ++i) {
      dwdt_.species[i] = dKdt_.species[i];
    }
  }
  w_rec_ += dwdt_;
  CompleteFromPrim(equation_, reconstruction, w_rec_);
}

template <typename EulerEquation>
void CharacteristicsReconstruction<EulerEquation>::Reconstruct(
    CompleteArray& reconstruction, const CompleteArray& q0,
    const CharacteristicsArray& dKdx, Duration dt, double dx, Direction dir,
    Side side) noexcept {
  PrimFromComplete(equation_, w_rec_array_, q0);
  const int ix = static_cast<int>(dir);
  const Array1d u = w_rec_array_.velocity.row(ix);
  const Array1d c = q0.speed_of_sound;
  const int sign = (side == Side::Upper) - (side == Side::Lower);
  const Array1d lambda = Array1d::Constant(sign * dt.count() / dx);
  const Array1d dx_half = Array1d::Constant(sign * 0.5 * dx);
  dKdt_array_.minus =
      dx_half * dKdx.minus * (Array1d::Constant(1.0) - lambda * (u - c));
  dKdt_array_.plus =
      dx_half * dKdx.plus * (Array1d::Constant(1.0) - lambda * (u + c));
  for (int i = 0; i < dKdt_array_.zero.rows(); ++i) {
    dKdt_array_.zero.row(i) =
        dx_half * dKdx.zero.row(i) * (Array1d::Constant(1.0) - lambda * u);
  }
  if constexpr (fub::euler::state_with_species<EulerEquation, Primitive>()) {
    for (int i = 0; i < dKdt_array_.species.rows(); ++i) {
      dKdt_array_.species.row(i) =
          dx_half * dKdx.species.row(i) * (Array1d::Constant(1.0) - lambda * u);
    }
  }
  const Array1d rho = w_rec_array_.density;
  const Array1d rhoc = rho * c;
  const Array1d c2 = c * c;
  dwdt_array_.pressure = 0.5 * (dKdt_array_.minus + dKdt_array_.plus);
  dwdt_array_.density = dwdt_array_.pressure / c2 + dKdt_array_.zero.row(0);
  dwdt_array_.velocity.row(ix) = 0.5 * (dKdt_array_.plus - dKdt_array_.minus) / rhoc;
  constexpr int Rank = EulerEquation::Rank();
  for (int i = 1; i < Rank; ++i) {
    const int iy = (ix + i) % Rank;
    dwdt_array_.velocity.row(iy) = dKdt_array_.zero.row(i);
  }
  if constexpr (fub::euler::state_with_species<EulerEquation, Primitive>()) {
    for (int i = 0; i < dKdt_array_.species.rows(); ++i) {
      dwdt_array_.species.row(i) = dKdt_array_.species.row(i);
    }
  }
  w_rec_array_ += dwdt_array_;
  CompleteFromPrim(equation_, reconstruction, w_rec_array_);
}

} // namespace fub

#endif