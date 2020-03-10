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

#include "fub/equations/perfect_gas/MusclHancockPrimMethod.hpp"
#include "fub/equations/perfect_gas/HllemMethod.hpp"

namespace fub::perfect_gas {

template <int Dim>
double
MusclHancockPrim<Dim>::ComputeStableDt(span<const Complete, 4>,
                                       double, Direction) {
  throw std::runtime_error{"Sequential version is not implemented"};
  return 0.0;
}

template <int Dim>
Array1d
MusclHancockPrim<Dim>::ComputeStableDt(span<const CompleteArray, 4> states,
                                       double dx, Direction dir) {
  Hllem hllem{equation_};
  return hllem.ComputeStableDt(states.template subspan<1, 2>(), dx, dir);
}

template <int Dim>
Array1d MusclHancockPrim<Dim>::ComputeStableDt(
    span<const CompleteArray, 4> states, Array1d face_fraction,
    span<const Array1d, 4> vols, double dx, Direction dir) {
  Hllem hllem{equation_};
  return hllem.ComputeStableDt(states.template subspan<1, 2>(), face_fraction,
                               vols.subspan<1, 2>(), dx, dir);
}


template <int Dim>
void MusclHancockPrim<Dim>::ComputeNumericFlux(
    Conservative&, span<const Complete, 4>, Duration ,
    double, Direction) {
  throw std::runtime_error{"Sequential version is not implemented"};
}

template <int Dim>
void MusclHancockPrim<Dim>::ComputeNumericFlux(
    ConservativeArray& flux, span<const CompleteArray, 4> stencil, Duration dt,
    double dx, Direction dir) {
  auto limiter = [](Array1d qL, Array1d qM, Array1d qR) -> Array1d {
    Array1d delta_q = 0.5 * (qR - qL);
    Array1d zeros = Array1d::Zero();
    Array1d ones = Array1d::Constant(1.0);
    MaskArray is_relevant = delta_q.abs() > 1e-12;
    Array1d delta_q_if_relevant = is_relevant.select(delta_q, ones);
    Array1d rL = is_relevant.select((qM - qL) / delta_q_if_relevant, zeros);
    Array1d rR = is_relevant.select((qR - qM) / delta_q_if_relevant, zeros);
    MaskArray is_positive = rL > 0.0 && rR > 0.0;
    Array1d r = is_positive.select(rL.min(rR), zeros);
    MaskArray less_than_one = r < 1.0;
    Array1d r_1 = 2 * r / (Array1d::Constant(1.0) + r);
    Array1d r_2 = 0.5 * (r + Array1d::Constant(1.0));
    Array1d sigma = less_than_one.select(r_1, r_2);
    return sigma * delta_q;
  };
  Array1d rho_ll = stencil[0].density;
  Array1d rho_l = stencil[1].density;
  Array1d rho_r = stencil[2].density;
  Array1d rho_rr = stencil[3].density;

  Array1d p_ll = stencil[0].pressure;
  Array1d p_l = stencil[1].pressure;
  Array1d p_r = stencil[2].pressure;
  Array1d p_rr = stencil[3].pressure;

  Array1d a_l = stencil[1].speed_of_sound;
  Array1d a_r = stencil[2].speed_of_sound;

  Array<double, Dim> u_ll;
  Array<double, Dim> u_l;
  Array<double, Dim> u_r;
  Array<double, Dim> u_rr;
  for (int i = 0; i < Dim; ++i) {
    u_ll.row(i) = stencil[0].momentum.row(i) / stencil[0].density;
    u_l.row(i) = stencil[1].momentum.row(i) / stencil[1].density;
    u_r.row(i) = stencil[2].momentum.row(i) / stencil[2].density;
    u_rr.row(i) = stencil[3].momentum.row(i) / stencil[3].density;
  }

  int d = static_cast<int>(dir);
  Array1d lambda = Array1d::Constant(dt.count() / dx);

  Array1d drho_l = limiter(rho_ll, rho_l, rho_r);
  Array1d dp_l = limiter(p_ll, p_l, p_r);
  Array<double, Dim> du_l;
  for (int i = 0; i < Dim; ++i) {
    du_l.row(i) = limiter(u_ll.row(i), u_l.row(i), u_r.row(i));
  }

  Array1d drho_r = limiter(rho_l, rho_r, rho_rr);
  Array1d dp_r = limiter(p_l, p_r, p_rr);
  Array<double, Dim> du_r;
  for (int i = 0; i < Dim; ++i) {
    du_r.row(i) = limiter(u_l.row(i), u_r.row(i), u_rr.row(i));
  }

  Array1d rec_rho_l =
      rho_l +
      0.5 * (drho_l - lambda * (u_l.row(d) * drho_l + rho_l * du_l.row(d)));
  Array1d rec_p_l =
      p_l + 0.5 * (dp_l - lambda * (rho_l * a_l * a_l * du_l.row(d) +
                                    u_l.row(d) * dp_l));
  Array<double, Dim> rec_u_l;
  for (int i = 0; i < Dim; ++i) {
    rec_u_l.row(i) =
        u_l.row(i) + 0.5 * (du_l.row(i) - lambda * (u_l.row(d) * du_l.row(i)));
  }
  rec_u_l.row(d) -= 0.5 * lambda * dp_l / rho_l;

  Array1d rec_rho_r =
      rho_r -
      0.5 * (drho_r + lambda * (u_r.row(d) * drho_r + rho_r * du_r.row(d)));
  Array1d rec_p_r =
      p_r - 0.5 * (dp_r + lambda * (rho_r * a_r * a_r * du_r.row(d) +
                                    u_r.row(d) * dp_r));
  Array<double, Dim> rec_u_r;
  for (int i = 0; i < Dim; ++i) {
    rec_u_r.row(i) =
        u_r.row(i) - 0.5 * (du_r.row(i) + lambda * (u_r.row(d) * du_r.row(i)));
  }
  rec_u_r.row(d) -= 0.5 * lambda * dp_r / rho_r;

  auto squaredNorm = [](Array<double, Dim> u) -> Array1d {
    Array1d norm2 = u.matrix().colwise().squaredNorm();
    return norm2;
  };

  std::array<CompleteArray, 2> w;
  w[0].density = rec_rho_l;
  for (int i = 0; i < Dim; ++i) {
    w[0].momentum.row(i) = rec_rho_l * rec_u_l.row(i);
  }
  w[0].pressure = rec_p_l;
  w[0].speed_of_sound = (equation_.gamma_array_ * rec_p_l / rec_rho_l).sqrt();
  w[0].energy = 0.5 * rec_rho_l * squaredNorm(rec_u_l) +
                rec_p_l * equation_.gamma_minus_1_inv_array_;

  FUB_ASSERT((w[0].density > 0.0).all());
  FUB_ASSERT((w[0].pressure > 0.0).all());

  w[1].density = rec_rho_r;
  for (int i = 0; i < Dim; ++i) {
    w[1].momentum.row(i) = rec_rho_r * rec_u_r.row(i);
  }
  w[1].pressure = rec_p_r;
  w[1].speed_of_sound = (equation_.gamma_array_ * rec_p_r / rec_rho_r).sqrt();
  w[1].energy = 0.5 * rec_rho_r * squaredNorm(rec_u_r) +
                rec_p_r * equation_.gamma_minus_1_inv_array_;

  FUB_ASSERT((w[1].density > 0.0).all());
  FUB_ASSERT((w[1].pressure > 0.0).all());

  Hllem hllem{equation_};
  hllem.ComputeNumericFlux(flux, span{w}, dt, dx, dir);

  FUB_ASSERT((!isnan(flux.density)).all());
  FUB_ASSERT((!isnan(flux.momentum)).all());
  FUB_ASSERT((!isnan(flux.energy)).all());
}

template <int Dim>
void MusclHancockPrim<Dim>::ComputeNumericFlux(
    ConservativeArray& flux, Array1d face_fractions,
    span<const CompleteArray, 4> stencil,
    span<const Array1d, 4> volume_fractions, Duration dt, double dx,
    Direction dir) {
  MaskArray mask = face_fractions > 0.0;
  if (!mask.any()) {
    return;
  }
  auto limiter = [](Array1d qL, Array1d qM, Array1d qR) -> Array1d {
    Array1d delta_q = 0.5 * (qR - qL);
    Array1d zeros = Array1d::Zero();
    Array1d ones = Array1d::Constant(1.0);
    MaskArray is_relevant = delta_q.abs() > 1e-12;
    Array1d delta_q_if_relevant = is_relevant.select(delta_q, ones);
    Array1d rL = is_relevant.select((qM - qL) / delta_q_if_relevant, zeros);
    Array1d rR = is_relevant.select((qR - qM) / delta_q_if_relevant, zeros);
    MaskArray is_positive = rL > 0.0 && rR > 0.0;
    Array1d r = is_positive.select(rL.min(rR), zeros);
    MaskArray less_than_one = r < 1.0;
    Array1d r_1 = 2 * r / (Array1d::Constant(1.0) + r);
    Array1d r_2 = 0.5 * (r + Array1d::Constant(1.0));
    Array1d sigma = less_than_one.select(r_1, r_2);
    return sigma * delta_q;
  };

  MaskArray left_mask = volume_fractions[0] > 0.0;
  MaskArray right_mask = volume_fractions[3] > 0.0;
  Array1d zero = Array1d::Zero();
  Array1d one = Array1d::Constant(1.0);

  Array1d rho_l = mask.select(stencil[1].density, one);
  Array1d rho_r = mask.select(stencil[2].density, one);
  Array1d rho_ll = left_mask.select(stencil[0].density, rho_l);
  Array1d rho_rr = right_mask.select(stencil[3].density, rho_r);

  Array<double, Dim> rho_u_ll = stencil[0].momentum;
  Array<double, Dim> rho_u_l = stencil[1].momentum;
  Array<double, Dim> rho_u_r = stencil[2].momentum;
  Array<double, Dim> rho_u_rr = stencil[3].momentum;

  Array1d p_l = mask.select(stencil[1].pressure, one);
  Array1d p_r = mask.select(stencil[2].pressure, one);
  Array1d p_ll = left_mask.select(stencil[0].pressure, p_l);
  Array1d p_rr = right_mask.select(stencil[3].pressure, p_r);

  Array1d a_l = mask.select(stencil[1].speed_of_sound, one);
  Array1d a_r = mask.select(stencil[2].speed_of_sound, one);

  Array<double, Dim> u_ll;
  Array<double, Dim> u_l;
  Array<double, Dim> u_r;
  Array<double, Dim> u_rr;
  for (int i = 0; i < Dim; ++i) {
    rho_u_l.row(i) = mask.select(rho_u_l.row(i), zero);
    rho_u_r.row(i) = mask.select(rho_u_r.row(i), zero);
    rho_u_ll.row(i) = left_mask.select(rho_u_ll.row(i), rho_u_l.row(i));
    rho_u_rr.row(i) = right_mask.select(rho_u_rr.row(i), rho_u_r.row(i));

    u_ll.row(i) = rho_u_ll.row(i) / rho_ll;
    u_l.row(i) = rho_u_l.row(i) / rho_l;
    u_r.row(i) = rho_u_r.row(i) / rho_r;
    u_rr.row(i) = rho_u_rr.row(i) / rho_rr;
  }

  int d = static_cast<int>(dir);
  Array1d lambda = Array1d::Constant(dt.count() / dx);

  Array1d drho_l = limiter(rho_ll, rho_l, rho_r);
  Array1d dp_l = limiter(p_ll, p_l, p_r);
  Array<double, Dim> du_l;
  for (int i = 0; i < Dim; ++i) {
    du_l.row(i) = limiter(u_ll.row(i), u_l.row(i), u_r.row(i));
  }

  Array1d drho_r = limiter(rho_l, rho_r, rho_rr);
  Array1d dp_r = limiter(p_l, p_r, p_rr);
  Array<double, Dim> du_r;
  for (int i = 0; i < Dim; ++i) {
    du_r.row(i) = limiter(u_l.row(i), u_r.row(i), u_rr.row(i));
  }

  Array1d rec_rho_l =
      rho_l +
      0.5 * (drho_l - lambda * (u_l.row(d) * drho_l + rho_l * du_l.row(d)));
  Array1d rec_p_l =
      p_l + 0.5 * (dp_l - lambda * (rho_l * a_l * a_l * du_l.row(d) +
                                    u_l.row(d) * dp_l));
  Array<double, Dim> rec_u_l;
  for (int i = 0; i < Dim; ++i) {
    rec_u_l.row(i) =
        u_l.row(i) + 0.5 * (du_l.row(i) - lambda * (u_l.row(d) * du_l.row(i)));
  }
  rec_u_l.row(d) -= 0.5 * lambda * dp_l / rho_l;

  Array1d rec_rho_r =
      rho_r -
      0.5 * (drho_r + lambda * (u_r.row(d) * drho_r + rho_r * du_r.row(d)));
  Array1d rec_p_r =
      p_r - 0.5 * (dp_r + lambda * (rho_r * a_r * a_r * du_r.row(d) +
                                    u_r.row(d) * dp_r));
  Array<double, Dim> rec_u_r;
  for (int i = 0; i < Dim; ++i) {
    rec_u_r.row(i) =
        u_r.row(i) - 0.5 * (du_r.row(i) + lambda * (u_r.row(d) * du_r.row(i)));
  }
  rec_u_r.row(d) -= 0.5 * lambda * dp_r / rho_r;

  auto squaredNorm = [](Array<double, Dim> u) -> Array1d {
    Array1d norm2 = u.matrix().colwise().squaredNorm();
    return norm2;
  };

  std::array<CompleteArray, 2> w;
  w[0].density = rec_rho_l;
  for (int i = 0; i < Dim; ++i) {
    w[0].momentum.row(i) = rec_rho_l * rec_u_l.row(i);
  }
  w[0].pressure = rec_p_l;
  w[0].speed_of_sound = (equation_.gamma_array_ * rec_p_l / rec_rho_l).sqrt();
  w[0].energy = 0.5 * rec_rho_l * squaredNorm(rec_u_l) +
                rec_p_l * equation_.gamma_minus_1_inv_array_;

  FUB_ASSERT((!mask || w[0].density > 0.0).all());
  FUB_ASSERT((!mask || w[0].pressure > 0.0).all());
  FUB_ASSERT((!mask || w[0].energy > 0.0).all());

  w[1].density = rec_rho_r;
  for (int i = 0; i < Dim; ++i) {
    w[1].momentum.row(i) = rec_rho_r * rec_u_r.row(i);
  }
  w[1].pressure = rec_p_r;
  w[1].speed_of_sound = (equation_.gamma_array_ * rec_p_r / rec_rho_r).sqrt();
  w[1].energy = 0.5 * rec_rho_r * squaredNorm(rec_u_r) +
                rec_p_r * equation_.gamma_minus_1_inv_array_;

  FUB_ASSERT((!mask || w[1].density > 0.0).all());
  FUB_ASSERT((!mask || w[1].pressure > 0.0).all());
  FUB_ASSERT((!mask || w[1].energy > 0.0).all());

  Hllem hllem{equation_};
  hllem.ComputeNumericFlux(flux, face_fractions, span{w},
                           volume_fractions.subspan<1, 2>(), dt, dx, dir);

  FUB_ASSERT((!isnan(flux.density)).all());
  FUB_ASSERT((!isnan(flux.momentum)).all());
  FUB_ASSERT((!isnan(flux.energy)).all());
}

} // namespace fub::perfect_gas