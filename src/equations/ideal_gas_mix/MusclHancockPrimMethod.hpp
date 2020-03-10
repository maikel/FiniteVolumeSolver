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

#include "fub/equations/ideal_gas_mix/MusclHancockPrimMethod.hpp"

namespace fub::ideal_gas {

template <int Rank>
MusclHancockPrimitive<Rank>::MusclHancockPrimitive(const IdealGasMix<Rank>& eq)
    : hll_(eq, EinfeldtSignalVelocities<IdealGasMix<Rank>>{}) {
  const int nspecies = eq.GetReactor().GetNSpecies();
  dpdx.mass_fractions.resize(nspecies);
  dpdt.mass_fractions.resize(nspecies);
  pL.mass_fractions.resize(nspecies);
  pM.mass_fractions.resize(nspecies);
  pR.mass_fractions.resize(nspecies);
  dpdx_array_.mass_fractions.resize(nspecies, kDefaultChunkSize);
  dpdt_array_.mass_fractions.resize(nspecies, kDefaultChunkSize);
  pL_array_.mass_fractions.resize(nspecies, kDefaultChunkSize);
  pM_array_.mass_fractions.resize(nspecies, kDefaultChunkSize);
  pR_array_.mass_fractions.resize(nspecies, kDefaultChunkSize);
}

namespace {
template <int Rank>
void CompleteFromPrim(IdealGasMix<Rank>& eq,
                      Complete<IdealGasMix<Rank>>& complete,
                      const Primitive<Rank>& prim) {
  FlameMasterReactor& reactor = eq.GetReactor();
  reactor.SetDensity(1.0);
  reactor.SetMassFractions(prim.mass_fractions);
  reactor.SetTemperature(prim.temperature);
  reactor.SetPressure(prim.pressure);
  eq.CompleteFromReactor(complete, prim.velocity);
}

template <int Rank>
void CompleteFromPrim(IdealGasMix<Rank>& eq,
                      CompleteArray<IdealGasMix<Rank>>& complete,
                      const PrimitiveArray<Rank>& prim) {
  FlameMasterReactor& reactor = eq.GetReactor();
  reactor.SetDensityArray(Array1d::Constant(1.0));
  reactor.SetMassFractionsArray(prim.mass_fractions);
  reactor.SetTemperatureArray(prim.temperature);
  reactor.SetPressureArray(prim.pressure);
  eq.CompleteFromReactor(complete, prim.velocity);
}

template <int Rank>
void CompleteFromPrim(IdealGasMix<Rank>& eq,
                      CompleteArray<IdealGasMix<Rank>>& complete,
                      const PrimitiveArray<Rank>& prim, MaskArray mask) {
  FlameMasterReactor& reactor = eq.GetReactor();
  reactor.SetDensityArray(Array1d::Constant(1.0));
  reactor.SetMassFractionsArray(prim.mass_fractions);
  reactor.SetTemperatureArray(prim.temperature);
  reactor.SetPressureArray(prim.pressure);
  eq.CompleteFromReactor(complete, prim.velocity, mask);
}

template <int Rank>
void ToPrim(Primitive<Rank>& p, const Complete<IdealGasMix<Rank>>& q) {
  p.pressure = q.pressure;
  p.velocity = q.momentum / q.density;
  p.temperature = q.temperature;
  std::transform(q.species.data(), q.species.data() + q.species.size(),
                 p.mass_fractions.data(),
                 [rho = q.density](double rhoY) { return rhoY / rho; });
}

template <int Rank>
void ToPrim(PrimitiveArray<Rank>& p,
            const CompleteArray<IdealGasMix<Rank>>& q) {
  p.pressure = q.pressure;
  Array1d rho_s = (q.density > 0.0).select(q.density, Array1d::Constant(1.0));
  for (int i = 0; i < Rank; ++i) {
    p.velocity.row(i) = q.momentum.row(i) / rho_s;
  }
  p.temperature = q.temperature;
  for (int i = 0; i < p.mass_fractions.rows(); ++i) {
    p.mass_fractions.row(i) = q.species.row(i) / rho_s;
  }
}

template <int Rank>
void ToPrim(PrimitiveArray<Rank>& p, const CompleteArray<IdealGasMix<Rank>>& q,
            MaskArray mask) {
  const Array1d zero = Array1d::Zero();
  p.pressure = mask.select(q.pressure, zero);
  Array1d rho_s = (q.density > 0.0).select(q.density, Array1d::Constant(1.0));
  for (int i = 0; i < Rank; ++i) {
    p.velocity.row(i) = mask.select(q.momentum.row(i) / rho_s, zero);
  }
  p.temperature = mask.select(q.temperature, zero);
  for (int i = 0; i < p.mass_fractions.rows(); ++i) {
    p.mass_fractions.row(i) = mask.select(q.species.row(i) / rho_s, zero);
  }
}

template <int Rank>
void PrimDerivatives(IdealGasMix<Rank>& eq, Primitive<Rank>& dt,
                     const Complete<IdealGasMix<Rank>>& q,
                     const Primitive<Rank>& dx, Direction dir) {
  FlameMasterReactor& reactor = eq.GetReactor();
  const double dTdx = dx.temperature;
  const double dpdx = dx.pressure;
  const int d = static_cast<int>(dir);
  const double dudx = dx.velocity[d];
  const auto& dYdx = dx.mass_fractions;

  const double p = q.pressure;
  const double rho = q.density;
  const double T = q.temperature;
  const double u = q.momentum[d] / q.density;
  const double Rhat = reactor.GetUniversalGasConstant();

  double c_v = q.c_p / q.gamma;
  double R = q.c_p - c_v;
  span<const double> M = reactor.GetMolarMasses();
  double sum_dXdx = dYdx[0] / M[0];
  for (int i = 1; i < dYdx.rows(); ++i) {
    sum_dXdx += dYdx[i] / M[i];
  }
  double dRdt = -Rhat * u * sum_dXdx;
  double dTdt = (-p * dudx - rho * u * c_v * dTdx) / (rho * c_v);
  double drhodx =
      dpdx / (R * T) - p / (R * T * T) * dTdx - rho * Rhat * sum_dXdx / R;
  double dpdt =
      -R * T * (u * drhodx + rho * dudx) + rho * R * dTdt + rho * T * dRdt;
  dt.velocity = -u * dx.velocity;
  dt.velocity[d] -= dpdx / rho;
  dt.mass_fractions = -u * dYdx;
  dt.pressure = dpdt;
  dt.temperature = dTdt;
}

template <int Rank>
void PrimDerivatives(IdealGasMix<Rank>& eq, PrimitiveArray<Rank>& dt,
                     const CompleteArray<IdealGasMix<Rank>>& q,
                     const PrimitiveArray<Rank>& dx, Direction dir) {
  FlameMasterReactor& reactor = eq.GetReactor();
  const Array1d dTdx = dx.temperature;
  const Array1d dpdx = dx.pressure;
  const int d = static_cast<int>(dir);
  const Array1d dudx = dx.velocity.row(d);
  const auto& dYdx = dx.mass_fractions;

  const Array1d p = q.pressure;
  const Array1d rho = q.density;
  const Array1d T = q.temperature;
  const Array1d u = q.momentum.row(d) / q.density;
  const Array1d Rhat = Array1d::Constant(reactor.GetUniversalGasConstant());

  const Array1d c_v = q.c_p / q.gamma;
  const Array1d R = q.c_p - c_v;
  span<const double> M = reactor.GetMolarMasses();
  Array1d sum_dXdx = dYdx.row(0) / M[0];
  for (int i = 1; i < dYdx.rows(); ++i) {
    sum_dXdx += dYdx.row(i) / M[i];
  }
  const Array1d dRdt = -Rhat * u * sum_dXdx;
  const Array1d dTdt = (-p * dudx - rho * u * c_v * dTdx) / (rho * c_v);
  const Array1d drhodx =
      dpdx / (R * T) - p / (R * T * T) * dTdx - rho * Rhat * sum_dXdx / R;
  const Array1d dpdt =
      -R * T * (u * drhodx + rho * dudx) + rho * R * dTdt + rho * T * dRdt;
  for (int i = 0; i < Rank; ++i) {
    dt.velocity.row(i) = -u * dx.velocity.row(i);
  }
  dt.velocity.row(d) -= dpdx / rho;
  for (int i = 0; i < dYdx.rows(); ++i) {
    dt.mass_fractions.row(i) = -u * dYdx.row(i);
  }
  dt.pressure = dpdt;
  dt.temperature = dTdt;
}

template <int Rank>
void PrimDerivatives(IdealGasMix<Rank>& eq, PrimitiveArray<Rank>& dt,
                     const CompleteArray<IdealGasMix<Rank>>& q,
                     const PrimitiveArray<Rank>& dx, Direction dir,
                     MaskArray mask) {
  const Array1d zero = Array1d::Constant(0.0);
  const Array1d one = Array1d::Constant(1.0);

  FlameMasterReactor& reactor = eq.GetReactor();
  const Array1d dTdx = dx.temperature;
  const Array1d dpdx = dx.pressure;
  const int d = static_cast<int>(dir);
  const Array1d dudx = dx.velocity.row(d);
  const auto& dYdx = dx.mass_fractions;

  const Array1d p = mask.select(q.pressure, zero);
  const Array1d rho = mask.select(q.density, zero);
  const Array1d rho_s = mask.select(q.density, one);
  const Array1d T = mask.select(q.temperature, zero);
  const Array1d T_s = mask.select(q.temperature, one);
  const Array1d u = mask.select(q.momentum.row(d) / rho_s, zero);
  const Array1d Rhat = Array1d::Constant(reactor.GetUniversalGasConstant());

  const Array1d c_v = mask.select(q.c_p / q.gamma, one);
  const Array1d R = mask.select(q.c_p - c_v, zero);
  const Array1d R_s = mask.select(q.c_p - c_v, one);
  span<const double> M = reactor.GetMolarMasses();
  Array1d sum_dXdx = dYdx.row(0) / M[0];
  for (int i = 1; i < dYdx.rows(); ++i) {
    sum_dXdx += dYdx.row(i) / M[i];
  }
  const Array1d dRdt = -u * Rhat * sum_dXdx;
  const Array1d dTdt = (-p * dudx - rho * u * c_v * dTdx) / (rho_s * c_v);
  const Array1d drhodx = dpdx / (R_s * T_s) - p / (R_s * T_s * T_s) * dTdx -
                         rho * Rhat * sum_dXdx / R_s;
  const Array1d dpdt =
      -R * T * (u * drhodx + rho * dudx) + rho * R * dTdt + rho * T * dRdt;
  for (int i = 0; i < Rank; ++i) {
    dt.velocity.row(i) = -u * dx.velocity.row(i);
  }
  dt.velocity.row(d) -= dpdx / rho_s;
  for (int i = 0; i < dYdx.rows(); ++i) {
    dt.mass_fractions.row(i) = -u * dYdx.row(i);
  }
  dt.pressure = dpdt;
  dt.temperature = dTdt;
}

template <int Rank>
void ApplyLimiter(Primitive<Rank>& dpdx, const Primitive<Rank>& pL,
                  const Primitive<Rank>& pM, const Primitive<Rank>& pR) {
  auto Limiter = [](double& q, double qL, double qM, double qR) {
    const double sL = qM - qL;
    const double sR = qR - qM;
    if (sR > 0) {
      q = std::max(0.0, std::min(sL, sR));
    } else {
      q = std::min(0.0, std::max(sL, sR));
    }
  };
  Limiter(dpdx.pressure, pL.pressure, pM.pressure, pR.pressure);
  Limiter(dpdx.temperature, pL.temperature, pM.temperature, pR.temperature);
  for (int i = 0; i < Rank; ++i) {
    Limiter(dpdx.velocity[i], pL.velocity[i], pM.velocity[i], pR.velocity[i]);
  }
  for (int i = 0; i < dpdx.mass_fractions.rows(); ++i) {
    Limiter(dpdx.mass_fractions[i], pL.mass_fractions[i], pM.mass_fractions[i],
            pR.mass_fractions[i]);
  }
}

template <int Rank>
void ApplyLimiter(PrimitiveArray<Rank>& dpdx, const PrimitiveArray<Rank>& pL,
                  const PrimitiveArray<Rank>& pM,
                  const PrimitiveArray<Rank>& pR) {
  auto dqdx = [](Array1d qL, Array1d /*qM*/, Array1d qR) -> Array1d {
    Array1d delta_q = 0.5 * (qR - qL);
    return delta_q;
  };
  auto Limiter = [](Array1d qL, Array1d qM, Array1d qR) -> Array1d {
    Array1d delta_q = 0.5 * (qR - qL);
    Array1d zeros = Array1d::Zero();
    Array1d ones = Array1d::Constant(1.0);
    MaskArray is_relevant = delta_q.abs() > 1e-12;
    Array1d delta_q_if_relevant = is_relevant.select(delta_q, ones);
    Array1d rL = is_relevant.select((qM - qL) / delta_q_if_relevant, zeros);
    Array1d rR = is_relevant.select((qR - qM) / delta_q_if_relevant, zeros);
    MaskArray is_positive = rL > 0.0 && rR > 0.0;
    Array1d r = is_positive.select(rL.min(rR), zeros);
    Array1d r_1 = 2 * r / (Array1d::Constant(1.0) + r);
    Array1d r_2 = 0.5 * (r + Array1d::Constant(1.0));
    MaskArray less_than_one = r < 1.0;
    Array1d sigma = less_than_one.select(r_1, r_2);
    return sigma;
  };
  dpdx.pressure = Limiter(pL.pressure, pM.pressure, pR.pressure) *
                  dqdx(pL.pressure, pM.pressure, pR.pressure);
  dpdx.temperature = Limiter(pL.temperature, pM.temperature, pR.temperature) *
                     dqdx(pL.temperature, pM.temperature, pR.temperature);
  for (int i = 0; i < Rank; ++i) {
    dpdx.velocity.row(i) =
        Limiter(pL.velocity.row(i), pM.velocity.row(i), pR.velocity.row(i)) *
        dqdx(pL.velocity.row(i), pM.velocity.row(i), pR.velocity.row(i));
  }
  Array1d limiter_rhoY =
      Limiter(pL.mass_fractions.row(0), pM.mass_fractions.row(0),
              pR.mass_fractions.row(0));
  for (int i = 0; i < dpdx.mass_fractions.rows(); ++i) {
    dpdx.mass_fractions.row(i) =
        limiter_rhoY * dqdx(pL.mass_fractions.row(i), pM.mass_fractions.row(i),
                            pR.mass_fractions.row(i));
  }
}

template <int Rank>
void ApplyLimiter(PrimitiveArray<Rank>& dpdx, const PrimitiveArray<Rank>& pL,
                  const PrimitiveArray<Rank>& pM,
                  const PrimitiveArray<Rank>& pR, MaskArray mask) {
  auto dqdx = [](Array1d qL, Array1d /*qM*/, Array1d qR, MaskArray mask) -> Array1d {
    Array1d delta_q = 0.5 * (qR - qL);
    Array1d masked_delta_q = mask.select(delta_q, 0.0);
    return masked_delta_q;
  };
  auto Limiter = [](Array1d qL, Array1d qM, Array1d qR, MaskArray mask) -> Array1d {
    Array1d delta_q = mask.select(0.5 * (qR - qL), 0.0);
    Array1d zeros = Array1d::Zero();
    Array1d ones = Array1d::Constant(1.0);
    MaskArray is_relevant = delta_q.abs() > 1e-12;
    Array1d delta_q_if_relevant = is_relevant.select(delta_q, ones);
    Array1d rL = is_relevant.select((qM - qL) / delta_q_if_relevant, zeros);
    Array1d rR = is_relevant.select((qR - qM) / delta_q_if_relevant, zeros);
    MaskArray is_positive = is_relevant && rL > 0.0 && rR > 0.0;
    Array1d r = is_positive.select(rL.min(rR), zeros);
    Array1d r_1 = 2 * r / (Array1d::Constant(1.0) + r);
    Array1d r_2 = 0.5 * (r + Array1d::Constant(1.0));
    MaskArray less_than_one = r < 1.0;
    Array1d sigma = less_than_one.select(r_1, r_2);
    return sigma;
  };
  dpdx.pressure = Limiter(pL.pressure, pM.pressure, pR.pressure, mask) *
                  dqdx(pL.pressure, pM.pressure, pR.pressure, mask);
  dpdx.temperature = Limiter(pL.temperature, pM.temperature, pR.temperature, mask) *
                     dqdx(pL.temperature, pM.temperature, pR.temperature, mask);
  for (int i = 0; i < Rank; ++i) {
    dpdx.velocity.row(i) =
        Limiter(pL.velocity.row(i), pM.velocity.row(i), pR.velocity.row(i), mask) *
        dqdx(pL.velocity.row(i), pM.velocity.row(i), pR.velocity.row(i), mask);
  }
  Array1d limiter_Y =
      Limiter(pL.mass_fractions.row(0), pM.mass_fractions.row(0),
              pR.mass_fractions.row(0), mask);
  for (int i = 0; i < dpdx.mass_fractions.rows(); ++i) {
    dpdx.mass_fractions.row(i) =
        limiter_Y * dqdx(pL.mass_fractions.row(i), pM.mass_fractions.row(i),
                            pR.mass_fractions.row(i), mask);
  }
}

} // namespace

template <int Rank>
void MusclHancockPrimitive<Rank>::ComputeNumericFlux(
    Conservative& flux, span<const Complete, 4> stencil, Duration dt, double dx,
    Direction dir) {
  const double dt_over_dx = dt.count() / dx;

  ToPrim(pL, stencil[0]);
  ToPrim(pM, stencil[1]);
  ToPrim(pR, stencil[2]);
  ApplyLimiter(dpdx, pL, pM, pR);
  PrimDerivatives(GetEquation(), dpdt, stencil[1], dpdx, dir);

  pR.pressure =
      pM.pressure + 0.5 * (dt_over_dx * dpdt.pressure + dpdx.pressure);
  pR.velocity =
      pM.velocity + 0.5 * (dt_over_dx * dpdt.velocity + dpdx.velocity);
  pR.temperature =
      pM.temperature + 0.5 * (dt_over_dx * dpdt.temperature + dpdx.temperature);
  pR.mass_fractions =
      pM.mass_fractions +
      0.5 * (dt_over_dx * dpdt.mass_fractions + dpdx.mass_fractions);
  CompleteFromPrim(GetEquation(), stencil_[0], pR);

  ToPrim(pL, stencil[1]);
  ToPrim(pM, stencil[2]);
  ToPrim(pR, stencil[3]);
  ApplyLimiter(dpdx, pL, pM, pR);
  PrimDerivatives(GetEquation(), dpdt, stencil[2], dpdx, dir);
  pL.pressure =
      pM.pressure + 0.5 * (dt_over_dx * dpdt.pressure - dpdx.pressure);
  pL.velocity =
      pM.velocity + 0.5 * (dt_over_dx * dpdt.velocity - dpdx.velocity);
  pL.temperature =
      pM.temperature + 0.5 * (dt_over_dx * dpdt.temperature - dpdx.temperature);
  pL.mass_fractions =
      pM.mass_fractions +
      0.5 * (dt_over_dx * dpdt.mass_fractions - dpdx.mass_fractions);

  CompleteFromPrim(GetEquation(), stencil_[1], pL);
  hll_.ComputeNumericFlux(flux, stencil_, dt, dx, dir);
}

template <int Rank>
void MusclHancockPrimitive<Rank>::ComputeNumericFlux(
    ConservativeArray& flux, span<const CompleteArray, 4> stencil, Duration dt,
    double dx, Direction dir) {
  const int nspecies = GetEquation().GetReactor().GetNSpecies();
  const Array1d dt_over_dx = Array1d::Constant(dt.count() / dx);

  ToPrim(pL_array_, stencil[0]);
  ToPrim(pM_array_, stencil[1]);
  ToPrim(pR_array_, stencil[2]);
  ApplyLimiter(dpdx_array_, pL_array_, pM_array_, pR_array_);
  PrimDerivatives(GetEquation(), dpdt_array_, stencil[1], dpdx_array_, dir);
  pR_array_.pressure =
      pM_array_.pressure +
      0.5 * (dpdx_array_.pressure + dt_over_dx * dpdt_array_.pressure);
  for (int i = 0; i < Rank; ++i) {
    pR_array_.velocity.row(i) =
        pM_array_.velocity.row(i) +
        0.5 * (dpdx_array_.velocity.row(i) +
               dt_over_dx * dpdt_array_.velocity.row(i));
  }
  pR_array_.temperature =
      pM_array_.temperature +
      0.5 * (dpdx_array_.temperature + dt_over_dx * dpdt_array_.temperature);
  for (int i = 0; i < nspecies; ++i) {
    pR_array_.mass_fractions.row(i) =
        pM_array_.mass_fractions.row(i) +
        0.5 * (dpdx_array_.mass_fractions.row(i) +
               dt_over_dx * dpdt_array_.mass_fractions.row(i));
  }
  CompleteFromPrim(GetEquation(), stencil_array_[0], pR_array_);

  ToPrim(pL_array_, stencil[1]);
  ToPrim(pM_array_, stencil[2]);
  ToPrim(pR_array_, stencil[3]);
  ApplyLimiter(dpdx_array_, pL_array_, pM_array_, pR_array_);
  PrimDerivatives(GetEquation(), dpdt_array_, stencil[2], dpdx_array_, dir);
  pL_array_.pressure =
      pM_array_.pressure +
      0.5 * (-dpdx_array_.pressure + dt_over_dx * dpdt_array_.pressure);
  for (int i = 0; i < Rank; ++i) {
    pL_array_.velocity.row(i) =
        pM_array_.velocity.row(i) +
        0.5 * (-dpdx_array_.velocity.row(i) +
               dt_over_dx * dpdt_array_.velocity.row(i));
  }
  pL_array_.temperature =
      pM_array_.temperature -
      0.5 * (-dpdx_array_.temperature + dt_over_dx * dpdt_array_.temperature);
  for (int i = 0; i < nspecies; ++i) {
    pL_array_.mass_fractions.row(i) =
        pM_array_.mass_fractions.row(i) +
        0.5 * (-dpdx_array_.mass_fractions.row(i) +
               dt_over_dx * dpdt_array_.mass_fractions.row(i));
  }
  CompleteFromPrim(GetEquation(), stencil_array_[1], pL_array_);
  hll_.ComputeNumericFlux(flux, stencil_array_, dt, dx, dir);
}

template <int Rank>
void MusclHancockPrimitive<Rank>::ComputeNumericFlux(
    ConservativeArray& flux, Array1d face_fractions,
    span<const CompleteArray, 4> stencil,
    span<const Array1d, 4> volume_fractions, Duration dt, double dx,
    Direction dir) {
  MaskArray mask = face_fractions > 0.0;
  if (!mask.any()) {
    return;
  }

  const int nspecies = GetEquation().GetReactor().GetNSpecies();
  const Array1d dt_over_dx = Array1d::Constant(dt.count() / dx);

  MaskArray left_mask = volume_fractions[0] > 0.0;
  MaskArray right_mask = volume_fractions[3] > 0.0;

  ToPrim(pL_array_, stencil[0], left_mask);
  ToPrim(pM_array_, stencil[1], mask);
  ToPrim(pR_array_, stencil[2], mask);
  ApplyLimiter(dpdx_array_, pL_array_, pM_array_, pR_array_, left_mask && mask);
  PrimDerivatives(GetEquation(), dpdt_array_, stencil[1], dpdx_array_, dir,
                  mask);

  pR_array_.pressure =
      pM_array_.pressure +
      0.5 * (dpdx_array_.pressure + dt_over_dx * dpdt_array_.pressure);
  for (int i = 0; i < Rank; ++i) {
    pR_array_.velocity.row(i) =
        pM_array_.velocity.row(i) +
        0.5 * (dpdx_array_.velocity.row(i) +
               dt_over_dx * dpdt_array_.velocity.row(i));
  }
  pR_array_.temperature =
      pM_array_.temperature +
      0.5 * (dpdx_array_.temperature + dt_over_dx * dpdt_array_.temperature);
  for (int i = 0; i < nspecies; ++i) {
    pR_array_.mass_fractions.row(i) =
        pM_array_.mass_fractions.row(i) +
        0.5 * (dpdx_array_.mass_fractions.row(i) +
               dt_over_dx * dpdt_array_.mass_fractions.row(i));
  }
  CompleteFromPrim(GetEquation(), stencil_array_[0], pR_array_, mask);

  ToPrim(pL_array_, stencil[1], mask);
  ToPrim(pM_array_, stencil[2], mask);
  ToPrim(pR_array_, stencil[3], right_mask);
  ApplyLimiter(dpdx_array_, pL_array_, pM_array_, pR_array_, mask && right_mask);
  PrimDerivatives(GetEquation(), dpdt_array_, stencil[2], dpdx_array_, dir,
                  mask);

  pL_array_.pressure =
      pM_array_.pressure +
      0.5 * (-dpdx_array_.pressure + dt_over_dx * dpdt_array_.pressure);
  for (int i = 0; i < Rank; ++i) {
    pL_array_.velocity.row(i) =
        pM_array_.velocity.row(i) +
        0.5 * (-dpdx_array_.velocity.row(i) +
               dt_over_dx * dpdt_array_.velocity.row(i));
  }
  pL_array_.temperature =
      pM_array_.temperature -
      0.5 * (-dpdx_array_.temperature + dt_over_dx * dpdt_array_.temperature);
  for (int i = 0; i < nspecies; ++i) {
    pL_array_.mass_fractions.row(i) =
        pM_array_.mass_fractions.row(i) +
        0.5 * (-dpdx_array_.mass_fractions.row(i) +
               dt_over_dx * dpdt_array_.mass_fractions.row(i));
  }
  CompleteFromPrim(GetEquation(), stencil_array_[1], pL_array_, mask);

  hll_.ComputeNumericFlux(flux, face_fractions, stencil_array_,
                          volume_fractions.template subspan<1, 2>(), dt, dx,
                          dir);
}

template <int Rank>
double MusclHancockPrimitive<Rank>::ComputeStableDt(
    span<const Complete, 4> states, double dx, Direction dir) noexcept {
  return hll_.ComputeStableDt(states.template subspan<1, 2>(), dx, dir);
}

template <int Rank>
Array1d MusclHancockPrimitive<Rank>::ComputeStableDt(
    span<const CompleteArray, 4> states, double dx, Direction dir) noexcept {
  return hll_.ComputeStableDt(states.template subspan<1, 2>(), dx, dir);
}

template <int Rank>
Array1d MusclHancockPrimitive<Rank>::ComputeStableDt(
    span<const CompleteArray, 4> states, Array1d face_fraction,
    span<const Array1d, 4> volume_fraction, double dx, Direction dir) noexcept {
  return hll_.ComputeStableDt(states.template subspan<1, 2>(), face_fraction,
                              volume_fraction.template subspan<1, 2>(), dx,
                              dir);
}

} // namespace fub::ideal_gas
