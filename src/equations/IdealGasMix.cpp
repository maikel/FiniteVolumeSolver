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

#include "fub/equations/IdealGasMix.hpp"

namespace fub {

template <int Dim>
void IdealGasMix<Dim>::Flux(Conservative& flux, const Complete& state,
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
  flux.species = velocity * state.species;
}

template <int Dim>
void IdealGasMix<Dim>::Flux(ConservativeArray& flux, const CompleteArray& state,
                            Direction dir) const noexcept {
  const int d0 = static_cast<int>(dir);
  const Array1d velocity = state.momentum.row(d0) / state.density;
  flux.density = state.momentum.row(d0);
  for (int d = 0; d < Dim; ++d) {
    flux.momentum.row(d) = velocity * state.momentum.row(d);
  }
  flux.momentum.row(d0) += state.pressure;
  flux.energy = velocity * (state.energy + state.pressure);
  for (int s = 0; s < flux.species.rows(); ++s) {
    flux.species.row(s) = velocity * state.species.row(s);
  }
}

namespace {
double ComputeSpeedOfSound(const FlameMasterReactor& reactor) {
  const double gamma = reactor.GetCp() / reactor.GetCv();
  FUB_ASSERT(1.0 < gamma && gamma <= 2.0);
  const double p = reactor.GetPressure();
  FUB_ASSERT(0.0 < p);
  const double rho = reactor.GetDensity();
  FUB_ASSERT(0.0 < rho);
  const double a = std::sqrt(gamma * p / rho);
  FUB_ASSERT(0.0 < a);
  return a;
}

auto ComputeSpeedOfSoundArray(const FlameMasterReactor& reactor) {
  const Array1d cp = reactor.GetCpArray();
  const Array1d gamma = cp / reactor.GetCvArray();
  const Array1d p = reactor.GetPressureArray();
  const Array1d rho = reactor.GetDensityArray();
  const Array1d a = (gamma * p / rho).sqrt();
  return std::make_tuple(a, cp, gamma);
}
} // namespace

template <int Dim>
void IdealGasMix<Dim>::SetReactorStateFromComplete(const Complete& state) {
  reactor_.SetDensity(state.density);
  reactor_.SetMassFractions(state.species);
  reactor_.SetTemperature(state.temperature);
}

template <int Dim>
void IdealGasMix<Dim>::CompleteFromReactor(
    Complete& state, const Eigen::Array<double, Dim, 1>& velocity) const {
  state.density = reactor_.GetDensity();
  state.momentum = state.density * velocity;
  const double rhoE_internal = state.density * reactor_.GetInternalEnergy();
  const double rhoE_kin = KineticEnergy(state.density, state.momentum);
  state.energy = rhoE_internal + rhoE_kin;
  state.pressure = reactor_.GetPressure();
  span<const double> Y = reactor_.GetMassFractions();
  std::transform(Y.begin(), Y.end(), &state.species[0],
                 [rho = state.density](double Yi) { return rho * Yi; });
  state.speed_of_sound = ComputeSpeedOfSound(reactor_);
  state.temperature = reactor_.GetTemperature();
  state.c_p = reactor_.GetCp();
  state.gamma = reactor_.GetCp() / reactor_.GetCv();
}

template <int Dim>
void IdealGasMix<Dim>::CompleteFromReactor(
    CompleteArray& state, const Array<double, Dim>& velocity) const {
  state.density = reactor_.GetDensityArray();
  for (int i = 0; i < Dim; ++i) {
    state.momentum.row(i) = state.density * velocity.row(i);
  }
  const Array1d rhoE_internal =
      state.density * reactor_.GetInternalEnergyArray();
  const Array1d rhoE_kin = KineticEnergy(state.density, state.momentum);
  state.energy = rhoE_internal + rhoE_kin;
  state.pressure = reactor_.GetPressureArray();
  const ArrayXd& Y = reactor_.GetMassFractionsArray();
  for (int i = 0; i < Y.rows(); ++i) {
    state.species.row(i) = state.density * Y.row(i);
  }
  state.temperature = reactor_.GetTemperatureArray();
  std::tie(state.speed_of_sound, state.c_p, state.gamma) =
      ComputeSpeedOfSoundArray(reactor_);
}

template <int Dim>
void IdealGasMix<Dim>::CompleteFromCons(Complete& complete,
                                        const ConservativeBase& cons) {
  reactor_.SetDensity(cons.density);
  reactor_.SetMassFractions(cons.species);
  const double rhoE_kin = KineticEnergy(cons.density, cons.momentum);
  const double e_internal = (cons.energy - rhoE_kin) / cons.density;
  reactor_.SetTemperature(300);
  reactor_.SetInternalEnergy(e_internal);
  FUB_ASSERT(reactor_.GetTemperature() > 0.0);
  FUB_ASSERT(cons.density == reactor_.GetDensity());
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  complete.species = cons.species;
  complete.pressure = reactor_.GetPressure();
  complete.speed_of_sound = ComputeSpeedOfSound(reactor_);
  complete.temperature = reactor_.GetTemperature();
  complete.c_p = reactor_.GetCp();
  complete.gamma = reactor_.GetCp() / reactor_.GetCv();
}

template <int Dim>
void IdealGasMix<Dim>::CompleteFromCons(CompleteArray& complete,
                                        const ConservativeArrayBase& cons) {
  reactor_.SetDensityArray(cons.density);
  reactor_.SetMassFractionsArray(cons.species);
  reactor_.SetTemperatureArray(Array1d::Constant(300));
  const Array1d rhoE_kin = KineticEnergy(cons.density, cons.momentum);
  const Array1d e_internal = (cons.energy - rhoE_kin) / cons.density;
  reactor_.SetInternalEnergyArray(e_internal);
  complete.density = cons.density;
  complete.momentum = cons.momentum;
  complete.energy = cons.energy;
  complete.species = cons.species;
  complete.pressure = reactor_.GetPressureArray();
  complete.temperature = reactor_.GetTemperatureArray();
  std::tie(complete.speed_of_sound, complete.c_p, complete.gamma) =
      ComputeSpeedOfSoundArray(reactor_);
}

template <int Dim>
void IdealGasMix<Dim>::CompleteFromCons(CompleteArray& complete,
                                        const ConservativeArrayBase& cons,
                                        MaskArray mask) {
  reactor_.SetDensityArray(cons.density);
  reactor_.SetMassFractionsArray(cons.species);
  reactor_.SetTemperatureArray(Array1d::Constant(300));
  const Array1d rhoE_kin = KineticEnergy(cons.density, cons.momentum);
  const Array1d e_internal = (cons.energy - rhoE_kin) / cons.density;
  reactor_.SetInternalEnergyArray(e_internal, mask);
  complete.density = mask.select(cons.density, Array1d::Zero());
  for (int i = 0; i < Dim; ++i) {
    complete.momentum.row(i) =
        mask.select(cons.momentum.row(i), Array1d::Zero());
  }
  complete.energy = mask.select(cons.energy, Array1d::Zero());
  for (int i = 0; i < cons.species.rows(); ++i) {
    complete.species.row(i) = mask.select(cons.species.row(i), Array1d::Zero());
  }
  complete.pressure = mask.select(reactor_.GetPressureArray(), Array1d::Zero());
  complete.temperature =
      mask.select(reactor_.GetTemperatureArray(), Array1d::Zero());
  std::tie(complete.speed_of_sound, complete.c_p, complete.gamma) =
      ComputeSpeedOfSoundArray(reactor_);
  complete.speed_of_sound =
      mask.select(complete.speed_of_sound, Array1d::Zero());
  complete.c_p = mask.select(complete.c_p, Array1d::Zero());
  complete.gamma = mask.select(complete.gamma, Array1d::Zero());
}

template class IdealGasMix<1>;
template class IdealGasMix<2>;
template class IdealGasMix<3>;

void Rotate(Conservative<IdealGasMix<2>>& rotated,
            const Conservative<IdealGasMix<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation,
            const IdealGasMix<2>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.species = state.species;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Rotate(Complete<IdealGasMix<2>>& rotated,
            const Complete<IdealGasMix<2>>& state,
            const Eigen::Matrix<double, 2, 2>& rotation,
            const IdealGasMix<2>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.pressure = state.pressure;
  rotated.speed_of_sound = state.speed_of_sound;
  rotated.temperature = state.temperature;
  rotated.c_p = state.c_p;
  rotated.gamma = state.gamma;
  rotated.species = state.species;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Rotate(Conservative<IdealGasMix<3>>& rotated,
            const Conservative<IdealGasMix<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation,
            const IdealGasMix<3>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.species = state.species;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Rotate(Complete<IdealGasMix<3>>& rotated,
            const Complete<IdealGasMix<3>>& state,
            const Eigen::Matrix<double, 3, 3>& rotation,
            const IdealGasMix<3>&) {
  rotated.density = state.density;
  rotated.energy = state.energy;
  rotated.pressure = state.pressure;
  rotated.speed_of_sound = state.speed_of_sound;
  rotated.temperature = state.temperature;
  rotated.c_p = state.c_p;
  rotated.gamma = state.gamma;
  rotated.species = state.species;
  rotated.momentum = (rotation * state.momentum.matrix()).array();
}

void Reflect(Complete<IdealGasMix<1>>& reflected,
             const Complete<IdealGasMix<1>>& state,
             const Eigen::Matrix<double, 1, 1>& normal, const IdealGasMix<1>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.species = state.species;
  reflected.temperature = state.temperature;
  reflected.c_p = state.c_p;
  reflected.gamma = state.gamma;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
}

void Reflect(Complete<IdealGasMix<2>>& reflected,
             const Complete<IdealGasMix<2>>& state,
             const Eigen::Vector2d& normal, const IdealGasMix<2>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.species = state.species;
  reflected.temperature = state.temperature;
  reflected.c_p = state.c_p;
  reflected.gamma = state.gamma;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
}

void Reflect(Complete<IdealGasMix<3>>& reflected,
             const Complete<IdealGasMix<3>>& state,
             const Eigen::Vector3d& normal, const IdealGasMix<3>&) {
  reflected.density = state.density;
  reflected.energy = state.energy;
  reflected.pressure = state.pressure;
  reflected.speed_of_sound = state.speed_of_sound;
  reflected.species = state.species;
  reflected.temperature = state.temperature;
  reflected.c_p = state.c_p;
  reflected.gamma = state.gamma;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
}

template <int Dim>
std::array<double, 2> EinfeldtSignalVelocities<IdealGasMix<Dim>>::
operator()(const IdealGasMix<Dim>&, const Complete& left, const Complete& right,
           Direction dir) const noexcept {
  FUB_ASSERT(left.density > 0.0);
  FUB_ASSERT(right.density > 0.0);
  const double rhoL = left.density;
  const double rhoR = right.density;
  const double rhoUL = left.momentum[int(dir)];
  const double rhoUR = right.momentum[int(dir)];
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
std::array<Array1d, 2> EinfeldtSignalVelocities<IdealGasMix<Dim>>::
operator()(const IdealGasMix<Dim>&, const CompleteArray& left,
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

template struct EinfeldtSignalVelocities<IdealGasMix<1>>;
template struct EinfeldtSignalVelocities<IdealGasMix<2>>;
template struct EinfeldtSignalVelocities<IdealGasMix<3>>;

namespace ideal_gas {
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
  for (int i = 0; i < Rank; ++i) {
    p.velocity.row(i) = q.momentum.row(i) / q.density;
  }
  p.temperature = q.temperature;
  for (int i = 0; i < p.mass_fractions.rows(); ++i) {
    p.mass_fractions.row(i) = q.species.row(i) / q.density;
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
  double dRdt = Rhat * u * sum_dXdx;
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
  const Array1d dRdt = Rhat * u * sum_dXdx;
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
void MinModLimiter(Primitive<Rank>& dpdx, const Primitive<Rank>& pL,
                   const Primitive<Rank>& pM, const Primitive<Rank>& pR) {
  auto MinMod = [](double& q, double qL, double qM, double qR) {
    const double sL = qM - qL;
    const double sR = qR - qM;
    if (sR > 0) {
      q = std::max(0.0, std::min(sL, sR));
    } else {
      q = std::min(0.0, std::max(sL, sR));
    }
  };
  MinMod(dpdx.pressure, pL.pressure, pM.pressure, pR.pressure);
  MinMod(dpdx.temperature, pL.temperature, pM.temperature, pR.temperature);
  for (int i = 0; i < Rank; ++i) {
    MinMod(dpdx.velocity[i], pL.velocity[i], pM.velocity[i], pR.velocity[i]);
  }
  for (int i = 0; i < dpdx.mass_fractions.rows(); ++i) {
    MinMod(dpdx.mass_fractions[i], pL.mass_fractions[i], pM.mass_fractions[i],
           pR.mass_fractions[i]);
  }
}

template <int Rank>
void MinModLimiter(PrimitiveArray<Rank>& dpdx, const PrimitiveArray<Rank>& pL,
                   const PrimitiveArray<Rank>& pM,
                   const PrimitiveArray<Rank>& pR) {
  auto MinMod = [](auto&& cons, Array1d qL, Array1d qM, Array1d qR) {
    const Array1d sL = qM - qL;
    const Array1d sR = qR - qM;
    const Array1d zero = Array1d::Zero();
    cons = (sR > 0).select(zero.max(sL.min(sR)), zero.min(sL.max(sR)));
  };
  MinMod(dpdx.pressure, pL.pressure, pM.pressure, pR.pressure);
  MinMod(dpdx.temperature, pL.temperature, pM.temperature, pR.temperature);
  for (int i = 0; i < Rank; ++i) {
    MinMod(dpdx.velocity.row(i), pL.velocity.row(i), pM.velocity.row(i),
           pR.velocity.row(i));
  }
  for (int i = 0; i < dpdx.mass_fractions.rows(); ++i) {
    MinMod(dpdx.mass_fractions.row(i), pL.mass_fractions.row(i),
           pM.mass_fractions.row(i), pR.mass_fractions.row(i));
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
  MinModLimiter(dpdx, pL, pM, pR);
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
  MinModLimiter(dpdx, pL, pM, pR);
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
  MinModLimiter(dpdx_array_, pL_array_, pM_array_, pR_array_);
  PrimDerivatives(GetEquation(), dpdt_array_, stencil[1], dpdx_array_, dir);
  pR_array_.pressure =
      pM_array_.pressure +
      0.5 * (dt_over_dx * dpdt_array_.pressure + dpdx_array_.pressure);
  for (int i = 0; i < Rank; ++i) {
    pR_array_.velocity.row(i) =
        pM_array_.velocity.row(i) +
        0.5 * (dt_over_dx * dpdt_array_.velocity.row(i) +
               dpdx_array_.velocity.row(i));
  }
  pR_array_.temperature =
      pM_array_.temperature +
      0.5 * (dt_over_dx * dpdt_array_.temperature + dpdx_array_.temperature);
  for (int i = 0; i < nspecies; ++i) {
    pR_array_.mass_fractions.row(i) =
        pM_array_.mass_fractions.row(i) +
        0.5 * (dt_over_dx * dpdt_array_.mass_fractions.row(i) +
               dpdx_array_.mass_fractions.row(i));
  }
  CompleteFromPrim(GetEquation(), stencil_array_[0], pR_array_);

  ToPrim(pL_array_, stencil[1]);
  ToPrim(pM_array_, stencil[2]);
  ToPrim(pR_array_, stencil[3]);
  MinModLimiter(dpdx_array_, pL_array_, pM_array_, pR_array_);
  PrimDerivatives(GetEquation(), dpdt_array_, stencil[2], dpdx_array_, dir);
  pL_array_.pressure =
      pM_array_.pressure +
      0.5 * (dt_over_dx * dpdt_array_.pressure - dpdx_array_.pressure);
  for (int i = 0; i < Rank; ++i) {
    pL_array_.velocity.row(i) =
        pM_array_.velocity.row(i) +
        0.5 * (dt_over_dx * dpdt_array_.velocity.row(i) -
               dpdx_array_.velocity.row(i));
  }
  pL_array_.temperature =
      pM_array_.temperature +
      0.5 * (dt_over_dx * dpdt_array_.temperature - dpdx_array_.temperature);
  for (int i = 0; i < nspecies; ++i) {
    pL_array_.mass_fractions.row(i) =
        pM_array_.mass_fractions.row(i) +
        0.5 * (dt_over_dx * dpdt_array_.mass_fractions.row(i) -
               dpdx_array_.mass_fractions.row(i));
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
  const int nspecies = GetEquation().GetReactor().GetNSpecies();
  const Array1d dt_over_dx = Array1d::Constant(dt.count() / dx);
  if ((volume_fractions[0] < 1.0).any() || (volume_fractions[3].any() < 1.0)) {
    hll_.ComputeNumericFlux(flux, face_fractions, stencil.template subspan<1, 2>(),
                                   volume_fractions.template subspan<1, 2>(), dt, dx,
                                   dir);
    return;
  }

  ToPrim(pL_array_, stencil[0]);
  ToPrim(pM_array_, stencil[1]);
  ToPrim(pR_array_, stencil[2]);
  MinModLimiter(dpdx_array_, pL_array_, pM_array_, pR_array_);
  PrimDerivatives(GetEquation(), dpdt_array_, stencil[1], dpdx_array_, dir);
  Array1d zero = Array1d::Zero();
  MaskArray mask0 = (volume_fractions[0] > 0.0);
  dpdx_array_.pressure = mask0.select(dpdx_array_.pressure, zero);
  dpdt_array_.pressure = mask0.select(dpdt_array_.pressure, zero);
  dpdx_array_.temperature = mask0.select(dpdx_array_.temperature, zero);
  dpdt_array_.temperature = mask0.select(dpdt_array_.temperature, zero);
  for (int i = 0; i < Rank; ++i) {
    dpdx_array_.velocity.row(i) =
        mask0.select(dpdx_array_.velocity.row(i), zero);
    dpdt_array_.velocity.row(i) =
        mask0.select(dpdt_array_.velocity.row(i), zero);
  }
  for (int i = 0; i < nspecies; ++i) {
    dpdx_array_.mass_fractions.row(i) =
        mask0.select(dpdx_array_.mass_fractions.row(i), zero);
    dpdt_array_.mass_fractions.row(i) =
        mask0.select(dpdt_array_.mass_fractions.row(i), zero);
  }
  pR_array_.pressure =
      pM_array_.pressure +
      0.5 * (dt_over_dx * dpdt_array_.pressure + dpdx_array_.pressure);
  for (int i = 0; i < Rank; ++i) {
    pR_array_.velocity.row(i) =
        pM_array_.velocity.row(i) +
        0.5 * (dt_over_dx * dpdt_array_.velocity.row(i) +
               dpdx_array_.velocity.row(i));
  }
  pR_array_.temperature =
      pM_array_.temperature +
      0.5 * (dt_over_dx * dpdt_array_.temperature + dpdx_array_.temperature);
  for (int i = 0; i < nspecies; ++i) {
    pR_array_.mass_fractions.row(i) =
        pM_array_.mass_fractions.row(i) +
        0.5 * (dt_over_dx * dpdt_array_.mass_fractions.row(i) +
               dpdx_array_.mass_fractions.row(i));
  }
  CompleteFromPrim(GetEquation(), stencil_array_[0], pR_array_);

  ToPrim(pL_array_, stencil[1]);
  ToPrim(pM_array_, stencil[2]);
  ToPrim(pR_array_, stencil[3]);
  MinModLimiter(dpdx_array_, pL_array_, pM_array_, pR_array_);
  PrimDerivatives(GetEquation(), dpdt_array_, stencil[2], dpdx_array_, dir);
  MaskArray mask3 = (volume_fractions[3] > 0.0);
  dpdx_array_.pressure = mask3.select(dpdx_array_.pressure, zero);
  dpdt_array_.pressure = mask3.select(dpdt_array_.pressure, zero);
  dpdx_array_.temperature = mask3.select(dpdx_array_.temperature, zero);
  dpdt_array_.temperature = mask3.select(dpdt_array_.temperature, zero);
  for (int i = 0; i < Rank; ++i) {
    dpdx_array_.velocity.row(i) =
        mask3.select(dpdx_array_.velocity.row(i), zero);
    dpdt_array_.velocity.row(i) =
        mask3.select(dpdt_array_.velocity.row(i), zero);
  }
  for (int i = 0; i < nspecies; ++i) {
    dpdx_array_.mass_fractions.row(i) =
        mask3.select(dpdx_array_.mass_fractions.row(i), zero);
    dpdt_array_.mass_fractions.row(i) =
        mask3.select(dpdt_array_.mass_fractions.row(i), zero);
  }
  pL_array_.pressure =
      pM_array_.pressure +
      0.5 * (dt_over_dx * dpdt_array_.pressure - dpdx_array_.pressure);
  for (int i = 0; i < Rank; ++i) {
    pL_array_.velocity.row(i) =
        pM_array_.velocity.row(i) +
        0.5 * (dt_over_dx * dpdt_array_.velocity.row(i) -
               dpdx_array_.velocity.row(i));
  }
  pL_array_.temperature =
      pM_array_.temperature +
      0.5 * (dt_over_dx * dpdt_array_.temperature - dpdx_array_.temperature);
  for (int i = 0; i < nspecies; ++i) {
    pL_array_.mass_fractions.row(i) =
        pM_array_.mass_fractions.row(i) +
        0.5 * (dt_over_dx * dpdt_array_.mass_fractions.row(i) -
               dpdx_array_.mass_fractions.row(i));
  }
  CompleteFromPrim(GetEquation(), stencil_array_[1], pL_array_);
  hll_.ComputeNumericFlux(flux, face_fractions, stencil_array_,
                          volume_fractions.template subspan<1, 2>(), dt, dx,
                          dir);
}

template <int Rank>
double MusclHancockPrimitive<Rank>::ComputeStableDt(
    span<const Complete, 4> states, double dx, Direction dir) noexcept {
  return 0.5 * hll_.ComputeStableDt(states.template subspan<1, 2>(), dx, dir);
}

template <int Rank>
Array1d MusclHancockPrimitive<Rank>::ComputeStableDt(
    span<const CompleteArray, 4> states, double dx, Direction dir) noexcept {
  return 0.5 * hll_.ComputeStableDt(states.template subspan<1, 2>(), dx, dir);
}

template <int Rank>
Array1d MusclHancockPrimitive<Rank>::ComputeStableDt(
    span<const CompleteArray, 4> states, Array1d face_fraction,
    span<const Array1d, 4> volume_fraction, double dx, Direction dir) noexcept {
  return 0.5 * hll_.ComputeStableDt(
                   states.template subspan<1, 2>(), face_fraction,
                   volume_fraction.template subspan<1, 2>(), dx, dir);
}

template class MusclHancockPrimitive<1>;
template class MusclHancockPrimitive<2>;
template class MusclHancockPrimitive<3>;
} // namespace ideal_gas

template class FluxMethod<ideal_gas::MusclHancockPrimitive<1>>;
template class FluxMethod<ideal_gas::MusclHancockPrimitive<2>>;
template class FluxMethod<ideal_gas::MusclHancockPrimitive<3>>;

} // namespace fub
