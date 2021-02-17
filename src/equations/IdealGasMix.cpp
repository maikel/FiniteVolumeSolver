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
  if (state.density <= 0)
    return;
  //  FUB_ASSERT(state.density > 0);
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

template <int Dim>
void IdealGasMix<Dim>::Flux(ConservativeArray& flux, const CompleteArray& state,
                            MaskArray mask, Direction dir) const noexcept {
  const int d0 = static_cast<int>(dir);
  const Array1d rho = mask.select(state.density, 1.0);
  const Array1d rho_u = mask.select(state.momentum.row(d0), 0.0);
  const Array1d velocity = rho_u / rho;
  const Array1d pressure = mask.select(state.pressure, 0.0);
  const Array1d energy = mask.select(state.energy, 0.0);
  flux.density = rho_u;
  for (int d = 0; d < Dim; ++d) {
    const Array1d rho_u_i = mask.select(state.momentum.row(d), 0.0);
    flux.momentum.row(d) = velocity * rho_u_i;
  }
  flux.momentum.row(d0) += pressure;
  flux.energy = velocity * (energy + pressure);
  for (int s = 0; s < flux.species.rows(); ++s) {
    const Array1d Yi = mask.select(state.species.row(s), 0.0);
    flux.species.row(s) = velocity * Yi;
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

auto ComputeSpeedOfSoundArray(const FlameMasterReactor& reactor,
                              MaskArray mask) {
  const Array1d cp = mask.select(reactor.GetCpArray(), 0.0);
  const Array1d gamma = cp / reactor.GetCvArray();
  const Array1d p = mask.select(reactor.GetPressureArray(), 0.0);
  const Array1d rho = mask.select(reactor.GetDensityArray(), 1.0);
  const Array1d a = (gamma * p / rho).sqrt();
  return std::make_tuple(a, cp, gamma);
}
} // namespace

template <int Dim>
Array<double, Dim, 1>
IdealGasMix<Dim>::Velocity(const ConservativeBase& cons) noexcept {
  Array<double, Dim, 1> velocity = cons.momentum / cons.density;
  return velocity;
}

template <int Dim>
Array<double, Dim>
IdealGasMix<Dim>::Velocity(const ConservativeArrayBase& cons) noexcept {
  Array<double, Dim> velocity;
  for (int d = 0; d < Dim; ++d) {
    velocity.row(d) = cons.momentum.row(d) / cons.density;
  }
  return velocity;
}

template <int Dim>
Array<double, Dim> IdealGasMix<Dim>::Velocity(const ConservativeArrayBase& cons,
                                              MaskArray mask) noexcept {
  Array<double, Dim> velocity;
  Array1d density = mask.select(cons.density, 1.0);
  for (int d = 0; d < Dim; ++d) {
    Array1d momentum = mask.select(cons.momentum.row(d), 0.0);
    velocity.row(d) = momentum / density;
  }
  return velocity;
}

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
  FUB_ASSERT((state.density > 0.0).all());
  FUB_ASSERT((state.pressure > 0.0).all());
  const ArrayXd& Y = reactor_.GetMassFractionsArray();
  for (int i = 0; i < Y.rows(); ++i) {
    state.species.row(i) = state.density * Y.row(i);
  }
  FUB_ASSERT((state.species.colwise().sum().isApprox(state.density)));
  state.temperature = reactor_.GetTemperatureArray();
  std::tie(state.speed_of_sound, state.c_p, state.gamma) =
      ComputeSpeedOfSoundArray(reactor_);
}

template <int Dim>
void IdealGasMix<Dim>::CompleteFromReactor(CompleteArray& state,
                                           const Array<double, Dim>& velocity,
                                           MaskArray mask) const {
  state.density = mask.select(reactor_.GetDensityArray(), 0.0);
  for (int i = 0; i < Dim; ++i) {
    state.momentum.row(i) = state.density * velocity.row(i);
  }
  const Array1d rhoE_internal =
      state.density * reactor_.GetInternalEnergyArray();
  const Array1d rhoE_kin = KineticEnergy(state.density, state.momentum, mask);
  state.energy = rhoE_internal + rhoE_kin;
  state.pressure = mask.select(reactor_.GetPressureArray(), 0.0);
  const ArrayXd& Y = reactor_.GetMassFractionsArray();
  for (int i = 0; i < Y.rows(); ++i) {
    state.species.row(i) = state.density * Y.row(i);
  }
  state.temperature = mask.select(reactor_.GetTemperatureArray(), 0.0);
  std::tie(state.speed_of_sound, state.c_p, state.gamma) =
      ComputeSpeedOfSoundArray(reactor_, mask);
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
  const Array<double, Dim, 1> velocity = Velocity(cons);
  CompleteFromReactor(complete, velocity);
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
  const Array<double, Dim> velocity = Velocity(cons);
  CompleteFromReactor(complete, velocity);
}

template <int Dim>
void IdealGasMix<Dim>::CompleteFromCons(CompleteArray& complete,
                                        const ConservativeArrayBase& cons,
                                        MaskArray mask) {
  const Array1d zero = Array1d::Zero();
  if (!mask.any()) {
    complete.density = zero;
    complete.momentum = Array<double, Dim>::Zero();
    complete.energy = zero;
    complete.species.fill(0.0);
    complete.pressure = zero;
    complete.temperature = zero;
    complete.speed_of_sound = zero;
    complete.c_p = zero;
    complete.gamma = zero;
    return;
  }
  Array1d rho_s = mask.select(cons.density, 1.0);
  FUB_ASSERT((rho_s > 0.0).all());
  reactor_.SetDensityArray(rho_s);
  reactor_.SetMassFractionsArray(cons.species, mask);
  reactor_.SetTemperatureArray(Array1d::Constant(300));
  Array1d rhoE = mask.select(cons.energy, 0.0);
  const Array1d rhoE_kin = KineticEnergy(cons.density, cons.momentum, mask);
  const Array1d e_internal = (rhoE - rhoE_kin) / rho_s;
  reactor_.SetInternalEnergyArray(e_internal, mask);
  const Array<double, Dim> velocity = Velocity(cons, mask);
  CompleteFromReactor(complete, velocity, mask);
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

} // namespace fub
