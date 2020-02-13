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

}
