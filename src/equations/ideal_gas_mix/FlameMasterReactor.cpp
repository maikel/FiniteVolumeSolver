// Copyright (c) 2016 Philip Berndt
// Copyright (c) 2018 Maikel Nadolski
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

#include "fub/equations/ideal_gas_mix/FlameMasterReactor.hpp"
#include "fub/ode_solver/RadauSolver.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

namespace fub {
namespace {
void UpdateMassFractionsFromMoles(FlameMasterState& state) noexcept {
  for (std::size_t i = 0; i < static_cast<std::size_t>(state.nSpecies); i++) {
    const std::ptrdiff_t j = static_cast<std::ptrdiff_t>(i);
    state.massFractions[i] =
        state.moles[j] * state.molarMasses[i] / state.density;
  }
}

void UpdateMolesFromMassFractions(FlameMasterState& state) noexcept {
  for (std::size_t i = 0; i < static_cast<std::size_t>(state.nSpecies); i++) {
    const std::ptrdiff_t j = static_cast<std::ptrdiff_t>(i);
    state.moles[j] =
        state.massFractions[i] * state.density / state.molarMasses[i];
  }
}
//
// void UpdateMassFractionsFromMoles(FlameMasterArrayState& state,
//                                  const FlameMasterState& x) noexcept {
//  for (int i = 0; i < state.massFractions.rows(); i++) {
//    state.massFractions.row(i) =
//        state.moles.row(i) * x.molarMasses[i] / state.density;
//  }
//}

void UpdateMolesFromMassFractions(FlameMasterArrayState& state,
                                  const FlameMasterState& x) noexcept {
  for (int i = 0; i < state.moles.rows(); i++) {
    state.moles.row(i) =
        state.massFractions.row(i) * state.density / x.molarMasses[i];
    state.molesStorage.row(i + 1) = state.moles.row(i);
  }
}

/// Update function for the thermodynamic state
void UpdateThermoState(const FlameMasterMechanism& mechanism,
                       FlameMasterState& state, double temperature) {
  if (state.thermoTemp == temperature) {
    return;
  }
  span<double> h = state.enthalpies;
  span<double> cp = state.heat_capacities_at_constant_pressure;
  span<double> s = state.entropies;
  mechanism.ComputeThermoData(h, cp, temperature, s);
  state.thermoTemp = temperature;
}

void UpdateThermoState(const FlameMasterMechanism& mechanism,
                       FlameMasterArrayState& state, Array1d temperature) {
  if ((state.thermoTemp == temperature).all()) {
    return;
  }
  ArrayXd& h = state.enthalpies;
  ArrayXd& cp = state.heat_capacities_at_constant_pressure;
  mechanism.ComputeThermoData(h, cp, temperature);
  state.thermoTemp = temperature;
}

void UpdateThermoState(const FlameMasterMechanism& mechanism,
                       FlameMasterState& state) {
  UpdateThermoState(mechanism, state, *state.temperature);
}

int OdeRhsAdvance(double /* t */, span<const double> XD, span<double> Xdot,
                  const FlameMasterMechanism& mechanism,
                  FlameMasterState& state) {
  span<double> h = state.enthalpies;
  span<double> cp = state.heat_capacities_at_constant_pressure;
  double temperature = XD[0];

  // We need to work in a local variable rather than in XD, because
  // the QSS code might alter XD beyond its size (nSpecies), assuming it
  // had space for nSpeciesEffective
  FUB_ASSERT(XD.size() > 0);
  std::size_t size = static_cast<std::size_t>(XD.size() - 1);
  span<double> concentrations((double*)alloca(sizeof(double) * size),
                              static_cast<std::ptrdiff_t>(size));
  std::copy(XD.begin() + 1, XD.end(), concentrations.begin());

  // Compute the pressure from the ideal gas law
  double concentrationSum = 0;
  for (int i = 1; i < 1 + state.nSpecies; i++) {
    concentrationSum += XD[i];
  }
  double pressure = concentrationSum * state.gasConstant * temperature;

  // Compute species' derivatives
  mechanism.ComputeProductionRates(Xdot.subspan(1),                 // cdot
                                   state.reaction_rates,            // w
                                   state.rate_coefficients,         // k
                                   concentrations,                  // c
                                   state.third_body_concentrations, // M
                                   temperature, // Temperature
                                   pressure     // Pressure
  );

  // Compute temperature derivative
  //
  // The idea is that
  //
  //   dU/dt         = 0
  // ↔ d/dt ∫ c_v dT = 0
  // ↔ -c_v dT/dt    = ∫ d/dt c_v dT
  // ↔ -c_v dT/dt    = ∂/∂t U(T, t)
  // ↔ -c_v dT/dt    = \sum u_i \dot Y_i
  //
  // So this derivative implicitly solves the algebraic constraint
  //  U = const

  // Compute thermodynamic data
  // if (XD.data() != state.molesStorage.data()) {
  //   std::transform(XD.begin() + 1, XD.end(), state.moles.begin(),
  //                  [=](double ci) { return ci / concentrationSum; });
  //   *state.temperature = XD[0];
  // }
  UpdateThermoState(mechanism, state, temperature);
  double& tempDt = Xdot[0];
  tempDt = 0;
  double meanCp = 0;
  double meanMolarMass = 0;

  const double R = state.gasConstant;
  span<const double> M{state.molarMasses};
  const double rho = state.density;

  for (int i = 0; i < state.nSpecies; i++) {
    tempDt += (h[i] - R / M[i] * temperature) * (Xdot[i + 1] * M[i]);
    meanCp += cp[i] * XD[i + 1] * M[i] / rho;
    meanMolarMass += XD[i + 1] * M[i] / concentrationSum;
  }
  tempDt /= rho * (R / meanMolarMass - meanCp);

  return 0;
}

int OdeRhsSetP(span<double> dydp, span<const double> y, double p,
               FlameMasterReactor& reactor) {
  // We integrate dH/dp = V <=> dT = 1/cp 1/rho dp from p0 to the desired value
  // of p
  // Reset to stored temperature
  reactor.SetTemperature(y[0]);

  // Calculate derivative
  reactor.SetPressure(p);
  double cp_rho = reactor.GetDensity() * reactor.GetCp();
  dydp[0] = 1. / cp_rho;

  // Also integrate velocity du = -sqrt(1/(gamma p rho)) dp
  dydp[1] = -std::sqrt(reactor.GetCv() / (cp_rho * p));

  return 0;
}

void setEnergiesHOrU(FlameMasterReactor& reactor, double target, double dTtol,
                     bool HOrU) {
  double dT{};

  double Tnew = reactor.GetTemperature();
  double Unew = HOrU ? reactor.GetEnthalpy() : reactor.GetInternalEnergy();
  double Cvnew = HOrU ? reactor.GetCp() : reactor.GetCv();

  double Utop = Unew;
  double Ubot = Unew;
  double Ttop = Tnew;
  double Tbot = Tnew;

  bool unstablePhase = false;
  double UConvErr{};

  // Newton iteration
  // This is exactly like Cantera's setState_UV implementation,
  // except that we do not check for upper/lower temperature limits
  for (int i = 0; i < 1000; i++) {
    const double Told = Tnew;
    const double Uold = Unew;

    double cvd = Cvnew;
    if (cvd < 0) {
      unstablePhase = true;
    }

    dT = std::max(-100., std::min(100., (target - Uold) / cvd));
    Tnew = Told + dT;

    // This limits the step size to make the algorithm cvergent
    // See Cantera for details
    if ((dT > 0 && unstablePhase) || (dT <= 0 && !unstablePhase)) {
      if (Ubot < target && Tnew < (0.75 * Tbot + 0.25 * Told)) {
        dT = .75 * (Tbot - Told);
      }
    } else if (Utop > target && Tnew > (.75 * Ttop + .25 * Told)) {
      dT = .75 * (Ttop - Told);
      Tnew = Told + dT;
    }

    // Set the new temperature, but try to stay in the stable region
    // with cv > 0
    for (int its = 0; its < 10; its++) {
      Tnew = Told + dT;
      reactor.SetTemperature(Tnew);

      Unew = HOrU ? reactor.GetEnthalpy() : reactor.GetInternalEnergy();
      Cvnew = HOrU ? reactor.GetCp() : reactor.GetCv();

      if (Cvnew < 0) {
        dT *= .25;
      } else {
        break;
      }
    }

    if (Unew == target) {
      return;
    } else if (Unew > target && (Utop < target || Unew < Utop)) {
      Utop = Unew;
      Ttop = Tnew;
    } else if (Unew < target && (Ubot > target || Unew > Ubot)) {
      Ubot = Unew;
      Tbot = Tnew;
    }

    // Check for convergence
    double Uerr = target - Unew;
    double acvd = std::max(std::abs(cvd), 1e-5);
    double denom = std::max(std::abs(target), acvd * dTtol);
    UConvErr = std::abs(Uerr / denom);
    if (UConvErr < 1e-5 * dTtol || std::abs(dT) < dTtol) {
      return;
    }
  }

  std::ostringstream errInfo;
  errInfo << "Failed to converge to given " << (HOrU ? "enthalpy" : "energy")
          << ". Missed target " << target << "J/m^3 by "
          << (target -
              (HOrU ? reactor.GetEnthalpy() : reactor.GetInternalEnergy()))
          << "J/m^3 (relerr = " << UConvErr
          << "). Stopped Newton iteration with dT = " << dT << ".";

  throw FlameMasterReactorException(errInfo.str());
}
} // namespace

FlameMasterReactor::FlameMasterReactor(const FlameMasterMechanism& mechanism)
    : mechanism_{&mechanism}, state_{}, ode_solver_{
                                            std::make_unique<RadauSolver>()} {
  state_.nSpecies = mechanism_->getNSpecs();
  state_.nSpeciesEffective = mechanism_->getNSpecies();
  if (state_.nSpecies == 0) {
    state_.nSpecies = state_.nSpeciesEffective;
  }
  state_.nReactions = mechanism_->getNReactions();
  state_.nThirdBodyReactions = mechanism_->getNThirdBodyReactions();
  state_.gasConstant = mechanism_->getUniversalGasConstant();

  // Allocate memory for the mass fractions and mole count
  // We store the temperature at the beginning of that vector such that we can
  // pass it directly to the ODE solver
  state_.massFractions.resize(static_cast<std::size_t>(state_.nSpecies));
  state_.molesStorage.resize(
      static_cast<std::size_t>(state_.nSpeciesEffective + 1));
  state_.temperature = state_.molesStorage.data();
  state_.moles = make_span(state_.molesStorage).subspan(1);

  array_state_.massFractions.resize(state_.nSpecies, kDefaultChunkSize);
  array_state_.moles.resize(state_.nSpecies, kDefaultChunkSize);
  array_state_.molesStorage.resize(state_.nSpeciesEffective + 1,
                                   kDefaultChunkSize);

  // Computational space for calls to the chemistry implementation
  state_.production_rates.resize(
      static_cast<std::size_t>(state_.nSpeciesEffective));
  state_.reaction_rates.resize(static_cast<std::size_t>(state_.nReactions));
  state_.rate_coefficients.resize(static_cast<std::size_t>(state_.nReactions));
  state_.third_body_concentrations.resize(
      static_cast<std::size_t>(state_.nThirdBodyReactions));
  state_.enthalpies.resize(static_cast<std::size_t>(state_.nSpeciesEffective));
  state_.heat_capacities_at_constant_pressure.resize(
      static_cast<std::size_t>(state_.nSpeciesEffective));
  state_.entropies.resize(static_cast<std::size_t>(state_.nSpeciesEffective));

  array_state_.production_rates.resize(state_.nSpeciesEffective,
                                       kDefaultChunkSize);
  array_state_.reaction_rates.resize(state_.nReactions, kDefaultChunkSize);
  array_state_.rate_coefficients.resize(state_.nReactions, kDefaultChunkSize);
  array_state_.third_body_concentrations.resize(state_.nThirdBodyReactions,
                                                kDefaultChunkSize);
  array_state_.enthalpies.resize(state_.nSpeciesEffective, kDefaultChunkSize);
  array_state_.heat_capacities_at_constant_pressure.resize(
      state_.nSpeciesEffective, kDefaultChunkSize);

  array_state_.entropies.resize(state_.nSpeciesEffective, kDefaultChunkSize);

  // Load the names of the species
  state_.speciesNames = mechanism_->getSpeciesNames();

  // get the molar masses of the species
  state_.molarMasses.resize(static_cast<std::size_t>(state_.nSpeciesEffective));
  mechanism_->getMolarMass(state_.molarMasses);

  // Initialize the state with something quite boring.
  state_.massFractions[0] = 1;
  UpdateMolesFromMassFractions(state_);
  SetTemperature(300.);
  SetPressure(101325.);
  array_state_.massFractions.row(0) = Array1d::Constant(1.0);
  SetTemperatureArray(Array1d::Constant(300.0));
  SetPressureArray(Array1d::Constant(101325.0));

  UpdateThermoState(*mechanism_, state_);
  UpdateThermoState(*mechanism_, array_state_, array_state_.temperature);
}

FlameMasterReactor::FlameMasterReactor(const FlameMasterReactor& other)
    : mechanism_{other.mechanism_},
      state_(other.state_), array_state_{other.array_state_} {
  state_.temperature = state_.molesStorage.data();
  state_.moles = make_span(state_.molesStorage).subspan(1);
  ode_solver_ = other.ode_solver_->Clone();
}

FlameMasterReactor::FlameMasterReactor(FlameMasterReactor&& other) noexcept
    : mechanism_{std::move(other.mechanism_)}, state_(std::move(other.state_)),
      array_state_(std::move(other.array_state_)),
      ode_solver_(std::move(other.ode_solver_)) {
  state_.temperature = state_.molesStorage.data();
  state_.moles = make_span(state_.molesStorage).subspan(1);
}

FlameMasterReactor& FlameMasterReactor::
operator=(FlameMasterReactor&& other) noexcept {
  mechanism_ = std::move(other.mechanism_);
  state_ = std::move(other.state_);
  ode_solver_ = std::move(other.ode_solver_);
  state_.temperature = state_.molesStorage.data();
  state_.moles = make_span(state_.molesStorage).subspan(1);
  return *this;
}

FlameMasterReactor& FlameMasterReactor::
operator=(const FlameMasterReactor& other) {
  FlameMasterReactor tmp(other);
  return *this = std::move(tmp);
}

struct AdvanceSystem {
  FlameMasterReactor* reactor;
  FlameMasterState* state;
  void operator()(span<double> dXdt, span<const double> X, double t) {
    OdeRhsAdvance(t, X, dXdt, reactor->getMechanism(), *state);
  }
};

struct AdvanceFeedback {
  FlameMasterReactor* reactor;
  // FlameMasterState* state;
  function_ref<void(span<const double>, double, FlameMasterReactor*)> feedback;

  void operator()(span<const double> y, double t) const noexcept {
    // if (y.data() != state->molesStorage.data()) {
    // FUB_ASSERT(y.size() == state->molesStorage.size());
    // reactor->setMoleFractions(y.subspan(1));
    // reactor->SetTemperature(y[0]);
    // }
    feedback(y, t, reactor);
  }
};

void FlameMasterReactor::Advance(double dt) {
  if (dt == 0) {
    return;
  }
  FUB_ASSERT(dt > 0);

  UpdateMolesFromMassFractions(state_);

  // Do the actual computation
  ode_solver_->Integrate(AdvanceSystem{this, &state_},
                         make_span(state_.molesStorage), 0.0, dt);

  // Retrieve mass fractions from the result vector
  UpdateMassFractionsFromMoles(state_);
}

void FlameMasterReactor::Advance(
    double dt,
    function_ref<void(span<const double>, double, FlameMasterReactor*)>
        feedback) {
  if (dt == 0) {
    return;
  }
  FUB_ASSERT(dt > 0);

  UpdateMolesFromMassFractions(state_);

  // Do the actual computation
  ode_solver_->Integrate(AdvanceSystem{this, &state_},
                         make_span(state_.molesStorage), 0.0, dt,
                         AdvanceFeedback{this, feedback});

  // Retrieve mass fractions from the result vector
  UpdateMassFractionsFromMoles(state_);
}

struct AdvanceAndFindMaxdTData {
  double oldTime;
  double oldTemperature;
  double maxdT;
  double argmaxdT;
};

double FlameMasterReactor::AdvanceAndFindMaxdT(double dt) {
  struct AdvanceAndFindMaxdT {
    AdvanceAndFindMaxdTData* data;

    void operator()(span<const double> y, double t) const noexcept {
      const double newTemperature = y[0];
      const double dT =
          (data->oldTemperature - newTemperature) / (data->oldTime - t);
      data->oldTime = t;
      if (data->maxdT < dT) {
        data->maxdT = dT;
        data->argmaxdT = t;
      }
    }
  };
  AdvanceAndFindMaxdTData data{0., GetTemperature(), -9000, -1};

  // Do the actual computation
  ode_solver_->Integrate(AdvanceSystem{this, &state_},
                         make_span(state_.molesStorage), 0.0, dt,
                         AdvanceAndFindMaxdT{&data});

  // Retrieve mass fractions from the result vector
  UpdateMassFractionsFromMoles(state_);
  return data.argmaxdT;
}

//
// Thermodynamic properties
//
// getters and setters
//

span<const double> FlameMasterReactor::GetMoleFractions() {
  UpdateMolesFromMassFractions(state_);
  return span<const double>(state_.molesStorage).subspan(1);
}

span<const double> FlameMasterReactor::GetCps() const {
  return state_.heat_capacities_at_constant_pressure;
}

double FlameMasterReactor::GetEnthalpy() const {
  return MeanY(state_.enthalpies);
}

Array1d FlameMasterReactor::GetEnthalpyArray() const {
  return MeanY(array_state_.enthalpies);
}

span<const double> FlameMasterReactor::GetEnthalpies() {
  UpdateThermoState(*mechanism_, state_, GetTemperature());
  return state_.enthalpies;
}

double FlameMasterReactor::GetCp() const {
  return MeanY(state_.heat_capacities_at_constant_pressure);
}

Array1d FlameMasterReactor::GetCpArray() const {
  return MeanY(array_state_.heat_capacities_at_constant_pressure);
}

double FlameMasterReactor::GetCv() const {
  return GetCp() - GetUniversalGasConstant() / GetMeanMolarMass();
}

Array1d FlameMasterReactor::GetCvArray() const {
  const Array1d mean_molar_mass = GetMeanMolarMassArray();
  const Array1d mean_molar_mass_s = (mean_molar_mass > 0.0).select(mean_molar_mass, 1.0);
  const Array1d R = Array1d::Constant(GetUniversalGasConstant());
  return (mean_molar_mass > 0.0).select(GetCpArray() - R / mean_molar_mass_s, 0.0);
}

double FlameMasterReactor::GetMeanMolarMass() const {
  const double minv = std::inner_product(
      state_.massFractions.begin(), state_.massFractions.end(),
      state_.molarMasses.begin(), 0.0, std::plus<>{}, std::divides<>{});
  return 1.0 / minv;
}

Array1d FlameMasterReactor::GetMeanMolarMassArray() const {
  Array1d minv = Array1d::Zero();
  for (int i = 0; i < array_state_.massFractions.rows(); ++i) {
    minv += array_state_.massFractions.row(i) / state_.molarMasses[i];
  }
  minv = (minv > 0.0).select(minv, 1.0);
  Array1d mean_molar_mass = (minv > 0.0).select(Array1d::Constant(1.0) / minv, 0.0);
  return mean_molar_mass;
}

double FlameMasterReactor::GetEntropy() const {
  span<const double> s = state_.entropies;
  if (s[0] == 0) {
    throw FlameMasterReactorException(
        "The mechanism you are using does not support entropy. (As indicated "
        "by the first species's entropy being zero.)");
  }
  // Compute entropy
  //
  // For entropy in ideal gas mixtures,
  //  c_p dT = T dS + RT/p dP
  // Also, there is entropy of mixing.

  const double meanM = GetMeanMolarMass();
  const double R = GetUniversalGasConstant();
  // Initialize entropy with the pressure part
  double entropy = MeanY(s) - R * std::log(GetPressure() / 101325.) / meanM;

  // Sum up the species' entropies and entropy of mixing
  double inv = 0;
  double moleSum = 0;
  const double eps = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < state_.nSpecies; i++) {
    inv += (state_.moles[i] > eps) ? state_.moles[i] * std::log(state_.moles[i])
                                   : 0;
    moleSum += state_.moles[i];
  }
  inv = inv / moleSum - std::log(moleSum);
  entropy -= GetUniversalGasConstant() * inv / meanM;
  return entropy;
}

void FlameMasterReactor::SetTemperatureArray(Array1d temperature) {
  array_state_.temperature = temperature;
  UpdateThermoState(*mechanism_, array_state_, temperature);
}

void FlameMasterReactor::SetEnthalpy(double enthalpy, double dTtol) {
  setEnergiesHOrU(*this, enthalpy, dTtol, true);
}

double FlameMasterReactor::GetInternalEnergy() const {
  double intEnergy = GetEnthalpy();
  double minv = 0;
  for (std::size_t i = 0; i < static_cast<std::size_t>(state_.nSpecies); i++) {
    minv += state_.massFractions[i] / state_.molarMasses[i];
  }
  intEnergy -= GetUniversalGasConstant() * minv * GetTemperature();
  return intEnergy;
}

Array1d FlameMasterReactor::GetInternalEnergyArray() const {
  Array1d intEnergy = GetEnthalpyArray();
  Array1d minv = Array1d::Zero();
  for (int i = 0; i < array_state_.massFractions.rows(); ++i) {
    minv += array_state_.massFractions.row(i) / state_.molarMasses[i];
  }
  intEnergy -= GetUniversalGasConstant() * minv * GetTemperatureArray();
  return intEnergy;
}

void FlameMasterReactor::SetInternalEnergy(double energy, double dTtol) {
  setEnergiesHOrU(*this, energy, dTtol, false);
}

void FlameMasterReactor::SetInternalEnergyArray(Array1d target, double dTtol) {
  Array1d dT = Array1d::Zero();

  Array1d Tnew = GetTemperatureArray();
  Array1d Unew = GetInternalEnergyArray();
  Array1d Cvnew = GetCvArray();

  Array1d Utop = Unew;
  Array1d Ubot = Unew;
  Array1d Ttop = Tnew;
  Array1d Tbot = Tnew;

  Array<bool, 1> unstablePhase = Array<bool, 1>::Constant(false);
  Array1d UConvErr = Array1d::Zero();

  // Newton iteration
  // This is exactly like Cantera's setState_UV implementation,
  // except that we do not check for upper/lower temperature limits
  for (int i = 0; i < 1000; i++) {
    const Array1d Told = Tnew;
    const Array1d Uold = Unew;

    Array1d cvd = Cvnew;
    unstablePhase = (cvd < 0.0);

    dT = ((target - Uold) / cvd).max(-100.).min(+100.);
    Tnew = Told + dT;
    dT = ((dT > 0.0 && unstablePhase) || (dT <= 0.0 && !unstablePhase))
             .select((Ubot < target && Tnew < (0.75 * Tbot + 0.25 * Told))
                         .select(0.75 * (Tbot - Told), dT),
                     (Utop > target && Tnew > (.75 * Ttop + .25 * Told))
                         .select(0.75 * (Ttop - Told), dT));

    // Set the new temperature, but try to stay in the stable region
    // with cv > 0
    for (int its = 0; its < 10; its++) {
      Tnew = Told + dT;
      SetTemperatureArray(Tnew);
      Unew = GetInternalEnergyArray();
      Cvnew = GetCvArray();
      if ((Cvnew < 0.0).any()) {
        dT = (Cvnew < 0.0).select(0.25 * dT, dT);
      } else {
        break;
      }
    }

    if ((Unew == target).all()) {
      return;
    }

    Array<bool, 1> update_top =
        (Unew > target && (Utop < target || Unew < Utop));
    Utop = update_top.select(Unew, Utop);
    Ttop = update_top.select(Tnew, Ttop);

    Array<bool, 1> update_bot =
        (Unew < target && (Ubot > target || Unew > Ubot));
    Ubot = update_bot.select(Unew, Ubot);
    Tbot = update_bot.select(Tnew, Tbot);

    // Check for convergence
    Array1d Uerr = target - Unew;
    Array1d acvd = cvd.abs().max(1e-5);
    Array1d denom = target.abs().max(acvd * dTtol);
    UConvErr = Uerr.abs();
    if ((UConvErr < denom * (1e-5 * dTtol) || dT.abs() < dTtol).all()) {
      return;
    }
  }

  std::ostringstream errInfo;
  errInfo << "Failed to converge to given energy. "
          << "Missed target [" << target << "] J/m^3 by ["
          << target - GetInternalEnergyArray()
          << "] J/m^3 (relerr = " << UConvErr
          << "). Stopped Newton iteration with dT = " << dT << ".";

  throw FlameMasterReactorException(errInfo.str());
}

void FlameMasterReactor::SetInternalEnergyArray(Array1d target, MaskArray mask,
                                                double dTtol) {
  Array1d dT = Array1d::Zero();

  target = mask.select(target, 0.0);
  Array1d Tnew = mask.select(GetTemperatureArray(), 0.0);
  Array1d Unew = mask.select(GetInternalEnergyArray(), 0.0);
  Array1d Cvnew = mask.select(GetCvArray(), 0.0);

  Array1d Utop = Unew;
  Array1d Ubot = Unew;
  Array1d Ttop = Tnew;
  Array1d Tbot = Tnew;

  Array<bool, 1> unstablePhase = Array<bool, 1>::Constant(false);
  Array1d UConvErr = Array1d::Zero();

  // Newton iteration
  // This is exactly like Cantera's setState_UV implementation,
  // except that we do not check for upper/lower temperature limits
  for (int i = 0; i < 1000; i++) {
    const Array1d Told = Tnew;
    const Array1d Uold = Unew;

    Array1d cvd = mask.select(Cvnew, 1.0);
    unstablePhase = (cvd < 0.0);

    dT = ((target - Uold) / cvd).max(-100.).min(+100.);
    Tnew = Told + dT;
    dT = ((dT > 0.0 && unstablePhase) || (dT <= 0.0 && !unstablePhase))
             .select((Ubot < target && Tnew < (0.75 * Tbot + 0.25 * Told))
                         .select(0.75 * (Tbot - Told), dT),
                     (Utop > target && Tnew > (.75 * Ttop + .25 * Told))
                         .select(0.75 * (Ttop - Told), dT));

    // Set the new temperature, but try to stay in the stable region
    // with cv > 0
    for (int its = 0; its < 10; its++) {
      Tnew = mask.select(Told + dT, 300.0);
      SetTemperatureArray(Tnew);
      Unew = mask.select(GetInternalEnergyArray(), 0.0);
      Cvnew = mask.select(GetCvArray(), 0.0);
      if ((Cvnew < 0.0).any()) {
        dT = (Cvnew < 0.0).select(0.25 * dT, dT);
      } else {
        break;
      }
    }

    if ((Unew == target).all()) {
      return;
    }

    Array<bool, 1> update_top =
        (Unew > target && (Utop < target || Unew < Utop));
    Utop = update_top.select(Unew, Utop);
    Ttop = update_top.select(Tnew, Ttop);

    Array<bool, 1> update_bot =
        (Unew < target && (Ubot > target || Unew > Ubot));
    Ubot = update_bot.select(Unew, Ubot);
    Tbot = update_bot.select(Tnew, Tbot);

    // Check for convergence
    Array1d Uerr = mask.select(target - Unew, 0.0);
    Array1d acvd = cvd.abs().max(1e-5);
    Array1d denom = mask.select(target.abs().max(acvd * dTtol), 1.0);
    UConvErr = Uerr.abs();
    if ((!mask || UConvErr < denom * (1e-5 * dTtol) || dT.abs() < dTtol).all()) {
      return;
    }
  }

  std::ostringstream errInfo;
  errInfo << "Failed to converge to given energy. "
          << "Missed target [" << target << "] J/m^3 with mask [" << mask
          << "] by [" << target - GetInternalEnergyArray()
          << "] J/m^3 (relerr = " << UConvErr
          << "). Stopped Newton iteration with dT = " << dT << ".";

  throw FlameMasterReactorException(errInfo.str());
}

//
// Species
//

double FlameMasterReactor::MeanX(span<const double> quantity) const {
  double meanWeight = GetMeanMolarMass();
  double retval = 0;
  for (std::size_t i = 0; i < static_cast<std::size_t>(state_.nSpecies); i++) {
    retval += quantity[static_cast<std::ptrdiff_t>(i)] *
              state_.massFractions[i] * meanWeight / state_.molarMasses[i];
  }
  return retval;
}

double FlameMasterReactor::MeanY(span<const double> quantity) const {
  double retval = 0;
  for (int i = 0; i < state_.nSpecies; i++) {
    retval += quantity[i] * state_.massFractions[static_cast<std::size_t>(i)];
  }
  return retval;
}

Array1d FlameMasterReactor::MeanY(const ArrayXd& quantity) const {
  Array1d retval = Array1d::Zero();
  for (int i = 0; i < array_state_.massFractions.rows(); i++) {
    retval += quantity.row(i) * array_state_.massFractions.row(i);
  }
  return retval;
}

void FlameMasterReactor::SetMassFractions(std::string newMassFractions) {
  std::fill(state_.massFractions.begin(), state_.massFractions.end(), 0.0);
  std::istringstream inputParser(newMassFractions);

  inputParser >> std::skipws;
  inputParser.imbue(std::locale("C"));

  double speciesSum = 0;

  // Parse mass fraction string in the format
  //  species: fraction, ...
  //
  while (inputParser.good()) {
    std::string speciesName;
    double speciesValue;

    while (inputParser.peek() != ':' && inputParser.good())
      speciesName += (char)inputParser.get();
    if (inputParser.get() != ':')
      throw FlameMasterReactorException("Invalid mass fraction string: " +
                                        newMassFractions);
    inputParser >> speciesValue;
    while (inputParser.peek() == ',' || inputParser.peek() == ' ')
      inputParser.get();

    bool found = false;
    for (std::size_t i = 0; i < static_cast<std::size_t>(state_.nSpecies);
         i++) {
      if (speciesName == GetSpeciesName(static_cast<int>(i))) {
        state_.massFractions[i] = speciesValue;
        speciesSum += speciesValue;
        found = true;
        break;
      }
    }

    if (!found) {
      throw FlameMasterReactorException("Unknown species: " + speciesName);
    }
  }
  if (!inputParser.eof())
    throw FlameMasterReactorException("Invalid mass fraction string: " +
                                      newMassFractions);

  for (std::size_t i = 0; i < static_cast<std::size_t>(state_.nSpecies); i++) {
    state_.massFractions[i] /= speciesSum;
  }
  UpdateMolesFromMassFractions(state_);
}

void FlameMasterReactor::SetMassFractions(span<const double> fractions) {
  std::transform(fractions.begin(), fractions.end(),
                 state_.massFractions.begin(),
                 [](double Y) { return std::max(0.0, Y); });
  const double sum = std::accumulate(state_.massFractions.begin(),
                                     state_.massFractions.end(), 0.0);
  if (sum != 0) {
    std::transform(state_.massFractions.begin(), state_.massFractions.end(),
                   state_.massFractions.begin(),
                   [sum](double frac) { return frac / sum; });
  }
  UpdateMolesFromMassFractions(state_);
}

void FlameMasterReactor::SetMassFractionsArray(const ArrayXd& newMassFractions, MaskArray mask) {
  Array1d sum = Array1d::Zero();
  Array1d Y0 = mask.select(newMassFractions.row(0), 1.0);
  array_state_.massFractions.row(0) = Y0.max(0.0);
  sum += array_state_.massFractions.row(0);
  for (int i = 1; i < array_state_.massFractions.rows(); ++i) {
    Array1d Yi = mask.select(newMassFractions.row(i), 0.0);
    array_state_.massFractions.row(i) = Yi.max(0.0);
    sum += array_state_.massFractions.row(i);
  }
  sum = (sum > 0.0).select(sum, 1.0);
  for (int i = 0; i < array_state_.massFractions.rows(); ++i) {
    array_state_.massFractions.row(i) /= sum;
  }
  UpdateMolesFromMassFractions(array_state_, state_);
}

void FlameMasterReactor::SetMassFractionsArray(
    const ArrayXd& newMassFractions) {
  Array1d sum = Array1d::Zero();
  for (int i = 0; i < array_state_.massFractions.rows(); ++i) {
    array_state_.massFractions.row(i) = newMassFractions.row(i).abs();
    sum += array_state_.massFractions.row(i);
  }
  sum = (sum > 0.0).select(sum, 1.0);
  for (int i = 0; i < array_state_.massFractions.rows(); ++i) {
    array_state_.massFractions.row(i) /= sum;
  }
  UpdateMolesFromMassFractions(array_state_, state_);
}

void FlameMasterReactor::SetTemperature(double temp) {
  *state_.temperature = temp;
  UpdateThermoState(*mechanism_, state_, temp);
}

void FlameMasterReactor::SetMoleFractions(span<const double> newMoleFractions) {
  SetMassFractions(newMoleFractions);
  double totalMoles = 0.0;
  for (std::size_t i = 0; i < static_cast<std::size_t>(state_.nSpecies); i++) {
    totalMoles += state_.massFractions[i] * state_.molarMasses[i];
  }
  for (std::size_t i = 0; i < static_cast<std::size_t>(state_.nSpecies); i++) {
    state_.massFractions[i] *= state_.molarMasses[i] / totalMoles;
  }
  UpdateMolesFromMassFractions(state_);
}

void FlameMasterReactor::SetMoleFractions(std::string newMoleFractions) {
  SetMassFractions(newMoleFractions);
  double totalMoles = 0.0;
  for (std::size_t i = 0; i < static_cast<std::size_t>(state_.nSpecies); i++) {
    totalMoles += state_.massFractions[i] * state_.molarMasses[i];
  }
  for (std::size_t i = 0; i < static_cast<std::size_t>(state_.nSpecies); i++) {
    state_.massFractions[i] *= state_.molarMasses[i] / totalMoles;
  }
  UpdateMolesFromMassFractions(state_);
}

span<const double> FlameMasterReactor::GetReactionRates() const {
  return state_.reaction_rates;
}

span<const double> FlameMasterReactor::GetProductionRates() {
  mechanism_->ComputeProductionRates(state_.production_rates,          // cdot
                                     state_.reaction_rates,            // w
                                     state_.rate_coefficients,         // k
                                     state_.moles,                     // c
                                     state_.third_body_concentrations, // M
                                     GetTemperature(), // Temperature
                                     GetPressure()     // Pressure
  );
  return state_.production_rates;
}

namespace {
struct SetIsentropicPSystem {
  FlameMasterReactor* reactor;
  void operator()(span<double> dydp, span<const double> y, double p) {
    OdeRhsSetP(dydp, y, p, *reactor);
  }
};
} // namespace

double FlameMasterReactor::SetPressureIsentropic(double pressure) {
  // Do nothing if the pressures match
  if (std::abs(pressure - GetPressure()) < 1e-5) {
    return 0;
  }

  state_.setPVector[0] = GetTemperature();
  state_.setPVector[1] = 0;

  // Do the actual computation
  ode_solver_->Integrate(SetIsentropicPSystem{this},
                         make_span(state_.setPVector), GetPressure(),
                         pressure - GetPressure());

  // Set the state of the reactor to the values from the result vector
  SetTemperature(state_.setPVector[0]);
  SetPressure(pressure);

  return state_.setPVector[1];
}

} // namespace fub
