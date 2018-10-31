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

#include "fub/solver/euler/Reactor.hpp"
#include "src/solver/ode_solver/Radau.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

namespace fub {
namespace {
void UpdateMassFractionsFromMoles(FlameMasterState& state) noexcept {
  for (int i = 0; i < state.nSpecies; i++) {
    state.massFractions[i] =
        state.moles[i] * state.molarMasses[i] / state.density;
  }
}

void UpdateMolesFromMassFractions(FlameMasterState& state) noexcept {
  for (int i = 0; i < state.nSpecies; i++) {
    state.moles[i] =
        state.massFractions[i] * state.density / state.molarMasses[i];
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
                       FlameMasterState& state) {
  UpdateThermoState(mechanism, state, *state.temperature);
}

int OdeRhsAdvance(double t, span<const double> XD, span<double> Xdot,
                  const FlameMasterMechanism& mechanism,
                  FlameMasterState& state) {
  span<double> h = state.enthalpies;
  span<double> cp = state.heat_capacities_at_constant_pressure;
  double temperature = XD[0];

  // We need to work in a local variable rather than in XD, because
  // the QSS code might alter XD beyond its size (nSpecies), assuming it
  // had space for nSpeciesEffective
  int size = XD.size() - 1;
  span<double> concentrations((double*)alloca(sizeof(double) * size), size);
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

void SetEnergiesHOrU(FlameMasterReactor& reactor, double target, double dTtol,
                     bool HOrU) {
  double dT;

  double Tnew = reactor.GetTemperature();
  double Unew = HOrU ? reactor.GetEnthalpy() : reactor.GetInternalEnergy();
  double Cvnew = HOrU ? reactor.GetCp() : reactor.GetCv();

  double Utop = Unew;
  double Ubot = Unew;
  double Uold = Unew;
  double Ttop = Tnew;
  double Tbot = Tnew;
  double Told = Tnew;

  bool unstablePhase = false;
  double Tunstable = -1;

  double UConvErr;

  // Newton iteration
  // This is exactly like Cantera's setState_UV implementation,
  // except that we do not check for upper/lower temperature limits
  for (int i = 0; i < 1000; i++) {
    Told = Tnew;
    Uold = Unew;

    double cvd = Cvnew;
    if (cvd < 0) {
      unstablePhase = true;
      Tunstable = Tnew;
    }

    dT = std::max(-100., std::min(100., (target - Uold) / cvd));
    Tnew = Told + dT;

    // This limits the step size to make the algorithm convergent
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
        Tunstable = Tnew;
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

FlameMasterReactor::FlameMasterReactor(
    std::unique_ptr<FlameMasterMechanism> mechanism)
    : mechanism_{std::move(mechanism)}, state_{} {
  state_.nSpecies = mechanism_->GetNSpecs();
  state_.nSpeciesEffective = mechanism_->GetNSpecies();
  if (state_.nSpecies == 0) {
    state_.nSpecies = state_.nSpeciesEffective;
  }
  state_.nReactions = mechanism_->GetNReactions();
  state_.nThirdBodyReactions = mechanism_->GetNThirdBodyReactions();
  state_.gasConstant = mechanism_->GetUniversalGasConstant();

  // Allocate memory for the mass fractions and mole count
  // We store the temperature at the beginning of that vector such that we can
  // pass it directly to the ODE solver
  state_.massFractions.resize(state_.nSpecies);
  state_.molesStorage.resize(state_.nSpeciesEffective + 1);
  state_.temperature = state_.molesStorage.data();
  state_.moles = span{state_.molesStorage}.subspan(1);

  // Computational space for calls to the chemistry implementation
  state_.production_rates.resize(state_.nSpeciesEffective);
  state_.reaction_rates.resize(state_.nReactions);
  state_.rate_coefficients.resize(state_.nReactions);
  state_.third_body_concentrations.resize(state_.nThirdBodyReactions);
  state_.enthalpies.resize(state_.nSpeciesEffective);
  state_.heat_capacities_at_constant_pressure.resize(state_.nSpeciesEffective);
  state_.entropies.resize(state_.nSpeciesEffective);

  // Load the names of the species
  state_.speciesNames = mechanism_->GetSpeciesNames();

  // Get the molar masses of the species
  state_.molarMasses.resize(state_.nSpeciesEffective);
  mechanism_->GetMolarMass(state_.molarMasses);

  // Initialize the state with something quite boring.
  state_.massFractions[0] = 1;
  SetTemperature(300.);
  SetPressure(101325.);

  UpdateThermoState(*mechanism_, state_);
  // mechanism_->ComputeProductionRates(state_.production_rates,          //
  // cdot
  //                                    state_.reaction_rates,            // w
  //                                    state_.rate_coefficients,         // k
  //                                    state_.moles,                     // c
  //                                    state_.third_body_concentrations, // M
  //                                    GetTemperature(), // Temperature
  //                                    GetPressure()     // Pressure
  // );
}

FlameMasterReactor::FlameMasterReactor(const FlameMasterReactor& other)
    : mechanism_{other.mechanism_->Clone()}, state_{other.state_} {
  state_.temperature = state_.molesStorage.data();
  state_.moles = span{state_.molesStorage}.subspan(1);
}

FlameMasterReactor::FlameMasterReactor(FlameMasterReactor&& other) noexcept
    : mechanism_{std::move(other.mechanism_)}, state_{std::move(other.state_)} {
  state_.temperature = state_.molesStorage.data();
  state_.moles = span{state_.molesStorage}.subspan(1);
}

struct AdvanceSystem {
  const FlameMasterMechanism* mechanism;
  FlameMasterState* state;
  int operator()(span<double> dXdt, span<const double> X, double t) {
    return OdeRhsAdvance(t, X, dXdt, *mechanism, *state);
  }
};

struct AdvanceFeedback {
  FlameMasterReactor* reactor;
  FlameMasterState* state;
  function_ref<int(double, FlameMasterReactor*)> feedback;

  int operator()(span<const double> y, double t) const noexcept {
    if (y.data() != state->molesStorage.data()) {
      FUB_ASSERT(y.size() == state->molesStorage.size());
      std::copy(y.begin(), y.end(), state->molesStorage.begin());
    }
    UpdateMassFractionsFromMoles(*state);
    return feedback(t, reactor);
  }
};

void FlameMasterReactor::Advance(double dt) {
  if (dt == 0) {
    return;
  }
  FUB_ASSERT(dt > 0);

  // Do the actual computation
  Radau::integrate(AdvanceSystem{mechanism_.get(), &state_}, state_.moles, 0.0,
                   dt);

  // Retrieve mass fractions from the result vector
  UpdateMassFractionsFromMoles(state_);
}

void FlameMasterReactor::Advance(
    double dt, function_ref<int(double, FlameMasterReactor*)> feedback) {
  if (dt == 0) {
    return;
  }
  FUB_ASSERT(dt > 0);

  // Do the actual computation
  Radau::integrate_feedback(AdvanceSystem{mechanism_.get(), &state_},
                            state_.moles, 0.0, dt,
                            AdvanceFeedback{this, &state_, feedback});

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

    int operator()(span<const double> y, double t) const noexcept {
      const double newTemperature = y[0];
      const double dT =
          (data->oldTemperature - newTemperature) / (data->oldTime - t);
      data->oldTime = t;
      if (data->maxdT < dT) {
        data->maxdT = dT;
        data->argmaxdT = t;
      }
      return 0;
    }
  };
  AdvanceAndFindMaxdTData data{0., GetTemperature(), -9000, -1};

  // Do the actual computation
  Radau::integrate_feedback(AdvanceSystem{mechanism_.get(), &state_},
                            state_.moles, 0.0, dt, AdvanceAndFindMaxdT{&data});

  // Retrieve mass fractions from the result vector
  UpdateMassFractionsFromMoles(state_);
  return data.argmaxdT;
}

//
// Thermodynamic properties
//
// Getters and setters
//

span<const double> FlameMasterReactor::GetCps() const {
  return state_.heat_capacities_at_constant_pressure;
}

double FlameMasterReactor::GetEnthalpy() const {
  return MeanY(state_.enthalpies);
}

double FlameMasterReactor::GetCp() const {
  return MeanY(state_.heat_capacities_at_constant_pressure);
}

double FlameMasterReactor::GetCv() const {
  return GetCp() - GetUniversalGasConstant() / GetMeanMolarMass();
}

double FlameMasterReactor::GetMeanMolarMass() const {
  double minv = 0;
  for (int i = 0; i < state_.nSpecies; i++) {
    minv += state_.massFractions[i] / state_.molarMasses[i];
  }
  return 1 / minv;
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

void FlameMasterReactor::SetEnthalpy(double enthalpy, double dTtol) {
  SetEnergiesHOrU(*this, enthalpy, dTtol, true);
}

double FlameMasterReactor::GetInternalEnergy() const {
  double intEnergy = GetEnthalpy();
  double minv = 0;
  for (int i = 0; i < state_.nSpecies; i++) {
    minv += state_.massFractions[i] / state_.molarMasses[i];
  }
  intEnergy -= GetUniversalGasConstant() * minv * GetTemperature();
  return intEnergy;
}

void FlameMasterReactor::SetInternalEnergy(double energy, double dTtol) {
  SetEnergiesHOrU(*this, energy, dTtol, false);
}

//
// Species
//

double FlameMasterReactor::MeanX(span<const double> quantity) const {
  double meanWeight = GetMeanMolarMass();
  double retval = 0;
  for (int i = 0; i < state_.nSpecies; i++) {
    retval += quantity[i] * state_.massFractions[i] * meanWeight /
              state_.molarMasses[i];
  }
  return retval;
}

double FlameMasterReactor::MeanY(span<const double> quantity) const {
  double retval = 0;
  for (int i = 0; i < state_.nSpecies; i++) {
    retval += quantity[i] * state_.massFractions[i];
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
    for (int i = 0; i < state_.nSpecies; i++) {
      if (speciesName == GetSpeciesName(i)) {
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

  for (int i = 0; i < state_.nSpecies; i++) {
    state_.massFractions[i] /= speciesSum;
  }
}

void FlameMasterReactor::SetMassFractions(span<const double> fractions) {
  std::copy(fractions.begin(), fractions.end(), state_.massFractions.begin());
  const double sum = std::accumulate(fractions.begin(), fractions.end(), 0.0);
  if (sum != 0) {
    std::transform(state_.massFractions.begin(), state_.massFractions.end(),
                   state_.massFractions.begin(),
                   [sum](double frac) { return frac / sum; });
  }
}

void FlameMasterReactor::SetTemperature(double temp) {
  *state_.temperature = temp;
  UpdateThermoState(*mechanism_, state_, temp);
}

void FlameMasterReactor::SetMoleFractions(span<const double> newMoleFractions) {
  SetMassFractions(newMoleFractions);
  double totalMoles = 0.0;
  for (int i = 0; i < state_.nSpecies; i++) {
    totalMoles += state_.massFractions[i] * state_.molarMasses[i];
  }
  for (int i = 0; i < state_.nSpecies; i++) {
    state_.massFractions[i] *= state_.molarMasses[i] / totalMoles;
  }
}

void FlameMasterReactor::SetMoleFractions(std::string newMoleFractions) {
  SetMassFractions(newMoleFractions);
  double totalMoles = 0.0;
  for (int i = 0; i < state_.nSpecies; i++) {
    totalMoles += state_.massFractions[i] * state_.molarMasses[i];
  }
  for (int i = 0; i < state_.nSpecies; i++) {
    state_.massFractions[i] *= state_.molarMasses[i] / totalMoles;
  }
}

span<const double> FlameMasterReactor::GetReactionRates() const {
  return state_.reaction_rates;
}

span<const double> FlameMasterReactor::GetProductionRates() {
  // We need to work in a local variable rather than in XD, because
  // the QSS code might alter XD beyond its size (nSpecies), assuming it
  // had space for nSpeciesEffective

  // Compute species' derivatives
  mechanism_->ComputeProductionRates(state_.production_rates,          // cdot
                                     state_.reaction_rates,            // w
                                     state_.rate_coefficients,         // k
                                     state_.moles,                     // c
                                     state_.third_body_concentrations, // M
                                     GetTemperature(), // Temperature
                                     GetPressure()     // Pressure
  );

  return state_.production_rates;
} // namespace fub

struct SetIsentropicPSystem {
  FlameMasterReactor* reactor;
  int operator()(span<double> dydp, span<const double> y, double p) {
    return OdeRhsSetP(dydp, y, p, *reactor);
  }
};

double FlameMasterReactor::SetPressureIsentropic(double pressure) {
  // Do nothing if the pressures match
  if (std::abs(pressure - GetPressure()) < 1e-5) {
    return 0;
  }

  state_.setPVector[0] = GetTemperature();
  state_.setPVector[1] = 0;

  // Do the actual computation
  Radau::integrate(SetIsentropicPSystem{this}, span{state_.setPVector},
                   GetPressure(), pressure - GetPressure());

  // Set the state of the reactor to the values from the result vector
  SetTemperature(state_.setPVector[0]);
  SetPressure(pressure);

  return state_.setPVector[1];
}

} // namespace fub