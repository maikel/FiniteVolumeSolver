#include "fub/core/mdspan.hpp"
#include "fub/core/span.hpp"
#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/mechanism/Zhao2008Dme.hpp"

extern "C" {
#include "TC_interface.h"
}

#include <numeric>

using namespace fub;
using namespace fub::ideal_gas;

void UpdateMassFractionsFromMoles(FlameMasterState& state) noexcept {
  for (int i = 0; i < state.nSpecies; i++) {
    state.massFractions[i] =
        state.moles[i] * state.molarMasses[i] / state.density;
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

template <typename Function>
void ComputeNumericJacobian(span<double> jac, Function f,
                            span<const double> x0) {
  const int dimensions = x0.size();
  DynamicMdSpan<double, 2> J(jac.data(), dimensions, dimensions);
  std::vector<double> y0(dimensions);
  f(y0, x0);
  std::vector<double> x(x0.begin(), x0.end());
  constexpr double tau = 1e-8;
  for (int i = 0; i < dimensions; ++i) {
    auto yi = jac.subspan(i * dimensions, dimensions);
    const double old_xi = x[i];
    x[i] += tau;
    f(yi, x);
    std::transform(
        y0.begin(), y0.end(), yi.begin(), yi.begin(),
        [=](double f0_xj, double fi_xj) { return (fi_xj - f0_xj) / tau; });
    x[i] = old_xi;
  }
}

void UpdateMolesFromMassFractions(FlameMasterState& state) noexcept {
  for (int i = 0; i < state.nSpecies; i++) {
    state.moles[i] =
        state.massFractions[i] * state.density / state.molarMasses[i];
  }
}

void SetTemperature(FlameMasterState& state, double temp) {
  *state.temperature = temp;
}

double GetTemperature(FlameMasterState& state) { return *state.temperature; }

double GetMeanMolarMass(const FlameMasterState& state) {
  return std::inner_product(state.massFractions.begin(),
                            state.massFractions.end(),
                            state.molarMasses.begin(), 0.0);
}

double GetMeanCp(const FlameMasterState& state) {
  return std::inner_product(
      state.massFractions.begin(), state.massFractions.end(),
      state.heat_capacities_at_constant_pressure.begin(), 0.0);
}

void SetPressure(FlameMasterState& state, double pressure) {
  const double M = GetMeanMolarMass(state);
  const double T = GetTemperature(state);
  const double R = 8314.462175;
  state.density = pressure * M / T / R;
}

int main() {
  Zhao2008Dme mech{};
  FlameMasterState state{};

  state.nSpecies = mech.getNSpecs();
  state.nSpeciesEffective = mech.getNSpecies();
  if (state.nSpecies == 0) {
    state.nSpecies = state.nSpeciesEffective;
  }
  state.nReactions = mech.getNReactions();
  state.nThirdBodyReactions = mech.getNThirdBodyReactions();
  state.gasConstant = mech.getUniversalGasConstant();

  // Allocate memory for the mass fractions and mole count
  // We store the temperature at the beginning of that vector such that we can
  // pass it directly to the ODE solver
  state.massFractions.resize(state.nSpecies);
  state.molesStorage.resize(state.nSpeciesEffective + 1);
  state.temperature = state.molesStorage.data();
  state.moles = make_span(state.molesStorage).subspan(1);

  // Computational space for calls to the chemistry implementation
  state.production_rates.resize(state.nSpeciesEffective);
  state.reaction_rates.resize(state.nReactions);
  state.rate_coefficients.resize(state.nReactions);
  state.third_body_concentrations.resize(state.nThirdBodyReactions);
  state.enthalpies.resize(state.nSpeciesEffective);
  state.heat_capacities_at_constant_pressure.resize(state.nSpeciesEffective);
  state.entropies.resize(state.nSpeciesEffective);

  // Load the names of the species
  state.speciesNames = mech.getSpeciesNames();

  // get the molar masses of the species
  state.molarMasses.resize(state.nSpeciesEffective);
  mech.getMolarMass(state.molarMasses);

  // Initialize the state with something quite boring.
  state.moles[Zhao2008Dme::sH2] = 2;
  state.moles[Zhao2008Dme::sO2] = 1;
  *state.temperature = 1100.;
  UpdateThermoState(mech, state, 1100.);
  const double total =
      std::accumulate(state.moles.begin(), state.moles.end(), 0.0);
  std::transform(state.moles.begin(), state.moles.end(), state.moles.begin(),
                 [rho = state.density, total](double X) { return X / total; });
  state.density = 1;
  UpdateMassFractionsFromMoles(state);
  SetPressure(state, 101325.);
  UpdateMolesFromMassFractions(state);

  int nvars = mech.getNSpecies() + 1;
  std::vector<double> jac(nvars * nvars);
  ComputeNumericJacobian(jac,
                         [&](span<double> y, span<const double> x) {
                           OdeRhsAdvance(0, x, y, mech, state);
                         },
                         state.molesStorage);

  DynamicMdSpan<double, 2> J(jac.data(), nvars, nvars);
  std::printf("%f %f %f %f %f\n", J(0, 0), J(1, 0), J(0, 1), J(2, 0), J(0, 2));
  std::printf("%f %f %f\n", J(1, 1), J(1, 2), J(2, 1));

  const char* mechfile = mech.GetMechanismResource();
  const char* thermofile = mech.GetThermoResource();
  FILE* mechin = fmemopen(const_cast<char*>(mechfile), strlen(mechfile), "r");
  FILE* thermoin =
      fmemopen(const_cast<char*>(thermofile), strlen(thermofile), "r");
  FUB_ASSERT(mechin && thermofile);
  constexpr int with_tab = 0;
  constexpr double tab_delta_T = 1.0;
  TC_initChemFromFile(mechin, thermoin, with_tab, tab_delta_T);

  std::copy(state.massFractions.begin(), state.massFractions.end(),
            state.moles.begin());
  TC_setDens(state.density);
  span<double> TandYs(state.molesStorage.data(), nvars);
  TC_getJacCVTYNanl(TandYs.data(), nvars - 1, jac.data());

  span<double> M = state.molarMasses;
  const double R = 8314.462175;
  const double Constant =
      state.density * (R / GetMeanMolarMass(state) - GetMeanCp(state));

  std::printf("%f %f %f %f %f\n", J(0, 0), J(1, 0) / state.density * M[0],
              J(0, 1) * state.density / M[0], J(2, 0) / M[1], J(0, 2) / M[1]);

  std::printf("%f %f %f\n", J(1, 1), J(1, 2) * M[0] / M[1],
              J(2, 1) * M[1] / M[0]);
}