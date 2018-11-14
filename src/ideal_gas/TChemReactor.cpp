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

#include "fub/ideal_gas/TChemReactor.hpp"
#include "fub/ode_solver/RadauSolver.hpp"

#include "TC_defs.h"
#include "TC_interface.h"

#include <algorithm>
#include <numeric>

#include <cmath>
#include <cstdio>

namespace fub {
namespace ideal_gas {
namespace {
double GetUSpecMs_(span<const double> TandY, span<double> ui) {
  TC_getUspecMs(TandY[0], ui.size(), ui.data());
  span<const double> Y = TandY.subspan(1);
  return std::inner_product(Y.begin(), Y.end(), ui.begin(), 0.0);
}

double GetTemperatureFromInternalEnergy_(double target, span<double> TandY,
                                         span<double> uis, double dTtol) {
  double dT = 0.0;
  double Tnew = 300.0;
  double Cvnew = 0.0;
  TandY[0] = Tnew;
  double Unew = GetUSpecMs_(TandY, uis);
  TC_getMs2CvMixMs(TandY.data(), TandY.size(), &Cvnew);

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
      TandY[0] = Tnew;
      Unew = GetUSpecMs_(TandY, uis);
      TC_getMs2CvMixMs(TandY.data(), TandY.size(), &Cvnew);

      if (Cvnew < 0) {
        Tunstable = Tnew;
        dT *= .25;
      } else {
        break;
      }
    }

    if (Unew == target) {
      return Tnew;
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
      return Tnew;
    }
  }
  throw std::runtime_error(
      "TChemKinetics: No Convergence in GetTemperatureFromInternalEnergy");
}

void UpdateMoleFractionsFromMassFractions(span<double> moles,
                                          span<const double> masses) {
  TC_getMs2Ml(masses.data(), masses.size(), moles.data());
}

void UpdateMassFractionsFromMoleFractions(span<double> masses,
                                          span<const double> moles) {
  TC_getMl2Ms(moles.data(), moles.size(), masses.data());
}

} // namespace

TChemReactor::TChemReactor(const TChemMechanism& mechanism) {
  if (TC_isInit_) {
    throw std::runtime_error("TChem2 is already initialized.");
  }
  const char* mechfile = mechanism.GetMechanismResource();
  const char* thermofile = mechanism.GetThermoResource();
  FILE* mechin = fmemopen(const_cast<char*>(mechfile), strlen(mechfile), "r");
  FILE* thermoin =
      fmemopen(const_cast<char*>(thermofile), strlen(thermofile), "r");
  FUB_ASSERT(mechin && thermofile);
  constexpr int with_tab = 1;
  constexpr double tab_delta_T = 1.0;
  int ec = TC_initChemFromFile(mechin, thermoin, with_tab, tab_delta_T);
  FUB_ASSERT(!ec);
  FUB_ASSERT(TC_isInit_);
  const int n_species = TC_getNspec();
  temperature_and_mass_fractions_.resize(n_species + 1);
  moles_.resize(n_species);
  cps_.resize(n_species);
  cvs_.resize(n_species);
  enthalpies_.resize(n_species);
  internal_energies_.resize(n_species);
  temperature_and_mass_fractions_[1] = 1.0;
  SetMassFractions(GetMassFractions());
  SetPressure(101325.0);
  SetTemperature(300);
  ode_solver_ = std::make_unique<RadauSolver>();
}

TChemReactor::~TChemReactor() noexcept {
  FUB_ASSERT(TC_isInit_);
  TC_reset();
  FUB_ASSERT(!TC_isInit_);
}

void TChemReactor::SetDensity(double rho) {
  temperature_and_mass_fractions_[0] = rho;
  double T = 0.0;
  TC_getTmixMs(temperature_and_mass_fractions_.data(),
               temperature_and_mass_fractions_.size(), &T);
  SetTemperature(T);
}

void TChemReactor::SetPressure(double p) {
  pressure_ = p;
  TC_setThermoPres(pressure_);
}

void TChemReactor::SetTemperature(double T) {
  temperature_ = T;
  temperature_and_mass_fractions_[0] = T;
  TC_getRhoMixMs(temperature_and_mass_fractions_.data(),
                 temperature_and_mass_fractions_.size(), &density_);
}

double TChemReactor::GetInternalEnergy() const {
  span<const double> mass_fractions = GetMassFractions();
  return std::inner_product(mass_fractions.begin(), mass_fractions.end(),
                            internal_energies_.begin(), 0.0);
}

void TChemReactor::SetInternalEnergy(double U, double tol) {
  const double T = GetTemperatureFromInternalEnergy_(
      U, temperature_and_mass_fractions_, internal_energies_, tol);
  SetTemperature(T);
}

void TChemReactor::SetMassFractions(span<const double> Y) {
  const double total = std::accumulate(Y.begin(), Y.end(), 0.0);
  assert(total > 0.0);
  std::transform(Y.begin(), Y.end(),
                 temperature_and_mass_fractions_.begin() + 1,
                 [=](double Yi) { return Yi / total; });
  TC_getRhoMixMs(temperature_and_mass_fractions_.data(),
                 temperature_and_mass_fractions_.size(), &density_);
  UpdateMoleFractionsFromMassFractions(moles_, GetMassFractions());
}

void TChemReactor::SetMoleFractions(span<const double> X) {
  const double total = std::accumulate(X.begin(), X.end(), 0.0);
  assert(total > 0.0);
  std::transform(X.begin(), X.end(),
                 temperature_and_mass_fractions_.begin() + 1,
                 [=](double Yi) { return Yi / total; });
  TC_getRhoMixMl(temperature_and_mass_fractions_.data(),
                 temperature_and_mass_fractions_.size(), &density_);
  std::copy(temperature_and_mass_fractions_.begin() + 1,
            temperature_and_mass_fractions_.end(), moles_.begin());
  span<double> masses = make_span(temperature_and_mass_fractions_).subspan(1);
  UpdateMassFractionsFromMoleFractions(masses, moles_);
}

double TChemReactor::GetCp() const {
  span<const double> cps = GetCps();
  span<const double> Y = GetMassFractions();
  return std::inner_product(Y.begin(), Y.end(), cps.begin(), 0.0);
}

double TChemReactor::GetCv() const {
  span<const double> cvs = GetCvs();
  span<const double> Y = GetMassFractions();
  return std::inner_product(Y.begin(), Y.end(), cvs.begin(), 0.0);
}

void TChemReactor::UpdateThermoState() {
  if (thermo_temperature_ != temperature_) {
    TC_setThermoPres(pressure_);
    TC_getCpSpecMs(temperature_, GetNSpecies(), cps_.data());
    TC_getCvSpecMs(temperature_, GetNSpecies(), cvs_.data());
    TC_getUspecMs(temperature_, GetNSpecies(), internal_energies_.data());
    TC_getHspecMs(temperature_, GetNSpecies(), enthalpies_.data());
    thermo_temperature_ = temperature_;
  }
}

double TChemReactor::GetGamma() const {
  const double cp = GetCp();
  const double cv = GetCv();
  return cp / cv;
}

double TChemReactor::GetSpeedOfSound() const {
  const double gamma = GetGamma();
  return std::sqrt(gamma * pressure_ / density_);
}

int TChemReactor::GetNSpecies() const { return TC_getNspec(); }

void TChemReactor::Advance(double dt) {
  auto source_term = [](span<double> dydt, span<const double> y, double t) {
    TC_getSrcCV(const_cast<double*>(y.data()), y.size(), dydt.data());
  };
  auto jacobian = [](span<double> jac, span<const double> TandYs, double t) {
    int Ns = TandYs.size() - 1;
    TC_getJacCVTYNanl(const_cast<double*>(TandYs.data()), Ns, jac.data());
  };
  TC_setDens(density_);
  ode_solver_->Integrate(source_term, temperature_and_mass_fractions_, 0, dt,
                         jacobian);
  UpdateMoleFractionsFromMassFractions(moles_, GetMassFractions());
  pressure_ = TC_getThermoPres();
  SetTemperature(temperature_and_mass_fractions_[0]);
  UpdateThermoState();
}

} // namespace ideal_gas
} // namespace fub