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

#ifndef FUB_IDEAL_GAS_TCHEM_REACTOR_HPP
#define FUB_IDEAL_GAS_TCHEM_REACTOR_HPP

#include "fub/core/span.hpp"
#include "fub/ode_solver/OdeSolver.hpp"

#include <memory>
#include <stdexcept>
#include <vector>

namespace fub {
namespace ideal_gas {

struct TChemNotInitialized : std::runtime_error {
  TChemNotInitialized()
      : std::runtime_error(
            "You have not initialized TChem before using this reactor.") {}
};

struct TChemReactorUniqueError : std::runtime_error {
  TChemReactorUniqueError()
      : std::runtime_error(
            "Only one Reactor at the time is allowed to exist.") {}
};

struct TChemMechanism {
  virtual ~TChemMechanism() = default;
  virtual int getNSpecies() const = 0;
  virtual const char* GetMechanismResource() const = 0;
  virtual const char* GetThermoResource() const = 0;
};

class TChemReactor {
public:
  using Feedback = OdeSolver::feedback_type;

  TChemReactor(const TChemMechanism& mechanism);
  ~TChemReactor() noexcept;

  TChemReactor(const TChemReactor&) = delete;
  TChemReactor& operator=(const TChemReactor&) = delete;
  TChemReactor(TChemReactor&&) = delete;
  TChemReactor& operator=(TChemReactor&&) = delete;

  void SetOdeSolver(std::unique_ptr<OdeSolver> solver);

  void Advance(double dt);
  void Advance(double dt, Feedback feedback);

  int GetNSpecies() const;
  double GetPressure() const { return pressure_; }
  double GetDensity() const { return density_; }
  double GetTemperature() const { return temperature_; }
  span<const double> GetMassFractions() const {
    return make_span(temperature_and_mass_fractions_).subspan(1);
  }
  span<const double> GetMoleFractions() const { return moles_; }
  span<const double> GetCps() const { return cps_; }
  span<const double> GetCvs() const { return cvs_; }
  span<const double> GetEnthalpies() const { return enthalpies_; }
  span<const double> GetInternalEnergies() const { return internal_energies_; }
  const OdeSolver& GetOdeSolver() const { return *ode_solver_; }

  double GetSpeedOfSound() const;
  double GetInternalEnergy() const;
  double GetCp() const;
  double GetCv() const; 
  double GetGamma() const;

  void SetPressure(double p);
  void SetDensity(double rho);
  void SetTemperature(double T);
  void SetInternalEnergy(double e, double tolerance = 1e-6);
  void SetMassFractions(span<const double> Y);
  void SetMoleFractions(span<const double> X);

  void UpdateThermoState();

  const char* GetSpeciesName(int i) const;

private:
  std::vector<double> temperature_and_mass_fractions_;
  std::vector<double> moles_;
  std::vector<double> internal_energies_;
  std::vector<double> enthalpies_;
  std::vector<double> cps_;
  std::vector<double> cvs_;
  double density_;
  double pressure_;
  double temperature_;
  double thermo_temperature_;
  std::unique_ptr<OdeSolver> ode_solver_;
};

} // namespace ideal_gas
} // namespace fub

#endif