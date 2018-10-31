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

#ifndef FUB_SOLVER_EULER_REACTOR_HPP
#define FUB_SOLVER_EULER_REACTOR_HPP

#include "fub/core/function_ref.hpp"
#include "fub/core/span.hpp"

#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace fub {

/// This class encapsulates all exceptions which can occur in the
/// FlameMasterReactor class
class FlameMasterReactorException : public std::runtime_error {
public:
  FlameMasterReactorException(const char* message)
      : std::runtime_error(message) {}

  FlameMasterReactorException(std::string message)
      : std::runtime_error(std::move(message)) {}
};

/// This abstract base class encapsulates the underlying chemistry for the
/// FlameMasterReactor.
struct FlameMasterMechanism {

  virtual ~FlameMasterMechanism() = default;

  virtual std::unique_ptr<FlameMasterMechanism> Clone() const = 0;

  virtual void ComputeProductionRates(span<double> cdot, span<double> w,
                                      span<double> k, span<double> c,
                                      span<double> M, double temp,
                                      double pressure) const = 0;

  virtual void ComputeThermoData(span<double> h, span<double> cp, double t,
                                 span<double> s) const = 0;

  virtual int GetNSpecies() const = 0;

  virtual int GetNReactions() const = 0;

  virtual int GetNThirdBodyReactions() const = 0;

  virtual std::vector<std::string> GetSpeciesNames() const = 0;

  virtual void GetMolarMass(span<double>) const = 0;

  virtual double GetUniversalGasConstant() const { return 8314.462175; }

  virtual int GetNSpecs() const { return GetNSpecies(); }
};

struct FlameMasterState {
  /// The number of species in the mechanism
  int nSpecies;
  /// Effective number of species for mechanisms with steady state species
  int nSpeciesEffective;
  /// The number of reactions in the mechanism
  int nReactions;
  /// The number of third body reactions in the mechanism
  int nThirdBodyReactions;
  /// Universal Gas Constant for this mechanism.
  double gasConstant;

  std::vector<double> massFractions;
  std::vector<double> molesStorage;

  /// Derived from massFractions, stores the actual mole counts during time
  /// Advancement
  span<double> moles;

  /// We make this a pointer because we want to make sure it is stored at the
  /// beginning of moles (for CVode)
  double* temperature;

  double density;

  /// Computational space for the reaction mechanism
  std::vector<double> production_rates;
  std::vector<double> reaction_rates;
  std::vector<double> rate_coefficients;
  std::vector<double> third_body_concentrations;
  std::vector<double> enthalpies;
  std::vector<double> heat_capacities_at_constant_pressure;
  std::vector<double> entropies;

  /// We store the names of the species here
  std::vector<std::string> speciesNames;

  /// We store the molar masses of the species here
  std::vector<double> molarMasses;

  /// Absolute integration tolerance
  double abstol;
  /// Relative integration tolerance
  double reltol;

  /// A vector containing temperature for setPressureIsentropic()
  std::array<double, 2> setPVector;

  /// Temperature at which the thermodynamic state was last evaluated
  double thermoTemp;
};

/// \brief A class mimicking the IdealGasMix / Reactor / ReactorNet interface of
/// Cantera, but with FlameMaster chemistry.
class FlameMasterReactor {
public:
  explicit FlameMasterReactor(std::unique_ptr<FlameMasterMechanism> mechanism);
  ~FlameMasterReactor() = default;

  FlameMasterReactor(const FlameMasterReactor& other);
  FlameMasterReactor(FlameMasterReactor&& other) noexcept;

  /// \brief Advance the reactor in time by dt and call a function for each
  /// internal timestep
  void Advance(double dt,
               function_ref<int(double, FlameMasterReactor*)> feedbackFun);

  /// \brief Advance the reactor in time by dt.
  ///
  /// \param[in] dt  The time step size.
  ///
  /// \throw FlameMasterReactorException  This exception may be thrown if the
  /// ode solver could not converge to a solution.
  void Advance(double dt);

  /// \brief Advance the reactor by one internal time step and return the time
  /// step size
  double Step();

  /// \brief Advance the reactor in time by dt and return the time where
  /// \f$dT/dt\f$ peaked.
  double AdvanceAndFindMaxdT(double dt);

  /// Set tolerances for the integration. Defaults are 1e-9 for reltol and
  /// 1e-15 for abstol.
  // void SetTolerances(double reltol, double abstol);

  ///@name Thermodynamic properties
  ///@{

  /// Returns the pressure of the mixture
  ///
  /// This used the ideal gas law to derive the pressure from
  /// temperature, the mass fractions and density.
  ///
  /// The unit for pressure is Pascal.
  double GetPressure() const {
    return GetDensity() * GetUniversalGasConstant() / GetMeanMolarMass() *
           GetTemperature();
  }

  /**
   * Set the pressure of the mixture
   *
   * This actually sets the density based on the current
   * temperature and mass fractions, using the ideal gas law
   *
   * The unit for pressure is Pascal.
   */
  void SetPressure(double pressure) {
    SetDensity(pressure * GetMeanMolarMass() / GetTemperature() /
               GetUniversalGasConstant());
  }

  /**
   * Adjust the pressure of the mixture isentropically
   *
   * Flamemaster does not offer entropy data. However, we still
   * have the relation dH = V dp + T dS, which reduces to
   * dH = V dp in the isentropic case. This function integrates
   * H to the desired p and ensures fixed entropy this way.
   *
   * The function simultaneously integrates velocity using the relation
   * du = - 1/(c \rho) dp = - sqrt(1/(\gamma p \rho)) dp.
   * This relation is fulfiled across rarefaction waves, and the call
   * allows to use the function to calculate the state behind such a
   * wave.
   */
  double SetPressureIsentropic(double pressure);

  /// Returns the density of the current mixture
  ///
  /// The units are \f$kg/m^3\f$
  double GetDensity() const { return state_.density; }

  /**
   * Sets the density of the current mixture
   *
   * The units are \f$kg/m^3\f$
   */
  void SetDensity(double density) { state_.density = density; }

  /**
   * Returns the temperature of the current mixture
   *
   * In Kelvin
   */
  double GetTemperature() const { return *(state_.temperature); }

  /**
   * Set the temperature of the current mixture
   *
   * In Kelvin
   */
  void SetTemperature(double temperature);

  /**
   * Return the specific heat capacity cv
   */
  double GetCv() const;

  /**
   * Return the specific heat capacity cp
   */
  double GetCp() const;

  /**
   * Return the specific heat capacities cp for all species
   */
  span<const double> GetCps() const;

  /**
   * Return the specific entropy
   */
  double GetEntropy() const;

  ///@}

  ///@name Mass-/Mole fractions
  ///@{

  /**
   * Return the name of the ith species
   */
  const std::string& GetSpeciesName(int i) const {
    return state_.speciesNames[i];
  }

  /**
   * Return the number of species in the mechanism
   */
  int GetNSpecies() const { return state_.nSpecies; }

  /**
   * Average a quantity over the mole fractions
   */
  double MeanX(span<const double> quantity) const;

  /**
   * Average a quantity over the mass fractions
   */
  double MeanY(span<const double> quantity) const;

  ///@{
  /**
   * Set the dimensionless mass fractions of the mixture.
   *
   * You can either supply them as a double* array, or as
   * a string of the form
   *   species: fraction, species: fraction, ..
   *
   */
  void SetMassFractions(span<const double> newMassFractions);
  void SetMassFractions(std::string massFractions);
  ///@}

  ///@{
  /**
   * Set the dimensionless mole fractions of the mixture
   *
   * You can either supply them as a double* array, or as
   * a string of the form
   *   species: fraction, species: fraction, ..
   */
  void SetMoleFractions(span<const double> newMoleFractions);
  void SetMoleFractions(std::string newMoleFractions);
  ///@}

  /**
   * Return the mass fractions of the species as a double* array
   */
  span<const double> GetMassFractions() const { return state_.massFractions; }

  /**
   * Get the molar mass of a single species
   *
   * Units are \f$kg/kmol\f$
   */
  const double GetSpeciesMolarMass(size_t i) const {
    return state_.molarMasses[i];
  }

  /**
   * Get the molar masses for all species
   */
  span<const double> GetMolarMasses() const { return state_.moles; }

  /**
   * Get the overall molar mass
   *
   * Units are \f$kg/kmol\f$
   */
  double GetMeanMolarMass() const;

  /**
   * Get the species' enthalpies
   */
  span<const double> GetEnthalpies();

  /**
   * Return the mix'es specific enthalpy
   *
   * The unit is \f$J/kg\f$
   */
  double GetEnthalpy() const;

  /**
   * Set the mix'es specific enthalpy
   *
   * The unit is \f$J/kg\f$. This function actually tries to match the
   * mixtures temperature such that the enthalpy is correct. It does so by
   * using the same Newton algorithm which is also used by Cantera.
   */
  void SetEnthalpy(double enthalpy, double dTtol);
  ///@}

  /**
   * Return the mix'es specific internal energy
   *
   * The unit is \f$J/kg\f$
   */
  double GetInternalEnergy() const;

  /**
   * Set the mix'es specific internal energy
   *
   * The unit is \f$J/kg\f$. This function actually tries to match the
   * mixtures temperature such that the energy is correct. It does so by
   * using the same Newton algorithm which is also used by Cantera.
   */
  void SetInternalEnergy(double energy, double dTtol = 1E-6);
  ///@}

  /// \brief Get the current reaction rates d/dt m, where m denotes the actual
  /// mole counts: Y = m * molarMasses / density
  span<const double> GetReactionRates() const;

  span<const double> GetProductionRates();

  /// \brief Retrieve the universal gas constant.
  ///
  /// Note that some mechanisms *can* override this, to use  normalized
  /// quantities!
  double GetUniversalGasConstant() const { return state_.gasConstant; }

private:
  std::unique_ptr<FlameMasterMechanism> mechanism_;
  FlameMasterState state_;
};

} // namespace fub

#endif
