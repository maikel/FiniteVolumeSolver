// Copyright (c) 2021 Maikel Nadolski
// Copyright (c) 2021 Christian Zenker
// Copyright (c) 2020 Rupert Klein
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

#ifndef FUB_EQUATIONS_PERFECT_GAS_MIX_PLENA_CONTROL_HPP
#define FUB_EQUATIONS_PERFECT_GAS_MIX_PLENA_CONTROL_HPP

/// \file
///
/// This file contains functions to control the states in the compressor and
/// turbine plena. They are relevant in studies regarding the efficieny of the
/// total gas turbine machine.

#include "fub/equations/perfect_gas_mix/ArrheniusKinetics.hpp"

namespace fub::perfect_gas_mix {

/// \brief Compute the ignition delay time for given mass fractions and
/// temperature
template <int Rank>
double IgnitionDelay(const PerfectGasMix<Rank>& equation,
                     const ArrheniusKineticsOptions& kinetics, double Y,
                     double T) {
  const double gm1 = equation.gamma - 1.0;
  double TT = T / (gm1 * kinetics.Q * std::max(Y, eps));
  double Ea0 = kinetics.EA * (1.0 / T);
  double B0 = kinetics.B * exp(kinetics.EA * (1.0 - (1.0 / T)));
  return ((TT / (B0 * Ea0)) * (1.0 + (2.0 + TT) / Ea0));
};

// /*
// ==================================================================================
// */

/// \brief Compute the derivative of the function IgnitionDelay.
///
/// This function is needed in the newton raphson iteration method.
///
/// \see TemperatureForIgnitionDelay
template <int Rank>
double dIgnitionDelay_dT(const PerfectGasMix<Rank>& equation,
                         const ArrheniusKineticsOptions& kinetics, double Y,
                         double T) {
  const double gm1 = equation.gamma - 1.0;
  double TQ = T / (gm1 * kinetics.Q * std::max(Y, eps));
  double TE = T / kinetics.EA;
  double A = 1.0 / (kinetics.B * kinetics.EA * gm1 * kinetics.Q * Y);
  double dtidT = A * kinetics.EA *
                 ((-1.0 + 2.0 * TE) * (1.0 + TE * (2.0 + TQ)) +
                  2.0 * TE * TE * (1.0 + TQ)) *
                 exp(-kinetics.EA * (1.0 - 1.0 / T));
  return dtidT;
};

// /*
// ==================================================================================
// */

/// \brief Find the temperature T for given mass fractions Y and ignition delay
/// time t_ign such that `t_ign = IgnitionDelay(Y, T)`.
///
template <int Rank>
double TemperatureForIgnitionDelay(const PerfectGasMix<Rank>& equation,
                                   const ArrheniusKineticsOptions& kinetics,
                                   double tign, double Y, double T0) {
  double T = T0;
  double ti = IgnitionDelay(equation, kinetics, Y, T);
  double dti = ti - tign;
  while (std::abs(dti) > eps) {
    T -= dti / dIgnitionDelay_dT(equation, kinetics, Y, T);
    ti = IgnitionDelay(Y, T);
    dti = ti - tign;
  }
  return T;
};

namespace gt {
/// \namespace
/// \brief This namespace includes all classes and functions related to the
/// efficiency computation for total gas turbine machine.

/// Compressor and turbine plena are represented by averaged pressure and
/// temperature values.
struct PlenumState {
  double pressure;
  double temperature;
};

template <int PlenumRank> class Control {
public:
  // Constructors

  /// Notify the control with a new turbine plenum value and update its current
  /// performance.
  ///
  /// The formula is given in GT_Notes equation (35).
  ///
  /// \return Returns the new performance value
  double
  UpdateTurbinePlenum(double mdot,
                      const PlenumState& turbine_plenum_current) noexcept {
    PlenumState T_half;
    T_half.pressure =
        0.5 * (turbine_plenum_old_.pressure + turbine_plenum_current.pressure);
    T_half.temperature = 0.5 * (turbine_plenum_old_.temperature +
                                turbine_plenum_current.temperature);
    const double gmoog = equation_.gamma_minus_one_over_gamma;
    const double cp = equation_.heat_capacity_at_constant_pressure;
    performance_turbine_ = mdot * cp * efficiency_turbine_ *
                           T_half.temperature *
                           (1 - std::pow(T_half.pressure, gm1og));
    turbine_plenum_old_ = turbine_plenum_current;
    return performance_turbine_;
  }

  /// The formula is given in GT_Notes equation (20).
  ///
  /// \return Returns the compressure pressure ratio
  double CompressorPressureRatio(double rpm) const noexcept {
    double pratio = pratiomean_ + c_0_ * pratiovar_ *
                                      std::atan(M_PI * (rpm / rpmmean_ - 1.0)) /
                                      M_PI;
    return pratio;
  }

  double UpdateCompressorPlenum(double flux_rho, Duration dt) {
    double rho_old = compressor_plenum_.pressure / compressor_plenum_.temperature / equation_.Rspec;
    UpdateCompressorFromRPM(current_rpm_);
    double rho_new = compressor_plenum_.pressure / compressor_plenum_.temperature / equation_.Rspec;
    double mdot_compressor = flux_rho * surface_area_tube_inlet_ + volume_compressor_plenum_ / dt.count() * (rho_new - rho_old);
  }

private:
  PerfectGasMix<PlenumRank> equation_;
  PlenumState turbine_plenum_old_;
  PlenumState compressor_plenum_; /// the 
  double efficiency_turbine_{0.9}; ///< efficiency from the turbine
  
  double surface_area_tube_inlet_{1.0};
  double surface_area_compressor_to_compressor_plenum_{8.0 * surface_area_tube_inlet_};
  double volume_compressor_plenum_{20.0 * surface_area_tube_inlet_ * length_tube_};
  double surface_area_turbine_plenum_to_turbine_{4.0 * surface_area_tube_inlet_};

  double performance_turbine_;
  double rpmmax_ = 60000.0 / 60.0;
  double rpmmin_ = 10000.0 / 60.0;
  double rpmmean_ = 0.5 * (rpmmax_ + rpmmin_);
  double pratiomax_ = 50.0;
  double pratiomin_ = 1.0;
  double pratiomean_ = 0.5 * (pratiomax_ + pratiomin_);
  double pratiovar_ = pratiomax_ - pratiomin_;
  double c_0_ = 22.0 / 17.0;

  double current_rpm_{rpmmin_};
};

} // namespace gt

/// \brief Compute the performance for the turbine, P_T.
double PerformanceTurbine(const PerfectGasMix<Rank>& equation, double mdot)

} // namespace fub::perfect_gas_mix

#endif