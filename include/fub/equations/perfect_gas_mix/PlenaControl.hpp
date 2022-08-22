// Copyright (c) 2021 Maikel Nadolski
// Copyright (c) 2021 Christian Zenker
// Copyright (c) 2021 Rupert Klein
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

#include "fub/equations/PerfectGasMix.hpp"
#include "fub/equations/perfect_gas_mix/ArrheniusKinetics.hpp"

#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include <range/v3/view/zip.hpp>

#include <cmath>

namespace fub::perfect_gas_mix {

/// \brief Compute the ignition delay time for given mass fractions and
/// temperature
template <int Rank>
double IgnitionDelay(const PerfectGasMix<Rank>& equation,
                     const ArrheniusKineticsOptions& kinetics, double Y,
                     double T, double eps) {
  const double gm1 = equation.gamma - 1.0;
  double TT = T / (gm1 * kinetics.Q * std::max(Y, eps));
  double Ea0 = kinetics.EA * (1.0 / T);
  double B0 = kinetics.B * exp(kinetics.EA * (1.0 - (1.0 / T)));
  return ((TT / (B0 * Ea0)) * (1.0 + (2.0 + TT) / Ea0));
}

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
                         double T, double eps) {
  const double gm1 = equation.gamma - 1.0;
  double TQ = T / (gm1 * kinetics.Q * std::max(Y, eps));
  double TE = T / kinetics.EA;
  double A = 1.0 / (kinetics.B * kinetics.EA * gm1 * kinetics.Q * Y);
  double dtidT = A * kinetics.EA *
                 ((-1.0 + 2.0 * TE) * (1.0 + TE * (2.0 + TQ)) +
                  2.0 * TE * TE * (1.0 + TQ)) *
                 exp(-kinetics.EA * (1.0 - 1.0 / T));
  return dtidT;
}

// /*
// ==================================================================================
// */

/// \brief Find the temperature T for given mass fractions Y and ignition delay
/// time t_ign such that `t_ign = IgnitionDelay(Y, T)`.
///
template <int Rank>
double TemperatureForIgnitionDelay(const PerfectGasMix<Rank>& equation,
                                   const ArrheniusKineticsOptions& kinetics,
                                   double tign, double Y, double T0,
                                   double eps) {
  double T = T0;
  double ti = IgnitionDelay(equation, kinetics, Y, T, eps);
  double dti = ti - tign;
  while (std::abs(dti) > eps) {
    T -= dti / dIgnitionDelay_dT(equation, kinetics, Y, T, eps);
    ti = IgnitionDelay(equation, kinetics, Y, T, eps);
    dti = ti - tign;
  }
  return T;
}

namespace gt {
/// \namespace
/// \brief This namespace includes all classes and functions related to the
/// efficiency computation for total gas turbine machine.

/// Compressor and turbine plena are represented by averaged pressure and
/// temperature values.
struct PlenumState {
  double pressure{};
  double temperature{};
  double power{};
  double mass_flow_in{};
  double mass_flow_out{};
  bool SEC_Mode{false};
};

struct ControlOptions {
  ControlOptions() = default;
  ControlOptions(const ProgramOptions& options);

  void Print(SeverityLogger& log);

  double efficiency_turbine{0.9}; ///< efficiency from the turbine

  double length_tube{1.0};              ///< the length of the tubes
  double surface_area_tube_inlet{1.0};  ///< initial surface area of the tube
  double surface_area_tube_outlet{4.0}; ///< final surface area of the tube
  /// surface area from the compressor to the compressor plenum
  double surface_area_compressor_to_compressor_plenum{8.0 *
                                                      surface_area_tube_inlet};
  /// volume of the compressor plenum
  double volume_compressor_plenum{20.0 * surface_area_tube_inlet * length_tube};
  /// surface area from the turbine plenum to the turbine
  double surface_area_turbine_plenum_to_turbine{4.0 * surface_area_tube_inlet};
  /// volume of the turbine plenum
  double volume_turbine_plenum{20.0 * surface_area_tube_inlet * length_tube};
  double efficiency_compressor{0.9};        ///< efficiency from the compressor
  double rpmmin = 10000.0 / 60.0;           ///< minimum speed of the compressor
  double rpmmax = 60000.0 / 60.0;           ///< maximum speed of the compressor
  double rpmmean = 0.5 * (rpmmax + rpmmin); ///< mean speed of the compressor
  double pratiomin = 1.0;  ///< minimum pressure from the compressor
  double pratiomax = 50.0; ///< maximum pressure from the compressor
  /// mean pressure from the compressor
  double pratiomean = 0.5 * (pratiomax + pratiomin);
  /// maximum pressure difference from the compressor
  double pratiovar = pratiomax - pratiomin;

  /// constant to compute the pressure from the speed of the compressor (see
  /// function: CompressorPressureRatio)
  double c_0 = 22.0 / 17.0;

  double inertial_moment{0.1}; ///< the initial moment of inertia of the rotor
  double mu{10.0};             ///< some constant for function ChangeRPM

  double target_pressure_compressor{
      6.5}; ///< the target pressure from the compressor

  /// tolerance for the relative difference in pressure to acitvate the SEC Mode
  double pressure_relDiffTolerance_SEC{1e-04};

  double Q{3.0 / 0.4}; ///< Needs to be equal to the ArrheniusKinetics Q.

  PlenumState initial_turbine_state{};
};

struct ControlState {
  double current_rpm{};
  PlenumState compressor{};
  PlenumState turbine{};
  double power_out{};
  double fuel_consumption{};
  double fuel_consumption_rate{};
  double efficiency{};
};

class Control {
public:
  // Constructors

  Control(const PerfectGasConstants& eq, const ControlOptions& options);

  // Modifier

  void UpdatePlena(double mdot_turbine,
                   const PlenumState& turbine_boundary_state, double flux_rho,
                   double flux_spec, double flux_rho_last, Duration dt);

  // Accessors

  std::shared_ptr<const ControlState> GetSharedState() const noexcept {
    return state_;
  }

  std::shared_ptr<ControlState> GetSharedState() noexcept { return state_; }

  const ControlOptions& GetOptions() const noexcept { return options_; }

private:
  PerfectGasConstants equation_;
  ControlOptions options_{};
  std::shared_ptr<ControlState> state_;
};

} // namespace gt

} // namespace fub::perfect_gas_mix

namespace boost::serialization {

template <typename Archive>
void serialize(Archive& ar, ::fub::perfect_gas_mix::gt::PlenumState& state,
               unsigned int) {
  // clang-format off
  ar & state.pressure;
  ar & state.temperature;
  ar & state.power;
  ar & state.mass_flow_in;
  ar & state.mass_flow_out;
  ar & state.SEC_Mode;
  // clang-format on
}

template <typename Archive>
void serialize(Archive& ar, ::fub::perfect_gas_mix::gt::ControlState& state,
               unsigned int) {
  // clang-format off
  ar & state.current_rpm;
  ar & state.compressor;
  ar & state.turbine;
  ar & state.power_out;
  ar & state.fuel_consumption;
  ar & state.fuel_consumption_rate;
  ar & state.efficiency;
  // clang-format on
}

} // namespace boost::serialization

#endif