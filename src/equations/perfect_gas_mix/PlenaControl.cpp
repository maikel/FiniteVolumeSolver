// Copyright (c) 2021 Maikel Nadolski
// Copyright (c) 2021 Christian Zenker
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

#include "fub/equations/perfect_gas_mix/PlenaControl.hpp"

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/map.hpp>

namespace fub::perfect_gas_mix::gt {

namespace {
/// Notify the control with a new turbine plenum value and update its current
/// power.
///
/// The formula is given in GT_Notes equation (35).
///
/// \return Returns the new power value for the turbine
void UpdateTurbinePlenum(ControlState& state, const PerfectGasConstants& eq,
                         const ControlOptions& options, double mdot,
                         const PlenumState& turbine_plenum_current,
                         double flux_rho_last, Duration dt) noexcept {
  const double pressure_half =
      0.5 * (state.turbine.pressure + turbine_plenum_current.pressure);
  const double T_half =
      0.5 * (state.turbine.temperature + turbine_plenum_current.temperature);
  const double gmoog = eq.gamma_minus_one_over_gamma;
  const double cp = eq.heat_capacity_at_constant_pressure;
  const double ootau = 1.0;
  double dt_over_tau = dt.count() * ootau;
  state.turbine.mass_flow_in =
      (1.0 - dt_over_tau) * state.turbine.mass_flow_in +
      dt_over_tau * flux_rho_last * options.surface_area_tube_outlet;
  state.turbine.mass_flow_out =
      (1.0 - dt_over_tau) * state.turbine.mass_flow_out + dt_over_tau * mdot;
  state.turbine.power = mdot * cp * options.efficiency_turbine * T_half *
                        (1.0 - std::pow(pressure_half, -gmoog));
  state.turbine.pressure = turbine_plenum_current.pressure;
  state.turbine.temperature = turbine_plenum_current.temperature;
}

/// The formula is given in GT_Notes equation (20).
///
/// \return Returns the compressure pressure ratio
double CompressorPressureRatio(const ControlOptions& options,
                               double rpm) noexcept {
  double pratio = options.pratiomean +
                  options.c_0 * options.pratiovar *
                      std::atan(M_PI * (rpm / options.rpmmean - 1.0)) / M_PI;
  return pratio;
}

/// \brief Compute the compressor state by given speed of the compressor.
void UpdateCompressorFromRPM(ControlState& state, const PerfectGasConstants& eq,
                             const ControlOptions& options) {
  double new_pressure = CompressorPressureRatio(options, state.current_rpm);
  state.compressor.pressure = new_pressure;
  const double gmoog = eq.gamma_minus_one_over_gamma;

  // see GT_Notes equation 21
  state.compressor.temperature = 1.0 + (std::pow(new_pressure, gmoog) - 1.0) /
                                           options.efficiency_compressor;
}

/// Notify the control with a new compressor plenum value and update its
/// current power.
///
/// The formula is given in GT_Notes equation (24).
///
/// \return Returns the new power value
void UpdateCompressorPlenum(ControlState& state, const PerfectGasConstants& eq,
                            const ControlOptions& options, double flux_rho,
                            double flux_spec, Duration dt) noexcept {
  double temperature_old = state.compressor.temperature;
  double rho_old =
      state.compressor.pressure / state.compressor.temperature * eq.ooRspec;

  // update the compressor_plenum state with the new speed of the compressor
  // so we get new pressure and temperature
  UpdateCompressorFromRPM(state, eq, options);

  double rho_new =
      state.compressor.pressure / state.compressor.temperature * eq.ooRspec;

  const double mdot = flux_rho * options.surface_area_tube_inlet +
                      options.volume_compressor_plenum / dt.count() *
                          (rho_new - rho_old); // equation (23) GT_Notes

  double T_v_half = 0.5 * (temperature_old + state.compressor.temperature);
  double T_ref = 1.0;
  const double cp = eq.heat_capacity_at_constant_pressure;

  double power_compressor_increase =
      mdot * (T_v_half - T_ref) * cp / options.efficiency_compressor;

  const double ootau = 1.0;
  double dt_over_tau = dt.count() * ootau;

  state.compressor.mass_flow_in =
      (1.0 - dt_over_tau) * state.compressor.mass_flow_in + dt_over_tau * mdot;

  state.compressor.mass_flow_out =
      (1.0 - dt_over_tau) * state.compressor.mass_flow_out +
      dt_over_tau * flux_rho * options.surface_area_tube_inlet;

  state.compressor.power = (1.0 - dt_over_tau) * state.compressor.power +
                           dt_over_tau * power_compressor_increase;
  state.fuel_consumption +=
      flux_spec * options.surface_area_tube_inlet * dt.count();

  state.fuel_consumption_rate =
      (1.0 - dt_over_tau) * state.fuel_consumption_rate +
      dt_over_tau * flux_spec * options.surface_area_tube_inlet;
}

/// \brief Compute the new rotation speed of the compressor.
///
/// see GT_notes equation (42 - 44)
void ChangeRPM(ControlState& state, const ControlOptions& options,
               Duration dt) noexcept {
  double pressure = state.compressor.pressure;
  const double pressure_relative_difference =
      (options.target_pressure_compressor - pressure) /
      options.target_pressure_compressor;
  const double exp_term = std::exp(-options.mu * pressure_relative_difference);
  const double comp_rate = 0.015 * (1.0 - exp_term);
  constexpr double pi2 = M_PI * M_PI;
  const double Ieff = 2.0 * pi2 * options.inertial_moment;

  const double power_netto = state.turbine.power - state.compressor.power;
  const double d_Erot_dt = (pressure < options.target_pressure_compressor)
                               ? comp_rate * power_netto
                               : 0.0;
  const double Erot =
      Ieff * state.current_rpm * state.current_rpm + d_Erot_dt * dt.count();
  const double rpm_new =
      std::clamp(std::sqrt(Erot / Ieff), options.rpmmin, options.rpmmax);
  // FUB_ASSERT(options.rpmmin <= rpm_new && options.rpmmax >= rpm_new);
  state.current_rpm = rpm_new;

  // values to track:
  state.power_out = power_netto - d_Erot_dt; // Power output

  // efficiency of the machine
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  state.efficiency = std::max(0.0, state.power_out) /
                     (options.Q * state.fuel_consumption_rate + eps);

  if (!state.compressor.SEC_Mode) {
    // compare the pressure relative difference with a given tolerance
    state.compressor.SEC_Mode =
        pressure_relative_difference < options.pressure_relDiffTolerance_SEC;
  } else {
    double fall_back_pressure = 0.9 * options.target_pressure_compressor;
    state.compressor.SEC_Mode = pressure > fall_back_pressure;
  }
}
} // namespace

#ifdef FUB_GT_CONTROL_GET_OPTION
#undef FUB_GT_CONTROL_GET_OPTION
#endif
#define FUB_GT_CONTROL_GET_OPTION(var) var = GetOptionOr(options, #var, var)
ControlOptions::ControlOptions(const ProgramOptions& options) {
  FUB_GT_CONTROL_GET_OPTION(efficiency_turbine);
  FUB_GT_CONTROL_GET_OPTION(length_tube);
  FUB_GT_CONTROL_GET_OPTION(surface_area_tube_inlet);
  FUB_GT_CONTROL_GET_OPTION(surface_area_tube_outlet);
  FUB_GT_CONTROL_GET_OPTION(surface_area_compressor_to_compressor_plenum);
  FUB_GT_CONTROL_GET_OPTION(volume_compressor_plenum);
  FUB_GT_CONTROL_GET_OPTION(surface_area_turbine_plenum_to_turbine);
  FUB_GT_CONTROL_GET_OPTION(volume_turbine_plenum);
  FUB_GT_CONTROL_GET_OPTION(efficiency_compressor);
  FUB_GT_CONTROL_GET_OPTION(rpmmin);
  FUB_GT_CONTROL_GET_OPTION(rpmmax);
  FUB_GT_CONTROL_GET_OPTION(pratiomin);
  FUB_GT_CONTROL_GET_OPTION(pratiomax);
  FUB_GT_CONTROL_GET_OPTION(c_0);
  FUB_GT_CONTROL_GET_OPTION(inertial_moment);
  FUB_GT_CONTROL_GET_OPTION(mu);
  FUB_GT_CONTROL_GET_OPTION(target_pressure_compressor);
  FUB_GT_CONTROL_GET_OPTION(pressure_relDiffTolerance_SEC);
  FUB_GT_CONTROL_GET_OPTION(Q);
  initial_turbine_state.pressure = GetOptionOr(
      options, "initial_turbine_pressure", initial_turbine_state.pressure);
  initial_turbine_state.temperature =
      GetOptionOr(options, "initial_turbine_temperature",
                  initial_turbine_state.temperature);

  rpmmean = 0.5 * (rpmmax + rpmmin);
  pratiomean = 0.5 * (pratiomax + pratiomin);
  pratiovar = pratiomax - pratiomin;
}
#undef FUB_GT_CONTROL_GET_OPTION

#ifdef FUB_GT_CONTROL_PRINT
#undef FUB_GT_CONTROL_PRINT
#endif
#define FUB_GT_CONTROL_PRINT(variable)                                         \
  BOOST_LOG(log) << " - " #variable " = " << variable;

void ControlOptions::Print(SeverityLogger& log) {
  FUB_GT_CONTROL_PRINT(efficiency_turbine);
  FUB_GT_CONTROL_PRINT(length_tube);
  FUB_GT_CONTROL_PRINT(surface_area_tube_inlet);
  FUB_GT_CONTROL_PRINT(surface_area_tube_outlet);
  FUB_GT_CONTROL_PRINT(surface_area_compressor_to_compressor_plenum);
  FUB_GT_CONTROL_PRINT(volume_compressor_plenum);
  FUB_GT_CONTROL_PRINT(surface_area_turbine_plenum_to_turbine);
  FUB_GT_CONTROL_PRINT(volume_turbine_plenum);
  FUB_GT_CONTROL_PRINT(efficiency_compressor);
  FUB_GT_CONTROL_PRINT(rpmmin);
  FUB_GT_CONTROL_PRINT(rpmmax);
  FUB_GT_CONTROL_PRINT(pratiomin);
  FUB_GT_CONTROL_PRINT(pratiomax);
  FUB_GT_CONTROL_PRINT(c_0);
  FUB_GT_CONTROL_PRINT(inertial_moment);
  FUB_GT_CONTROL_PRINT(mu);
  FUB_GT_CONTROL_PRINT(target_pressure_compressor);
  FUB_GT_CONTROL_PRINT(pressure_relDiffTolerance_SEC);
  FUB_GT_CONTROL_PRINT(Q);
  BOOST_LOG(log) << " - initial_turbine_temperature = "
                 << initial_turbine_state.temperature;
  BOOST_LOG(log) << " - initial_turbine_pressure = "
                 << initial_turbine_state.pressure;
}
#undef FUB_GT_CONTROL_PRINT

Control::Control(const PerfectGasConstants& eq, const ControlOptions& options)
    : equation_(eq), options_{options}, state_{
                                            std::make_shared<ControlState>()} {
  state_->current_rpm = options_.rpmmin;
  state_->fuel_consumption = 0.0;
  state_->power_out = 0.0;
  state_->turbine = options_.initial_turbine_state;
  UpdateCompressorFromRPM(*state_, equation_, options_);
}

void Control::UpdatePlena(double mdot_turbine,
                          const PlenumState& turbine_boundary_state,
                          double flux_rho, double flux_spec,
                          double flux_rho_last, Duration dt) {
  UpdateTurbinePlenum(*state_, equation_, options_, mdot_turbine,
                      turbine_boundary_state, flux_rho_last, dt);
  UpdateCompressorPlenum(*state_, equation_, options_, flux_rho, flux_spec, dt);
  ChangeRPM(*state_, options_, dt);
}

} // namespace fub::perfect_gas_mix::gt
