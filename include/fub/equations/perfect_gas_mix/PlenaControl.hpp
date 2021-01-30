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

#include "fub/equations/PerfectGasMix.hpp"
#include "fub/equations/perfect_gas_mix/ArrheniusKinetics.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/AverageState.hpp"
#include "fub/AMReX/multi_block/MultiBlockIntegratorContext2.hpp"

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
                                   double tign, double Y, double T0, double eps) {
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
  double pressure;
  double temperature;
};

template <int PlenumRank> class Control {
public:
  static constexpr int TubeRank = 1;
  // Constructors

  Control(const PerfectGasMix<PlenumRank>& eq) : equation_(eq) {
    compressor_plenum_ = std::make_shared<PlenumState>();
    compressor_plenum_->pressure = 2.0;
    compressor_plenum_->temperature = 1.0;
    turbine_plenum_old_.pressure = 1.0;
    turbine_plenum_old_.temperature = 1.0;
  }

  void UpdatePlena(double mdot_turbine,
                   const PlenumState& turbine_boundary_state, double flux_rho,
                   Duration dt) {
    UpdateTurbinePlenum(mdot_turbine, turbine_boundary_state);
    UpdateCompressorPlenum(flux_rho, dt);
    ChangeRPM(dt);
  }

  std::shared_ptr<const PlenumState> GetSharedCompressorState() const noexcept {
    return compressor_plenum_;
  }

private:
  /// Notify the control with a new turbine plenum value and update its current
  /// power.
  ///
  /// The formula is given in GT_Notes equation (35).
  ///
  /// \return Returns the new power value for the turbine
  void UpdateTurbinePlenum(double mdot,
                           const PlenumState& turbine_plenum_current) noexcept {
    PlenumState T_half;
    T_half.pressure =
        0.5 * (turbine_plenum_old_.pressure + turbine_plenum_current.pressure);
    T_half.temperature = 0.5 * (turbine_plenum_old_.temperature +
                                turbine_plenum_current.temperature);
    const double gmoog = equation_.gamma_minus_one_over_gamma;
    const double cp = equation_.heat_capacity_at_constant_pressure;
    power_turbine_ = mdot * cp * efficiency_turbine_ * T_half.temperature *
                     (1.0 - std::pow(T_half.pressure, gmoog));
    turbine_plenum_old_ = turbine_plenum_current;
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

  /// Notify the control with a new compressor plenum value and update its
  /// current power.
  ///
  /// The formula is given in GT_Notes equation (24).
  ///
  /// \return Returns the new power value
  void UpdateCompressorPlenum(double flux_rho, Duration dt) {
    double temperature_old = compressor_plenum_->temperature;
    double rho_old = compressor_plenum_->pressure /
                     compressor_plenum_->temperature * equation_.ooRspec;

    // update the compressor_plenum state with the new speed of the compressor
    // so we get new pressure and temperature
    UpdateCompressorFromRPM(current_rpm_);

    double rho_new = compressor_plenum_->pressure /
                     compressor_plenum_->temperature / equation_.Rspec;
    double mdot_compressor = flux_rho * surface_area_tube_inlet_ +
                             volume_compressor_plenum_ / dt.count() *
                                 (rho_new - rho_old); // equation (23) GT_Notes

    double T_v_half = 0.5 * (temperature_old + compressor_plenum_->temperature);
    double T_ref = 1.0;
    const double cp = equation_.heat_capacity_at_constant_pressure;

    double power_compressor_increase =
        mdot_compressor * (T_v_half - T_ref) / efficiency_compressor_ * cp;

    const double ootau = 1.0;
    double dt_over_tau = dt.count() * ootau;

    power_compressor_ = (1.0 - dt_over_tau) * power_compressor_ +
                        dt_over_tau * power_compressor_increase;
  }

  /// \brief Compute the compressor state by given speed of the compressor.
  void UpdateCompressorFromRPM(double rpm) {
    double new_pressure = CompressorPressureRatio(rpm);
    compressor_plenum_->pressure = new_pressure;
    const double gmoog = equation_.gamma_minus_one_over_gamma;
    compressor_plenum_->temperature =
        1.0 + std::pow(new_pressure, gmoog) / efficiency_compressor_;
  }

  /// \brief Compute the new rotation speed of the compressor.
  ///
  /// see GT_notes equation (42 - 44)
  void ChangeRPM(Duration dt) {
    double pressure = compressor_plenum_->pressure;
    const double pressure_ratio =
        (target_pressure_compressor_ - pressure) / target_pressure_compressor_;
    const double exp_term = std::exp(-mu_ * pressure_ratio);
    const double comp_rate = 0.015 * (1.0 - exp_term);
    constexpr double pi2 = M_PI * M_PI;
    const double Ieff = 2.0 * pi2 * inertial_moment_;

    double power_netto = power_turbine_ - power_compressor_;
    double d_Erot_dt = (pressure < target_pressure_compressor_)
                           ? comp_rate * power_netto
                           : 0.0;
    double Erot = Ieff * current_rpm_ * current_rpm_ + d_Erot_dt * dt.count();
    double rpm_new = std::sqrt(Erot / Ieff);

    FUB_ASSERT(rpmmin_ < rpm_new && rpmmax_ > rpm_new);
    current_rpm_ = rpm_new;

    // values to track:
    // double power_out = power_netto - d_Erot_dt; // Power output

    // efficiency of the machine
    // GT->Efficiency[0] = MAX_own(0.0, GT->PowerOut[0]) / (ud.Q *
    // GT->FuelConsumptionRate[0] + ud.eps_Machine);
  }

private:
  template <int Rank> friend class ControlFeedback;

  PerfectGasMix<PlenumRank> equation_;
  PlenumState turbine_plenum_old_{}; ///< the state of the turbine plenum
  std::shared_ptr<PlenumState> compressor_plenum_{};  ///< the state of the compressor plenum
  double efficiency_turbine_{0.9}; ///< efficiency from the turbine

  double length_tube_{1.0};                  ///< the length of the tubes
  double surface_area_tube_inlet_{1.0}; ///< initial surface area of the tube

  /// surface area from the compressor to the compressor plenum
  double surface_area_compressor_to_compressor_plenum_{
      8.0 * surface_area_tube_inlet_};
  /// volume of the compressor plenum
  double volume_compressor_plenum_{20.0 * surface_area_tube_inlet_ *
                                   length_tube_};
  /// surface area from the turbine plenum to the turbine
  double surface_area_turbine_plenum_to_turbine_{4.0 *
                                                 surface_area_tube_inlet_};
  /// volume of the turbine plenum
  double volume_turbine_plenum_{20.0 * surface_area_tube_inlet_ * length_tube_};

  double power_turbine_{};              ///< power from the turbine
  double power_compressor_{};           ///< power from the compressor
  double efficiency_compressor_{0.9}; ///< efficiency from the compressor
  double rpmmin_ = 10000.0 / 60.0;    ///< minimum speed of the compressor
  double rpmmax_ = 60000.0 / 60.0;    ///< maximum speed of the compressor
  double rpmmean_ = 0.5 * (rpmmax_ + rpmmin_); ///< mean speed of the compressor
  double pratiomin_ = 1.0;  ///< minimum pressure from the compressor
  double pratiomax_ = 50.0; ///< maximum pressure from the compressor
  /// mean pressure from the compressor
  double pratiomean_ = 0.5 * (pratiomax_ + pratiomin_);
  /// maximum pressure difference from the compressor
  double pratiovar_ = pratiomax_ - pratiomin_;

  /// constant to compute the pressure from the speed of the compressor (see
  /// function: CompressorPressureRatio)
  double c_0_ = 22.0 / 17.0;

  double current_rpm_{rpmmin_}; ///< current speed of the compressor

  double inertial_moment_{0.1}; ///< the initial moment of inertia of the rotor
  double target_pressure_compressor_{
      6.5};         ///< the target pressure from the compressor
  double mu_{10.0}; ///< some constant for GT_control
};

template <int PlenumRank> class ControlFeedback {
public:
  ControlFeedback(
      const PerfectGasMix<PlenumRank>& plenum_equation,
      const PerfectGasMix<1>& tube_equation,
      const Control<PlenumRank>& c)
      : plenum_equation_(plenum_equation), tube_equation_(tube_equation),
        control_(c) {}

  // Update the control object using the mass flows in the multi block
  // integrator context.
  void operator()(amrex::MultiBlockIntegratorContext2& context, Duration dt) {
    /* grab the mean global mass flux into the turbine */
    amrex::cutcell::IntegratorContext& plenum = context.Plena()[0];
    int coarsest_level = 0;
    const Duration time_point = plenum.GetTimePoint(coarsest_level);
    const ::amrex::MultiFab& scratch = plenum.GetScratch(coarsest_level);
    const ::amrex::MultiFab& fluxes_x =
        plenum.GetFluxes(coarsest_level, Direction::X);
    const ::amrex::Geometry& geom = plenum.GetGeometry(coarsest_level);
    const ::amrex::Box cells = geom.Domain();
    const int boundary_n = cells.bigEnd(0);
    const ::amrex::IntVect smallEnd{boundary_n, cells.smallEnd(1)};
    const ::amrex::IntVect bigEnd = cells.bigEnd();
    const ::amrex::Box right_cell_boundary{smallEnd, bigEnd};
    const ::amrex::Box faces_x = ::amrex::convert(right_cell_boundary, {1, 0});
    const ::amrex::IntVect fsmallEnd = faces_x.smallEnd();
    const ::amrex::IntVect fbigEnd = {faces_x.smallEnd(0), faces_x.bigEnd(1)};
    const ::amrex::Box right_face_boundary{fsmallEnd, fbigEnd, faces_x.ixType()};
    Conservative<PerfectGasMix<2>> fluxes_turbine(plenum_equation_);
    Conservative<PerfectGasMix<2>> cons_turbine(plenum_equation_);
    amrex::AverageState(fluxes_turbine, fluxes_x, geom, right_face_boundary);
    amrex::AverageState(cons_turbine, scratch, geom, right_cell_boundary);

    PlenumState turbine_boundary_state;
    turbine_boundary_state.pressure =
        euler::Pressure(plenum_equation_, cons_turbine);
    turbine_boundary_state.temperature =
        turbine_boundary_state.pressure / cons_turbine.density * plenum_equation_.ooRspec;

    const double mdot_turbine =
        fluxes_turbine.density *
        control_.surface_area_turbine_plenum_to_turbine_;

    // For each tube we get the mass flux into the first cell
    span<amrex::IntegratorContext> tubes = context.Tubes();
    // (rhou)_{I + 1/2}
    std::vector<double> flux_rho_tube(tubes.size());

    for (auto&& [tube_context, flux_tube] :
         ranges::view::zip(tubes, flux_rho_tube)) {
      const ::amrex::MultiFab& fluxes_x =
          tube_context.GetFluxes(coarsest_level, Direction::X);
      double local_f_rho = 0.0;
      amrex::ForEachFab(fluxes_x, [&](const ::amrex::MFIter& mfi) {
        ::amrex::IntVect iv{0, 0};
        if (mfi.tilebox().contains(iv)) {
          local_f_rho = fluxes_x[mfi](iv, 0);
        }
      });
      ::MPI_Allreduce(&local_f_rho, &flux_tube, 1, MPI_DOUBLE, MPI_SUM,
                      ::amrex::ParallelDescriptor::Communicator());
    }
    const double oosize = 1.0 / static_cast<double>(tubes.size());
    double rhou_tubes = oosize *
        std::accumulate(flux_rho_tube.begin(), flux_rho_tube.end(), 0.0);

    control_.UpdatePlena(mdot_turbine, turbine_boundary_state, rhou_tubes, dt);
  }

private:
  PerfectGasMix<PlenumRank> plenum_equation_;
  PerfectGasMix<1> tube_equation_;
  Control<PlenumRank> control_;
};

} // namespace gt

} // namespace fub::perfect_gas_mix

#endif