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

#include "fub/AMReX/AverageState.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/multi_block/MultiBlockIntegratorContext2.hpp"

#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include "fub/output/Hdf5Handle.hpp"

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

template <int PlenumRank> class ControlFeedback {
public:
  ControlFeedback(const PerfectGasMix<PlenumRank>& plenum_equation,
                  const PerfectGasMix<1>& tube_equation, const Control& c)
      : plenum_equation_(plenum_equation), tube_index_(tube_equation),
        control_(c) {}

  // Update the control object using the mass flows in the multi block
  // integrator context.
  void operator()(amrex::MultiBlockIntegratorContext2& context, Duration dt) {
    /* grab the mean global mass flux into the turbine */
    amrex::cutcell::IntegratorContext& plenum = context.Plena()[0];
    int coarsest_level = 0;
    // const Duration time_point = plenum.GetTimePoint(coarsest_level);
    const ::amrex::MultiFab& scratch = plenum.GetScratch(coarsest_level);
    const ::amrex::MultiFab& fluxes_x =
        plenum.GetFluxes(coarsest_level, Direction::X);
    const ::amrex::Geometry& geom = plenum.GetGeometry(coarsest_level);
    const ::amrex::Box cells = geom.Domain();
    const int boundary_n = cells.bigEnd(0);
    const ::amrex::IntVect smallEnd{boundary_n, cells.smallEnd(1)};
    const ::amrex::IntVect bigEnd = cells.bigEnd();
    const ::amrex::Box right_cell_boundary{smallEnd, bigEnd};
    const ::amrex::Box faces_x = ::amrex::convert(
        right_cell_boundary, ::amrex::IntVect{AMREX_D_DECL(1, 0, 0)});
    const ::amrex::IntVect fsmallEnd = faces_x.smallEnd();
    const ::amrex::IntVect fbigEnd = {faces_x.smallEnd(0), faces_x.bigEnd(1)};
    const ::amrex::Box right_face_boundary{fsmallEnd, fbigEnd,
                                           faces_x.ixType()};
    Conservative<PerfectGasMix<2>> fluxes_turbine(plenum_equation_);
    Conservative<PerfectGasMix<2>> cons_turbine(plenum_equation_);
    amrex::AverageState(fluxes_turbine, fluxes_x, geom, right_face_boundary);
    amrex::AverageState(cons_turbine, scratch, geom, right_cell_boundary);

    PlenumState turbine_boundary_state;
    turbine_boundary_state.pressure =
        euler::Pressure(plenum_equation_, cons_turbine);
    turbine_boundary_state.temperature = turbine_boundary_state.pressure /
                                         cons_turbine.density *
                                         plenum_equation_.ooRspec;

    const double mdot_turbine =
        fluxes_turbine.density *
        control_.GetOptions().surface_area_turbine_plenum_to_turbine;

    // For each tube we get the mass flux into the first cell
    span<amrex::IntegratorContext> tubes = context.Tubes();
    // (rhou)_{I + 1/2}
    std::vector<double> flux_rho_tube(tubes.size());
    std::vector<double> flux_species_tube(tubes.size());
    std::vector<double> flux_rho_last_tube(tubes.size());

    for (auto&& [tube_context, flux_rho, flux_spec, flux_rho_last] :
         ranges::view::zip(tubes, flux_rho_tube, flux_species_tube,
                           flux_rho_last_tube)) {
      const ::amrex::Geometry& tube_geom =
          tube_context.GetGeometry(coarsest_level);
      const ::amrex::Box face_domain = ::amrex::convert(
          tube_geom.Domain(), ::amrex::IntVect::TheDimensionVector(0));
      const ::amrex::IntVect first_face = face_domain.smallEnd();
      const ::amrex::IntVect last_face = face_domain.bigEnd();
      const ::amrex::MultiFab& fluxes_x =
          tube_context.GetFluxes(coarsest_level, Direction::X);
      double local_f_rho = 0.0;
      double local_f_rho_last = 0.0;
      double local_f_spec = 0.0;
      amrex::ForEachFab(fluxes_x, [&](const ::amrex::MFIter& mfi) {
        if (mfi.tilebox().contains(first_face)) {
          local_f_rho = fluxes_x[mfi](first_face, tube_index_.density);
          local_f_spec = fluxes_x[mfi](first_face, tube_index_.species[0]);
        }
        if (mfi.tilebox().contains(last_face)) {
          local_f_rho_last = fluxes_x[mfi](last_face, tube_index_.density);
        }
      });
      ::MPI_Allreduce(&local_f_rho, &flux_rho, 1, MPI_DOUBLE, MPI_SUM,
                      ::amrex::ParallelDescriptor::Communicator());
      ::MPI_Allreduce(&local_f_rho_last, &flux_rho_last, 1, MPI_DOUBLE, MPI_SUM,
                      ::amrex::ParallelDescriptor::Communicator());
      ::MPI_Allreduce(&local_f_spec, &flux_spec, 1, MPI_DOUBLE, MPI_SUM,
                      ::amrex::ParallelDescriptor::Communicator());
    }
    const double oosize = 1.0 / static_cast<double>(tubes.size());
    const double rhou_tubes =
        oosize *
        std::accumulate(flux_rho_tube.begin(), flux_rho_tube.end(), 0.0);
    const double specu_tubes =
        oosize * std::accumulate(flux_species_tube.begin(),
                                 flux_species_tube.end(), 0.0);
    const double rhou_last_tubes =
        oosize * std::accumulate(flux_rho_last_tube.begin(),
                                 flux_rho_last_tube.end(), 0.0);

    control_.UpdatePlena(mdot_turbine, turbine_boundary_state, rhou_tubes,
                         specu_tubes, rhou_last_tubes, dt);
  }

private:
  PerfectGasMix<PlenumRank> plenum_equation_;
  IndexMapping<PerfectGasMix<1>> tube_index_;
  Control control_;
};

class ControlOutput
    : public OutputAtFrequencyOrInterval<amrex::MultiBlockGriddingAlgorithm2> {
public:
  ControlOutput(const ProgramOptions& options,
                std::shared_ptr<const ControlState> control_state);

  void operator()(const amrex::MultiBlockGriddingAlgorithm2& grid) override;

private:
  std::string file_path_;
  std::shared_ptr<const ControlState> control_state_;
  std::map<std::string, function_ref<double(const ControlState&)>> fields_;
  std::vector<double> data_buffer_;

  void CreateHdf5Database();
  void WriteHdf5Database(span<const double> data, Duration time,
                         std::ptrdiff_t cycle);
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