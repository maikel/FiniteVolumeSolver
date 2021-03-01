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

#ifndef FUB_EQUATIONS_PERFECT_GAS_MIX_PLENA_CONTROL_FEEDBACK_HPP
#define FUB_EQUATIONS_PERFECT_GAS_MIX_PLENA_CONTROL_FEEDBACK_HPP

#include "fub/AMReX/multi_block/MultiBlockIntegratorContext2.hpp"
#include "fub/equations/perfect_gas_mix/PlenaControl.hpp"

#include "fub/AMReX/AverageState.hpp"

namespace fub::perfect_gas_mix::gt {

template <int PlenumRank> class ControlFeedback {
public:
  ControlFeedback(const PerfectGasMix<PlenumRank>& plenum_equation,
                  const PerfectGasMix<1>& tube_equation, const Control& c,
                  const ::amrex::Box boundary_box)
      : plenum_equation_(plenum_equation), tube_index_(tube_equation),
        control_(c), boundary_box_(boundary_box) {}

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
    const ::amrex::Box cells =
        !boundary_box_.isEmpty() ? boundary_box_ : geom.Domain();
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
  ::amrex::Box boundary_box_;
};

struct ControlFeedbackOptions {
  ControlFeedbackOptions() = default;

  ControlFeedbackOptions(const ProgramOptions& options) {
    FUB_GET_OPTION_VAR(options, mdot_correlation);
  }

  void Print(SeverityLogger& log) const {
    FUB_PRINT_OPTION_VAR(log, mdot_correlation);
  }

  double mdot_correlation{0.24};
};

template <> class ControlFeedback<1> {
public:
  ControlFeedback(const PerfectGasMix<1>& equation, const Control& control,
                  const ControlFeedbackOptions& options)
      : equation_{equation},
        tube_index_{equation}, control_{control}, options_{options} {}

  double UpdateTurbinePlenumState(PlenumState& turbine_old, PlenumState& turbine_new, double f_rho,
                                  double f_rhoe, Duration dt,
                                  const ControlOptions& options) const {
    const double pT = turbine_old.pressure;
    const double TT = turbine_old.temperature;
    const double mdotT = options_.mdot_correlation * pT / std::sqrt(TT);

    // Update plenum state
    const double rho = pT / TT * equation_.ooRspec;
    const double rhoe = pT * equation_.gamma_minus_one_inv;
    const double h = (rhoe + pT) / rho;

    const double rho_n =
        rho + (f_rho * options.surface_area_tube_outlet - mdotT) * dt.count() /
                  options.volume_turbine_plenum;
    const double rhoe_n =
        rhoe + (f_rhoe * options.surface_area_tube_outlet - mdotT * h) *
                   dt.count() / options.volume_turbine_plenum;

    turbine_new.pressure = rhoe_n / equation_.gamma_minus_one_inv;
    turbine_new.temperature = turbine_new.pressure / rho_n * equation_.ooRspec;

    return mdotT;
  }

  void operator()(amrex::IntegratorContext& tube_context, Duration dt) {
    // Get inflow and outflow fluxes from the combustion tube
    double flux_rho{};
    double flux_species{};
    double flux_rho_last{};
    double flux_rhoe_last{};
    constexpr int coarsest_level = 0;
    const ::amrex::Geometry& geom = tube_context.GetGeometry(coarsest_level);
    const ::amrex::Box face_domain = ::amrex::convert(
        geom.Domain(), ::amrex::IntVect::TheDimensionVector(0));
    const ::amrex::IntVect first_face = face_domain.smallEnd();
    const ::amrex::IntVect last_face = face_domain.bigEnd();
    const ::amrex::MultiFab& fluxes_x =
        tube_context.GetFluxes(coarsest_level, Direction::X);
    double local_f_rho = 0.0;
    double local_f_rho_last = 0.0;
    double local_f_rhoe_last = 0.0;
    double local_f_spec = 0.0;
    amrex::ForEachFab(fluxes_x, [&](const ::amrex::MFIter& mfi) {
      if (mfi.tilebox().contains(first_face)) {
        local_f_rho = fluxes_x[mfi](first_face, tube_index_.density);
        local_f_spec = fluxes_x[mfi](first_face, tube_index_.species[0]);
      }
      if (mfi.tilebox().contains(last_face)) {
        local_f_rho_last = fluxes_x[mfi](last_face, tube_index_.density);
        local_f_rhoe_last = fluxes_x[mfi](last_face, tube_index_.energy);
      }
    });
    ::MPI_Allreduce(&local_f_rho, &flux_rho, 1, MPI_DOUBLE, MPI_SUM,
                    ::amrex::ParallelDescriptor::Communicator());
    ::MPI_Allreduce(&local_f_rho_last, &flux_rho_last, 1, MPI_DOUBLE, MPI_SUM,
                    ::amrex::ParallelDescriptor::Communicator());
    ::MPI_Allreduce(&local_f_spec, &flux_species, 1, MPI_DOUBLE, MPI_SUM,
                    ::amrex::ParallelDescriptor::Communicator());
    ::MPI_Allreduce(&local_f_rhoe_last, &flux_rhoe_last, 1, MPI_DOUBLE, MPI_SUM,
                    ::amrex::ParallelDescriptor::Communicator());

    PlenumState turbine_boundary_state = control_.GetSharedState()->turbine;
    PlenumState new_turbine_boundary_state;
    const double mdot_turbine =
        UpdateTurbinePlenumState(turbine_boundary_state, new_turbine_boundary_state, flux_rho_last,
                                 flux_rhoe_last, dt, control_.GetOptions());

    control_.UpdatePlena(mdot_turbine, new_turbine_boundary_state, flux_rho,
                         flux_species, flux_rho_last, dt);
  }

private:
  PerfectGasMix<1> equation_;
  IndexMapping<PerfectGasMix<1>> tube_index_;
  Control control_;
  ControlFeedbackOptions options_;
};

} // namespace fub::perfect_gas_mix::gt

#endif