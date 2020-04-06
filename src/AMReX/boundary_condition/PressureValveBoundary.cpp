// Copyright (c) 2019 Maikel Nadolski
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

#include "fub/AMReX/boundary_condition/PressureValveBoundary.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp"

#include <utility>

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

namespace fub::amrex {

void PressureValveOptions::Print(SeverityLogger& log) {
  BOOST_LOG(log) << fmt::format("Pressure Valve '{}' Options:", prefix);
  BOOST_LOG(log) << fmt::format(" - equivalence_ratio = {} [-]",
                                equivalence_ratio);
  BOOST_LOG(log) << fmt::format(" - outer_pressure = {} [Pa]", outer_pressure);
  BOOST_LOG(log) << fmt::format(" - outer_temperature = {} [K]",
                                outer_temperature);
  BOOST_LOG(log) << fmt::format(
      " - pressure_value_which_opens_boundary = {} [Pa]",
      pressure_value_which_opens_boundary);
  BOOST_LOG(log) << fmt::format(
      " - pressure_value_which_closes_boundary = {} [Pa]",
      pressure_value_which_closes_boundary);
  BOOST_LOG(log) << fmt::format(" - oxygen_measurement_position = {} [m]",
                                oxygen_measurement_position);
  BOOST_LOG(log) << fmt::format(" - oxygen_measurement_criterium = {} [mole]",
                                oxygen_measurement_criterium);
  BOOST_LOG(log) << fmt::format(" - fuel_measurement_position = {} [m]",
                                fuel_measurement_position);
  BOOST_LOG(log) << fmt::format(" - fuel_measurement_criterium = {} [-]",
                                fuel_measurement_criterium);
  BOOST_LOG(log) << fmt::format(" - valve_efficiency = {} [-]",
                                valve_efficiency);
  BOOST_LOG(log) << fmt::format(" - open_at_interval = {} [s]",
                                open_at_interval.count());
  BOOST_LOG(log) << fmt::format(" - offset = {} [s]", offset.count());
  BOOST_LOG(log) << fmt::format(" --- Massflow part:");
  massflow_boundary.Print(log);
}

PressureValveOptions::PressureValveOptions(const ProgramOptions& opts) {
  prefix = GetOptionOr(opts, "prefix", prefix);
  outer_pressure = GetOptionOr(opts, "outer_pressure", outer_pressure);
  outer_temperature = GetOptionOr(opts, "outer_temperature", outer_temperature);
  pressure_value_which_opens_boundary =
      GetOptionOr(opts, "pressure_value_which_opens_boundary",
                  pressure_value_which_opens_boundary);
  pressure_value_which_closes_boundary =
      GetOptionOr(opts, "pressure_value_which_closes_boundary",
                  pressure_value_which_closes_boundary);
  oxygen_measurement_position = GetOptionOr(opts, "oxygen_measurement_position",
                                            oxygen_measurement_position);
  oxygen_measurement_criterium = GetOptionOr(
      opts, "oxygen_measurement_criterium", oxygen_measurement_criterium);
  equivalence_ratio = GetOptionOr(opts, "equivalence_ratio", equivalence_ratio);
  open_at_interval =
      Duration(GetOptionOr(opts, "open_at_interval", open_at_interval.count()));
  offset = Duration(GetOptionOr(opts, "offset", offset.count()));
  valve_efficiency = GetOptionOr(opts, "efficiency", valve_efficiency);
  fuel_measurement_criterium = GetOptionOr(opts, "fuel_measurement_criterium",
                                           fuel_measurement_criterium);
  fuel_measurement_position =
      GetOptionOr(opts, "fuel_measurement_position", fuel_measurement_position);
  massflow_boundary = GetOptions(opts, "massflow_boundary");
}

PressureValveBoundary::PressureValveBoundary(const IdealGasMix<1>& equation,
                                             PressureValveOptions options)
    : options_{std::move(options)}, equation_{equation},
      shared_valve_{std::make_shared<PressureValve>(PressureValve{})} {}

PressureValveBoundary::PressureValveBoundary(
    const IdealGasMix<1>& eq,
    const std::map<std::string, pybind11::object>& map)
    : PressureValveBoundary(eq, PressureValveOptions(map)) {}

const PressureValveOptions& PressureValveBoundary::GetOptions() const noexcept {
  return options_;
}

namespace {

double GetMeanPressure_(const GriddingAlgorithm& grid, IdealGasMix<1>& eq) {
  const PatchHierarchy& hier = grid.GetPatchHierarchy();
  const ::amrex::MultiFab& data = hier.GetPatchLevel(0).data;
  const int n_cells = hier.GetGeometry(0).Domain().length(0);
  double local_pressure = 0.0;
  Complete<IdealGasMix<1>> state(eq);
  ForEachFab(data, [&local_pressure, &eq, n_cells, &data,
                    &state](const ::amrex::MFIter& mfi) {
    const ::amrex::FArrayBox& fab = data[mfi];
    auto view = MakeView<const Complete<IdealGasMix<1>>>(fab, eq);
    ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i) {
      Load(state, view, {i});
      local_pressure += state.pressure / n_cells;
    });
  });
  double pressure = 0.0;
  MPI_Allreduce(&local_pressure, &pressure, 1, MPI_DOUBLE, MPI_SUM,
                ::amrex::ParallelContext::CommunicatorAll());
  return pressure;
}

std::vector<double> GatherMoles_(const GriddingAlgorithm& grid, double x,
                                 IdealGasMix<1>& eq) {
  const PatchHierarchy& hier = grid.GetPatchHierarchy();
  const int nlevel = hier.GetNumberOfLevels();
  std::size_t size = static_cast<std::size_t>(eq.GetReactor().GetNSpecies());
  std::vector<double> global_moles(size);
  std::vector<double> local_moles(size);
  for (int level = nlevel - 1; level >= 0; --level) {
    const ::amrex::MultiFab& data = hier.GetPatchLevel(level).data;
    const ::amrex::Geometry& geom = hier.GetGeometry(level);
    int local_found = 0;
    Complete<IdealGasMix<1>> state(eq);
    ForEachFab(data, [&local_moles, &local_found, &eq, &data, &geom, x,
                      &state](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& fab = data[mfi];
      auto view = MakeView<const Complete<IdealGasMix<1>>>(fab, eq);
      ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i) {
        ::amrex::IntVect iv(AMREX_D_DECL(int(i), 0, 0));
        double x_lo[AMREX_SPACEDIM]{};
        geom.LoNode(iv, x_lo);
        double x_hi[AMREX_SPACEDIM]{};
        geom.HiNode(iv, x_hi);
        const double lo = x_lo[0];
        const double hi = x_hi[0];
        if (lo <= x && x < hi) {
          local_found = 1;
          Load(state, view, {i});
          eq.SetReactorStateFromComplete(state);
          span<const double> moles = eq.GetReactor().GetMoleFractions();
          std::copy(moles.begin(), moles.end(), local_moles.begin());
        }
      });
    });
    int found = 0;
    MPI_Allreduce(&local_found, &found, 1, MPI_INT, MPI_SUM,
                  ::amrex::ParallelContext::CommunicatorAll());
    if (found) {
      const int nspecies = static_cast<int>(size);
      MPI_Allreduce(local_moles.data(), global_moles.data(), nspecies, MPI_DOUBLE,
                    MPI_SUM, ::amrex::ParallelContext::CommunicatorAll());
      break;
    }
  }
  return global_moles;
}

double ChangeState_(PressureValveState& state, const ::amrex::Geometry& geom,
                    const GriddingAlgorithm& grid, int level,
                    Duration& last_closed, Duration& last_fuel_change,
                    const PressureValveOptions& options, IdealGasMix<1>& eq) {
  const double mean_pressure = GetMeanPressure_(grid, eq);
  const double dx_2 = geom.CellSize(0);
  const double xlo = geom.ProbDomain().lo(0) + dx_2;
  const double xhi = geom.ProbDomain().hi(0) - dx_2;
  const Duration current_time = grid.GetPatchHierarchy().GetTimePoint(0);
  const Duration next_open_time = (last_closed.count() < 0.0)
                                      ? Duration(0.0)
                                      : last_closed + Duration(1e-3);
  const Duration next_fuel_time =
      (last_fuel_change.count() < 0.0)
          ? options.offset
          : last_fuel_change + options.open_at_interval;
  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", options.prefix);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", current_time.count());
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Level", level);
  switch (state) {
  case PressureValveState::closed:
    if (next_open_time < current_time &&
        mean_pressure < options.pressure_value_which_opens_boundary) {
      state = PressureValveState::open_air;
      BOOST_LOG(log) << "pressure valve opened for air!";
    }
    break;
  case PressureValveState::open_air:
  case PressureValveState::open_fuel:
    if (mean_pressure > options.pressure_value_which_closes_boundary) {
      last_closed = current_time;
      state = PressureValveState::closed;
      BOOST_LOG(log) << "pressure valve closed due to mean pressure (which is "
                     << mean_pressure << " [Pa])!";
    }
    break;
  }
  if (next_fuel_time < current_time && state == PressureValveState::open_air) {
    const double x_air =
        std::clamp(options.oxygen_measurement_position, xlo, xhi);
    std::vector<double> moles = GatherMoles_(grid, x_air, eq);
    const double sum = std::accumulate(moles.begin(), moles.end(), 0.0);
    FUB_ASSERT(sum >= 0.0);
    if (sum > 0.0) {
      std::transform(moles.begin(), moles.end(), moles.begin(),
                     [sum](double m) { return m / sum; });
    }
    if (options.oxygen_measurement_criterium < moles[Burke2012::sO2]) {
      state = PressureValveState::open_fuel;
      last_fuel_change = next_fuel_time;
      BOOST_LOG(log) << "pressure valve changed to fuel!";
    }
  } else if (state == PressureValveState::open_fuel) {
    const double x_fuel =
        std::clamp(options.fuel_measurement_position, xlo, xhi);
    std::vector<double> moles = GatherMoles_(grid, x_fuel, eq);
    const double equivalence_ratio =
        moles[Burke2012::sO2]
            ? 0.5 * moles[Burke2012::sH2] / moles[Burke2012::sO2]
            : 0.0;
    if (equivalence_ratio > options.fuel_measurement_criterium) {
      last_closed = current_time;
      state = PressureValveState::closed;
      BOOST_LOG(log)
          << "pressure valve closed due to measured fuel: equivalence ratio at "
          << x_fuel << " [m] is " << equivalence_ratio << " [-]!";
    }
  }

  return mean_pressure;
}

} // namespace

void PressureValveBoundary::FillBoundary(::amrex::MultiFab& mf,
                                         const GriddingAlgorithm& grid,
                                         int level, Direction dir) {
  if (dir == Direction::X) {
    FillBoundary(mf, grid, level, dir);
  }
}

void PressureValveBoundary::FillBoundary(::amrex::MultiFab& mf,
                                         const GriddingAlgorithm& grid,
                                         int level) {
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  const int ngrow = mf.nGrow(0);
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }

  // Change State Machine if neccessary
  ChangeState_(shared_valve_->state, geom, grid, level,
               shared_valve_->last_closed, shared_valve_->last_fuel, options_,
               equation_);

  ReflectiveBoundary closed(execution::openmp, equation_, Direction::X, 0);
  MassflowBoundary inflow_boundary(equation_, options_.massflow_boundary);
  Complete<IdealGasMix<1>> state{equation_};
  switch (shared_valve_->state) {
  case PressureValveState::closed:
    closed.FillBoundary(mf, geom);
    break;
  case PressureValveState::open_air:
  case PressureValveState::open_fuel:
    FlameMasterReactor& reactor = equation_.GetReactor();
    reactor.SetDensity(1.0);
    reactor.SetMoleFractions("O2:21,N2:79");
    reactor.SetTemperature(300.0);
    reactor.SetPressure(101325.0);
    equation_.CompleteFromReactor(state);
    inflow_boundary.ComputeBoundaryState(state, state);
    inflow_boundary.FillBoundary(mf, geom, state);
  }
  if (shared_valve_->state == PressureValveState::open_fuel) {
    ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& fab = mf[mfi];
      for (const ::amrex::Box& boundary : boundaries) {
        ::amrex::Box shifted = ::amrex::shift(boundary, 0, ngrow);
        if (!geom.Domain().intersects(shifted)) {
          continue;
        }
        ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
        if (!box_to_fill.isEmpty()) {
          auto states = MakeView<Complete<IdealGasMix<1>>>(fab, equation_,
                                                           mfi.growntilebox());
          ForEachIndex(box_to_fill, [this, &state, &states](std::ptrdiff_t i,
                                                            auto...) {
            std::array<std::ptrdiff_t, 1> dest{i};
            Load(state, states, dest);
            FUB_ASSERT(state.density > 0.0);
            FUB_ASSERT(state.pressure > 0.0);
            FUB_ASSERT(state.temperature > 0.0);
            equation_.SetReactorStateFromComplete(state);
            FlameMasterReactor& reactor = equation_.GetReactor();
            span<const double> mass_over_mole =
                reactor.GetMolarMasses(); // kg / kmol
            span<const double> old_moles = reactor.GetMoleFractions();
            std::vector<double> new_moles(old_moles.begin(), old_moles.end());
            new_moles[Burke2012::sH2] =
                2.0 * options_.equivalence_ratio * new_moles[Burke2012::sO2];
            const double difference_h2_moles =
                new_moles[Burke2012::sH2] - old_moles[Burke2012::sH2];
            const double new_density =
                state.density +
                difference_h2_moles * mass_over_mole[Burke2012::sH2];
            reactor.SetDensity(new_density);
            reactor.SetMoleFractions(new_moles);
            reactor.SetTemperature(state.temperature);
            Array<double, 1, 1> velocity = equation_.Velocity(state);
            equation_.CompleteFromReactor(state, velocity);
            Store(states, state, dest);
          });
        }
      }
    });
  }
}
} // namespace fub::amrex
