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

#include "fub/ext/ProgramOptions.hpp"

#include <utility>

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

namespace fub::amrex {

void PressureValveOptions::Print(SeverityLogger& log) {
  BOOST_LOG(log) << fmt::format("Pressure Valve '{}' Options:", channel);
  FUB_PRINT_OPTION_VAR(log, equivalence_ratio);
  FUB_PRINT_OPTION_VAR(log, pressure_value_which_closes_boundary);
  FUB_PRINT_OPTION_VAR(log, change_to_fuel_at_interval);
  FUB_PRINT_OPTION_VAR(log, change_to_fuel_time_offset);
  BOOST_LOG(log) << fmt::format(" --- Massflow part:");
  massflow_boundary.Print(log);
}

PressureValveOptions::PressureValveOptions(const ProgramOptions& options) {
  FUB_GET_OPTION_VAR(options, channel);
  FUB_GET_OPTION_VAR(options, equivalence_ratio);
  FUB_GET_OPTION_VAR(options, pressure_value_which_closes_boundary);
  FUB_GET_OPTION_VAR(options, change_to_fuel_at_interval);
  FUB_GET_OPTION_VAR(options, change_to_fuel_time_offset);
  massflow_boundary = GetOptions(options, "massflow_boundary");
}

PressureValveBoundary::PressureValveBoundary(const IdealGasMix<1>& equation,
                                             PressureValveOptions options)
    : options_{std::move(options)}, equation_{equation} {}

const PressureValveOptions& PressureValveBoundary::GetOptions() const noexcept {
  return options_;
}

const PressureValve& PressureValveBoundary::GetValve() const noexcept {
  return valve_;
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

[[maybe_unused]] std::vector<double>
GatherMoles_(const GriddingAlgorithm& grid, double x, IdealGasMix<1>& eq) {
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
      MPI_Allreduce(local_moles.data(), global_moles.data(), nspecies,
                    MPI_DOUBLE, MPI_SUM,
                    ::amrex::ParallelContext::CommunicatorAll());
      break;
    }
  }
  return global_moles;
}

double ChangeState_(PressureValveState& state, const ::amrex::Geometry&,
                    const GriddingAlgorithm& grid, int level,
                    Duration& last_fuel_change,
                    const PressureValveOptions& options, IdealGasMix<1>& eq) {
  const double mean_pressure = GetMeanPressure_(grid, eq);

  const Duration current_time = grid.GetPatchHierarchy().GetTimePoint(0);
  const Duration next_fuel_time =
      (last_fuel_change.count() < 0.0)
          ? options.change_to_fuel_time_offset
          : Duration(last_fuel_change + options.change_to_fuel_at_interval);
  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", options.channel);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", current_time.count());
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Level", level);
  switch (state) {
  case PressureValveState::closed:
    state = PressureValveState::open_air;
    BOOST_LOG(log) << "Pressure valve opened for air!";
    break;
  case PressureValveState::open_fuel:
    if (mean_pressure > options.pressure_value_which_closes_boundary) {
      state = PressureValveState::open_air;
      BOOST_LOG(log)
          << "Pressure valve changed to air due to high mean pressure (which is"
          << mean_pressure << " [Pa])!";
    }
    break;
  case PressureValveState::open_air:
    if (next_fuel_time < current_time) {
      state = PressureValveState::open_fuel;
      last_fuel_change = next_fuel_time;
      BOOST_LOG(log) << "Pressure valve changed to fuel!";
    }
  default:
    break;
  }
  return mean_pressure;
}

} // namespace

void PressureValveBoundary::FillBoundary(::amrex::MultiFab& mf,
                                         const GriddingAlgorithm& grid,
                                         int level, Direction dir) {
  if (dir == Direction::X) {
    FillBoundary(mf, grid, level);
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
  ChangeState_(valve_.state, geom, grid, level, valve_.last_fuel, options_,
               equation_);

  ReflectiveBoundary closed(execution::seq, equation_, Direction::X, 0);
  MassflowBoundary inflow_boundary(equation_, options_.massflow_boundary);
  Complete<IdealGasMix<1>> state{equation_};
  switch (valve_.state) {
  case PressureValveState::closed:
    closed.FillBoundary(mf, geom);
    break;
  case PressureValveState::open_air:
  case PressureValveState::open_fuel:
    inflow_boundary.FillBoundary(mf, grid, level);
  }
  if (valve_.state == PressureValveState::open_air) {
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
          ForEachIndex(
              box_to_fill, [this, &state, &states](std::ptrdiff_t i, auto...) {
                std::array<std::ptrdiff_t, 1> dest{i};
                Load(state, states, dest);
                FUB_ASSERT(state.density > 0.0);
                FUB_ASSERT(state.pressure > 0.0);
                FUB_ASSERT(state.temperature > 0.0);
                equation_.SetReactorStateFromComplete(state);
                FlameMasterReactor& reactor = equation_.GetReactor();
                reactor.SetDensity(1.0);
                reactor.SetMoleFractions("O2:21,N2:79");
                reactor.SetTemperature(300.0);
                reactor.SetPressure(state.pressure);
                Array<double, 1, 1> velocity = equation_.Velocity(state);
                equation_.CompleteFromReactor(state, velocity);
                Store(states, state, dest);
              });
        }
      }
    });
  } else if (valve_.state == PressureValveState::open_fuel) {
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
            span<const double> old_moles = reactor.GetMoleFractions();
            std::vector<double> new_moles(old_moles.begin(), old_moles.end());
            new_moles[Burke2012::sH2] =
                2.0 * options_.equivalence_ratio * new_moles[Burke2012::sO2];
            reactor.SetDensity(1.0);
            reactor.SetMoleFractions(new_moles);
            reactor.SetTemperature(state.temperature);
            reactor.SetPressure(state.pressure);
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
