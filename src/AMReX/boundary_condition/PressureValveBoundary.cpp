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

PressureValveBoundary::PressureValveBoundary(const IdealGasMix<1>& equation,
                                             PressureValveOptions options)
    : options_{std::move(options)}, equation_{equation},
      shared_valve_{std::make_shared<PressureValve>(PressureValve{})} {}

PressureValveOptions::PressureValveOptions(const ProgramOptions& opts)
{
    prefix = GetOptionOr(opts, "prefix", prefix);
    outer_pressure = GetOptionOr(opts, "outer_pressure", outer_pressure);
    outer_temperature =
        GetOptionOr(opts, "outer_temperature", outer_temperature);
    pressure_value_which_opens_boundary =
        GetOptionOr(opts, "pressure_value_which_opens_boundary",
                    pressure_value_which_opens_boundary);
    pressure_value_which_closes_boundary =
        GetOptionOr(opts, "pressure_value_which_closes_boundary",
                    pressure_value_which_closes_boundary);
    oxygen_measurement_position = GetOptionOr(
        opts, "oxygen_measurement_position", oxygen_measurement_position);
    oxygen_measurement_criterium = GetOptionOr(
        opts, "oxygen_measurement_criterium", oxygen_measurement_criterium);
    equivalence_ratio =
        GetOptionOr(opts, "equivalence_ratio", equivalence_ratio);
    open_at_interval = Duration(
        GetOptionOr(opts, "open_at_interval", open_at_interval.count()));
    offset = Duration(GetOptionOr(opts, "offset", offset.count()));
    valve_efficiency = GetOptionOr(opts, "efficiency", valve_efficiency);
    fuel_measurement_criterium = GetOptionOr(opts, "fuel_measurement_criterium",
                                             fuel_measurement_criterium);
    fuel_measurement_position = GetOptionOr(opts, "fuel_measurement_position",
                                            fuel_measurement_position);
}

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
  std::vector<double> moles(eq.GetReactor().GetNSpecies());
  std::vector<double> local_moles(moles.size());
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
      MPI_Allreduce(local_moles.data(), moles.data(), moles.size(), MPI_DOUBLE,
                    MPI_SUM, ::amrex::ParallelContext::CommunicatorAll());
      break;
    }
  }
  return moles;
}

template <typename GriddingAlgorithm>
int FindLevel(const ::amrex::Geometry& geom,
              const GriddingAlgorithm& gridding) {
  for (int level = 0; level < gridding.GetPatchHierarchy().GetNumberOfLevels();
       ++level) {
    if (geom.Domain() ==
        gridding.GetPatchHierarchy().GetGeometry(level).Domain()) {
      return level;
    }
  }
  return -1;
}

double ChangeState_(PressureValveState& state, const ::amrex::Geometry& geom,
                    const GriddingAlgorithm& grid, Duration& last_closed,
                    Duration& last_fuel_change,
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
  const int level = FindLevel(geom, grid);
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

int Sign_(double x) { return (x > 0) - (x < 0); }

void IsentropicExpansionWithoutDissipation_(IdealGasMix<1>& eq,
                                            Complete<IdealGasMix<1>>& dest,
                                            const Complete<IdealGasMix<1>>& src,
                                            double dest_pressure,
                                            double efficiency = 1.0) {
  double old_velocity = src.momentum[0] / src.density;
  const double h_before = src.energy / src.density + src.pressure / src.density;
  eq.SetReactorStateFromComplete(src);
  eq.GetReactor().SetPressureIsentropic(dest_pressure);
  eq.CompleteFromReactor(dest);
  const double h_after =
      dest.energy / dest.density + dest.pressure / dest.density;
  const double enthalpyDifference = h_before - h_after;
  const double u_border =
      Sign_(enthalpyDifference) *
      std::sqrt(efficiency * std::abs(enthalpyDifference) * 2 +
                old_velocity * old_velocity);
  dest.momentum[0] = dest.density * u_border;
  dest.energy += 0.5 * u_border * dest.momentum[0];
}

std::string GetMolesString_(const PressureValveOptions& options,
                            PressureValveState state) {
  if (state == PressureValveState::open_fuel) {
    return fmt::format("N2:79,O2:21,H2:{}", options.equivalence_ratio * 42.0);
  }
  return "N2:79,O2:21";
}
} // namespace

void PressureValveBoundary::FillBoundary(::amrex::MultiFab& mf,
                                         const ::amrex::Geometry& geom,
                                         Duration dt,
                                         const GriddingAlgorithm& grid,
                                         Direction dir) {
  if (dir == Direction::X) {
    FillBoundary(mf, geom, dt, grid, dir);
  }
}

void PressureValveBoundary::FillBoundary(::amrex::MultiFab& mf,
                                         const ::amrex::Geometry& geom,
                                         Duration /* dt */,
                                         const GriddingAlgorithm& grid) {
  const int ngrow = mf.nGrow(0);
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  Complete<IdealGasMix<1>> state{equation_};
  if (boundaries.isEmpty()) {
    return;
  }

  // Change State Machine if neccessary
  const double mean_pressure =
      ChangeState_(shared_valve_->state, geom, grid, shared_valve_->last_closed,
                   shared_valve_->last_fuel, options_, equation_);

  ReflectiveBoundary closed(execution::openmp, equation_, Direction::X, 0);
  switch (shared_valve_->state) {
  case PressureValveState::closed:
    closed.FillBoundary(mf, geom);
    break;
  case PressureValveState::open_air:
  case PressureValveState::open_fuel:
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
          ForEachIndex(box_to_fill, [this, mean_pressure, &state,
                                     &states](std::ptrdiff_t i, auto...) {
            std::array<std::ptrdiff_t, 1> dest{i};
            equation_.GetReactor().SetMoleFractions(
                GetMolesString_(options_, shared_valve_->state));
            equation_.GetReactor().SetTemperature(options_.outer_temperature);
            equation_.GetReactor().SetPressure(options_.outer_pressure);
            equation_.CompleteFromReactor(state);
            IsentropicExpansionWithoutDissipation_(equation_, state, state,
                                                   mean_pressure,
                                                   options_.valve_efficiency);
            Store(states, state, dest);
          });
        }
      }
    });
  }
}
} // namespace fub::amrex
