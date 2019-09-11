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

namespace fub::amrex {

namespace po = boost::program_options;

po::options_description
PressureValveOptions::GetCommandLineOptions(std::string prefix) {
  po::options_description desc{};
  PressureValveOptions opts{};
  if (prefix.empty()) {
    // clang-format off
    desc.add_options()
    ("outer_pressure", po::value<double>()->default_value(opts.outer_pressure), "The mean pressure value for the outer side [Pa]")
    ("outer_temperature", po::value<double>()->default_value(opts.outer_temperature), "The mean pressure value for the outer side [K]")
    ("pressure_value_which_opens_boundary", po::value<double>()->default_value(opts.pressure_value_which_opens_boundary), "The mean pressure value in the tube which opens the boundary. [Pa]")
    ("pressure_value_which_closes_boundary", po::value<double>()->default_value(opts.pressure_value_which_closes_boundary), "The mean pressure value in the tube which closes the boundary. [Pa]")
    ("oxygen_measurement_position", po::value<double>()->default_value(opts.oxygen_measurement_position), "The position within the tube where the oxygen concentration will be measured. [m]")
    ("oxygen_measurement_criterium", po::value<double>()->default_value(opts.oxygen_measurement_criterium), "The oxygen concentration which will trigger fuel instead of air inflow. [-]")
    ("equivalence_ratio", po::value<double>()->default_value(opts.equivalence_ratio), "The equivalence ratio of the fuel which will be used. [-]")
    ("open_at_interval", po::value<double>()->default_value(opts.open_at_interval.count()), "If set to non-zero value this pressure valve will only open up once in a specified time interval. [s]")
        ("valve_efficiency", po::value<double>()->default_value(opts.valve_efficiency), "Sets the efficiency of the enthalpy to velocity conversion [-]");
    // clang-format on
  } else {
    using fmt::format;
    std::string name = format("{}:outer_pressure", prefix);
    desc.add_options()(name.c_str(),
                       po::value<double>()->default_value(opts.outer_pressure),
                       "The mean pressure value for the outer side [Pa]");

    name = format("{}:outer_temperature", prefix);
    desc.add_options()(
        name.c_str(),
        po::value<double>()->default_value(opts.outer_temperature),
        "The mean pressure value for the outer side [K]");

    name = format("{}:pressure_value_which_opens_boundary", prefix);
    desc.add_options()(
        name.c_str(),
        po::value<double>()->default_value(
            opts.pressure_value_which_opens_boundary),
        "The mean pressure value in the tube which opens the boundary. [Pa]");

    name = format("{}:pressure_value_which_closes_boundary", prefix);
    desc.add_options()(
        name.c_str(),
        po::value<double>()->default_value(
            opts.pressure_value_which_closes_boundary),
        "The mean pressure value in the tube which closes the boundary. [Pa]");

    name = format("{}:oxygen_measurement_position", prefix);
    desc.add_options()(
        name.c_str(),
        po::value<double>()->default_value(opts.oxygen_measurement_position),
        "The position within the tube where the oxygen concentration will be "
        "measured. [m]");

    name = format("{}:oxygen_measurement_criterium", prefix);
    desc.add_options()(
        name.c_str(),
        po::value<double>()->default_value(opts.oxygen_measurement_criterium),
        "The oxygen concentration which will trigger fuel instead of air "
        "inflow. [-]");

    name = format("{}:equivalence_ratio", prefix);
    desc.add_options()(
        name.c_str(),
        po::value<double>()->default_value(opts.equivalence_ratio),
        "The equivalence ratio of the fuel which will be used. [-]");

    name = format("{}:open_interval", prefix);
    desc.add_options()(
        name.c_str(),
        po::value<double>()->default_value(opts.open_at_interval.count()),
        "If set to non-zero value this pressure valve will only open up once "
        "in a specified time interval. [s]");
  }
  return desc;
}

PressureValveBoundary::PressureValveBoundary(const IdealGasMix<1>& equation,
                                             PressureValveOptions options)
    : options_{std::move(options)}, equation_{equation},
      state_{PressureValveState::open_air} {}

PressureValveOptions::PressureValveOptions(const po::variables_map& map) {
  outer_pressure = map["outer_pressure"].as<double>();
  outer_temperature = map["outer_temperature"].as<double>();
  pressure_value_which_opens_boundary =
      map["pressure_value_which_opens_boundary"].as<double>();
  pressure_value_which_closes_boundary =
      map["pressure_value_which_closes_boundary"].as<double>();
  oxygen_measurement_position = map["oxygen_measurement_position"].as<double>();
  oxygen_measurement_criterium =
      map["oxygen_measurement_criterium"].as<double>();
  equivalence_ratio = map["equivalence_ratio"].as<double>();
  open_at_interval = Duration(map["open_at_interval"].as<double>());
  valve_efficiency = map["valve_efficiency"].as<double>();
}

PressureValveBoundary::PressureValveBoundary(const IdealGasMix<1>& eq,
                                             const po::variables_map& map)
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

double ChangeState_(PressureValveState& state, const ::amrex::Geometry& geom,
                    const GriddingAlgorithm& grid,
                    const PressureValveOptions& options, IdealGasMix<1>& eq) {
  const double mean_pressure = GetMeanPressure_(grid, eq);
  const double dx_2 = geom.CellSize(0);
  const double xlo = geom.ProbDomain().lo(0) + dx_2;
  const double xhi = geom.ProbDomain().hi(0) - dx_2;
  switch (state) {
  case PressureValveState::closed:
    if (mean_pressure < options.pressure_value_which_opens_boundary) {
      state = PressureValveState::open_air;
    }
    break;
  case PressureValveState::open_air:
  case PressureValveState::open_fuel:
    if (mean_pressure > options.pressure_value_which_closes_boundary) {
      state = PressureValveState::closed;
    }
    break;
  }
  if (state == PressureValveState::open_air) {
    const double x_air =
        std::clamp(options.oxygen_measurement_position, xlo, xhi);
    std::vector<double> moles = GatherMoles_(grid, x_air, eq);
    const double sum = std::accumulate(moles.begin(), moles.end(), 0.0);
    std::transform(moles.begin(), moles.end(), moles.begin(),
                   [sum](double m) { return m / sum; });
    if (options.oxygen_measurement_criterium < moles[Burke2012::sO2]) {
      state = PressureValveState::open_fuel;
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
      ChangeState_(state_, geom, grid, options_, equation_);

  ReflectiveBoundary closed(execution::openmp, equation_, Direction::X, 0);
  switch (state_) {
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
                GetMolesString_(options_, state_));
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
