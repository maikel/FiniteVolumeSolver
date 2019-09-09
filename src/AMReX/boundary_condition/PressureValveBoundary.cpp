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

#include <utility>

namespace fub::amrex {

namespace po = boost::program_options;

po::options_description PressureValveBoundary::GetProgramOptions() {
  po::options_description desc{};
  // clang-format off
  desc.add_options()
    ("fuel_moles", po::value<std::string>()->default_value("N2:79,O2:21,H2:42"), "The species mixture for the fuel")
    ("air_moles", po::value<std::string>()->default_value("N2:79,O2:21"), "The species mixture for the air buffer")
    ("outer_pressure", po::value<double>()->default_value(6 * 101325.0), "The mean pressure value for the outer side [Pa]")
    ("pressure_value_which_opens_boudnary", po::value<double>()->default_value(2 * 101325.0), "The mean pressure value in the tube which opens the boundary [Pa]")
    ("pressure_value_which_closes_boundary", po::value<double>()->default_value(6 * 101325.0), "The mean pressure value in the tube which closes the boundary [Pa]")
    ("air_buffer_length", po::value<double>()->default_value(1.0), "The air buffer length which will be used to flush the tube [m]")
    ("fuel_length", po::value<double>()->default_value(1.0), "The fuel buffer length, will be capped by tube length [m]")
    ("ignition_position", po::value<double>()->default_value(0.69), "The position in the tube where the detonation will be triggered [m]");
  // clang-format on
  return desc;
}

PressureValveBoundary::PressureValveBoundary(const IdealGasMix<1>& equation,
                                             PressureValveOptions options)
    : options_{std::move(options)}, equation_{equation},
      state_{PressureValveState::open_air} {}

namespace {
PressureValveOptions GetOptions_(const po::variables_map& map) {
  PressureValveOptions options;
  options.fuel_moles = map["fuel_moles"].as<std::string>();
  options.air_moles = map["air_moles"].as<std::string>();
  options.outer_pressure = map["outer_pressure"].as<double>();
  options.pressure_value_which_opens_boudnary =
      map["pressure_value_which_opens_boudnary"].as<double>();
  options.pressure_value_which_closes_boundary =
      map["pressure_value_which_closes_boundary"].as<double>();
  options.air_buffer_length = map["air_buffer_length"].as<double>();
  options.fuel_length = map["fuel_length"].as<double>();
  options.ignition_position = map["ignition_position"].as<double>();
  return options;
}
} // namespace

PressureValveBoundary::PressureValveBoundary(const IdealGasMix<1>& eq,
                                             const po::variables_map& map)
    : PressureValveBoundary(eq, GetOptions_(map)) {}

const PressureValveOptions& PressureValveBoundary::GetOptions() const noexcept {
  return options_;
}

namespace {
std::vector<double> AsMolesVector(const std::string& moles,
                                  IdealGasMix<1>& eq) {
  eq.GetReactor().SetMoleFractions(moles);
  span<const double> ms = eq.GetReactor().GetMoleFractions();
  std::vector<double> mole(ms.begin(), ms.end());
  return mole;
}

bool AlmostEqual(span<const double> x, span<const double> y)
{
  return std::equal(x.begin(), x.end(), y.begin(), [](double x, double y) {
    return std::abs(x - y) < 1e-5;
  });
}

double GetMeanPressure(const GriddingAlgorithm& grid, IdealGasMix<1>& eq) {
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

std::vector<double> GatherMoles(const GriddingAlgorithm& grid, double x,
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
      MPI_Allreduce(local_moles.data(), moles.data(), moles.size(), MPI_DOUBLE, MPI_SUM,
                    ::amrex::ParallelContext::CommunicatorAll());
      break;
    }
  }
  return moles;
}

void ChangeState(PressureValveState& state, ::amrex::MultiFab& data,
                 const ::amrex::Geometry& geom, const GriddingAlgorithm& grid,
                 const PressureValveOptions& options, IdealGasMix<1>& eq) {
  const double mean_pressure = GetMeanPressure(grid, eq);
  const double xlo = geom.ProbDomain().lo(0);
  switch (state) {
  case PressureValveState::closed:
    if (mean_pressure < options.pressure_value_which_opens_boudnary) {
      state = PressureValveState::open_air;
    }
    break;
  case PressureValveState::open_air:
  case PressureValveState::open_fuel:
  case PressureValveState::ignition:
    if (mean_pressure > options.pressure_value_which_closes_boundary) {
      state = PressureValveState::closed;
    }
    break;
  }
  if (state == PressureValveState::open_air) {
    const double x_air =
        std::min(xlo + options.air_buffer_length, geom.ProbDomain().hi(0) - 0.5 * geom.CellSize(0));
    const std::vector<double> moles_at_x_air = GatherMoles(grid, x_air, eq);
    const std::vector<double> moles_air = AsMolesVector(options.air_moles, eq);
    if (AlmostEqual(moles_at_x_air, moles_air)) {
      state = PressureValveState::open_fuel;
    }
  }
  if (state == PressureValveState::open_fuel) {
    const double x_fuel =
        std::min(xlo + options.fuel_length, geom.ProbDomain().hi(0) - 0.5 * geom.CellSize(0));
    const std::vector<double> moles_at_x_fuel = GatherMoles(grid, x_fuel, eq);
    const std::vector<double> moles_fuel = AsMolesVector(options.fuel_moles, eq);
    if (AlmostEqual(moles_at_x_fuel, moles_fuel)) {
      const double x_ignite = std::min(x_fuel, options.ignition_position);
      Complete<IdealGasMix<1>> complete(eq);
      ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
        ::amrex::FArrayBox& fab = data[mfi];
        auto view = MakeView<Complete<IdealGasMix<1>>>(fab, eq, mfi.tilebox());
        ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i) {
          const double x = geom.CellCenter(int(i), 0);
          if (x_ignite - 0.05 <= x && x < x_ignite) {
            Load(complete, view, {i});
            eq.SetReactorStateFromComplete(complete);
            const double d = std::clamp((x_ignite - x) / 0.05, 0.0, 1.0);
            const double T_old = complete.temperature;
            constexpr double T_lo = 300.0;
            constexpr double T_hi = 2000.0;
            const double T_ramp = d * T_lo + (1.0 - d) * T_hi;
            const double T = std::max(T_old, T_ramp);
            eq.GetReactor().SetTemperature(T);
            eq.CompleteFromReactor(complete,
                                   complete.momentum / complete.density);
            Store(view, complete, {i});
          } else if (x_ignite <= x && x < x_ignite + 0.05) {
            Load(complete, view, {i});
            eq.SetReactorStateFromComplete(complete);
            const double d = std::clamp((x - x_ignite) / 0.05, 0.0, 1.0);
            const double T_old = complete.temperature;
            constexpr double T_lo = 300.0;
            constexpr double T_hi = 2000.0;
            const double T_ramp = d * T_lo + (1.0 - d) * T_hi;
            const double T = std::max(T_old, T_ramp);
            eq.GetReactor().SetTemperature(T);
            eq.CompleteFromReactor(complete,
                                   complete.momentum / complete.density);
            Store(view, complete, {i});
          }
        });
      });
      state = PressureValveState::ignition;
    }
  }
}

int Sign(double x) { return (x > 0) - (x < 0); }

void IsentropicExpansionWithoutDissipation(IdealGasMix<1>& eq,
                                           Complete<IdealGasMix<1>>& dest,
                                           const Complete<IdealGasMix<1>>& src,
                                           double dest_pressure,
                                           const std::string moles,
                                           double efficiency = 1.0) {
  double old_velocity = src.momentum[0] / src.density;
  eq.SetReactorStateFromComplete(src);
  eq.GetReactor().SetMoleFractions(moles);
  eq.GetReactor().SetTemperature(src.temperature);
  eq.CompleteFromReactor(dest);
  const double h_before =
      dest.energy / dest.density + dest.pressure / dest.density;
  eq.GetReactor().SetPressureIsentropic(dest_pressure);
  eq.CompleteFromReactor(dest);
  const double h_after =
      dest.energy / dest.density + dest.pressure / dest.density;
  const double enthalpyDifference = h_before - h_after;
  const double u_border =
      Sign(enthalpyDifference) *
      std::sqrt(efficiency * std::abs(enthalpyDifference) * 2 +
                old_velocity * old_velocity);
  dest.momentum[0] = dest.density * u_border;
  dest.energy += 0.5 * u_border * dest.momentum[0];
}

const std::string& GetMolesString(const PressureValveOptions& options,
                                  PressureValveState state) {
  if (state == PressureValveState::open_fuel ||
      state == PressureValveState::ignition) {
    return options.fuel_moles;
  }
  return options.air_moles;
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
  ChangeState(state_, mf, geom, grid, options_, equation_);

  ReflectiveBoundary closed(execution::openmp, equation_, Direction::X, 0);
  switch (state_) {
  case PressureValveState::closed:
    closed.FillBoundary(mf, geom);
    break;
  case PressureValveState::open_air:
  case PressureValveState::ignition:
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
          ForEachIndex(box_to_fill,
                       [this, &state, &states](std::ptrdiff_t i, auto...) {
                         std::array<std::ptrdiff_t, 1> dest{i};
                         std::array<std::ptrdiff_t, 1> src{0};
                         Load(state, states, src);
                         IsentropicExpansionWithoutDissipation(
                             equation_, state, state, options_.outer_pressure,
                             GetMolesString(options_, state_));
                         Store(states, state, dest);
                       });
        }
      }
    });
  }
}
} // namespace fub::amrex
