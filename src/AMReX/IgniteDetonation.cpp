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

#include "fub/AMReX/IgniteDetonation.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/ForEach.hpp"
#include "fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp"

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

namespace fub::amrex {

IgniteDetonationOptions::IgniteDetonationOptions(
    const ProgramOptions& options) {
  FUB_GET_OPTION_VAR(options, channel);
  FUB_GET_OPTION_VAR(options, measurement_position);
  FUB_GET_OPTION_VAR(options, equivalence_ratio_criterium);
  FUB_GET_OPTION_VAR(options, ignite_position);
  FUB_GET_OPTION_VAR(options, offset);
  FUB_GET_OPTION_VAR(options, ignite_interval);
  FUB_GET_OPTION_VAR(options, temperature_low);
  FUB_GET_OPTION_VAR(options, temperature_high);
  FUB_GET_OPTION_VAR(options, ramp_width);
}

void IgniteDetonationOptions::Print(SeverityLogger& log) const {
  FUB_PRINT_OPTION_VAR(log, channel);
  FUB_PRINT_OPTION_VAR(log, measurement_position);
  FUB_PRINT_OPTION_VAR(log, equivalence_ratio_criterium);
  FUB_PRINT_OPTION_VAR(log, ignite_position);
  FUB_PRINT_OPTION_VAR(log, offset);
  FUB_PRINT_OPTION_VAR(log, ignite_interval);
  FUB_PRINT_OPTION_VAR(log, temperature_low);
  FUB_PRINT_OPTION_VAR(log, temperature_high);
  FUB_PRINT_OPTION_VAR(log, ramp_width);
}

IgniteDetonation::IgniteDetonation(
    const fub::IdealGasMix<1>& eq, int max_number_levels,
    const fub::amrex::IgniteDetonationOptions& opts)
    : equation_(eq), options_{opts},
      last_ignition_time_(static_cast<std::size_t>(max_number_levels),
                          Duration(std::numeric_limits<double>::min())) {}

Duration IgniteDetonation::ComputeStableDt(const IntegratorContext&, int) const
    noexcept {
  return Duration(std::numeric_limits<double>::max());
}

namespace {

std::vector<double> GatherMoles_(const IntegratorContext& grid, double x,
                                 IdealGasMix<1>& eq) {
  const PatchHierarchy& hier = grid.GetPatchHierarchy();
  const int nlevel = hier.GetNumberOfLevels();
  std::vector<double> moles(
      static_cast<std::size_t>(eq.GetReactor().GetNSpecies()));
  std::vector<double> local_moles(moles.size());
  for (int level = nlevel - 1; level >= 0; --level) {
    const ::amrex::MultiFab& data = grid.GetScratch(level);
    const ::amrex::Geometry& geom = grid.GetGeometry(level);
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
    // TODO: Specify on which communicator shall be acted
    // grid.GetMpiCommunicator() ?
    int found = 0;
    MPI_Allreduce(&local_found, &found, 1, MPI_INT, MPI_SUM,
                  ::amrex::ParallelContext::CommunicatorAll());
    if (found) {
      MPI_Allreduce(local_moles.data(), moles.data(),
                    static_cast<int>(moles.size()), MPI_DOUBLE, MPI_SUM,
                    ::amrex::ParallelContext::CommunicatorAll());
      break;
    }
  }
  return moles;
}
} // namespace

Duration IgniteDetonation::GetNextIgnitionTimePoint(int level) const noexcept {
  const std::size_t l = static_cast<std::size_t>(level);
  if (last_ignition_time_[l] < options_.offset) {
    return options_.offset;
  } else {
    return last_ignition_time_[l] + options_.ignite_interval;
  }
}

void IgniteDetonation::SetLastIgnitionTimePoint(int level,
                                                Duration t) noexcept {
  const std::size_t l = static_cast<std::size_t>(level);
  last_ignition_time_[l] = t;
}

Result<void, TimeStepTooLarge>
IgniteDetonation::AdvanceLevel(IntegratorContext& grid, int level,
                               fub::Duration /* dt */,
                               const ::amrex::IntVect& /* ngrow */) {
  Duration current_time = grid.GetTimePoint(level);
  Duration next_ignition_time = GetNextIgnitionTimePoint(level);
  if (next_ignition_time < current_time) {
    const ::amrex::Geometry& geom = grid.GetGeometry(level);
    const double dx_2 = geom.CellSize(0);
    const double xlo = geom.ProbDomain().lo(0) + dx_2;
    const double xhi = geom.ProbDomain().hi(0) - dx_2;
    const double x_fuel = std::clamp(options_.measurement_position, xlo, xhi);
    std::vector<double> moles = GatherMoles_(grid, x_fuel, equation_);
    const double equivalence_ratio =
        moles[Burke2012::sO2]
            ? 0.5 * moles[Burke2012::sH2] / moles[Burke2012::sO2]
            : 0.0;
    if (equivalence_ratio > options_.equivalence_ratio_criterium) {
      ::amrex::MultiFab& data = grid.GetScratch(level);
      const ::amrex::Geometry& geom = grid.GetGeometry(level);
      const double x_ignite = std::clamp(options_.ignite_position, xlo, xhi);
      const double T_lo = options_.temperature_low;
      const double T_hi = options_.temperature_high;
      const double width = options_.ramp_width;
      IdealGasMix<1>::Complete state(equation_);
      ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
        ::amrex::FArrayBox& fab = data[mfi];
        auto view =
            MakeView<Complete<IdealGasMix<1>>>(fab, equation_, mfi.tilebox());
        ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i) {
          const double x = geom.CellCenter(static_cast<int>(i), 0);
          if (x_ignite - width <= x && x < x_ignite) {
            const double d = (x_ignite - x) / width;
            const double T = d * T_lo + (1.0 - d) * T_hi;
            Load(state, view, {i});
            equation_.SetReactorStateFromComplete(state);
            equation_.GetReactor().SetTemperature(
                std::max(state.temperature, T));
            equation_.CompleteFromReactor(state);
            Store(view, state, {i});
          } else if (x_ignite <= x && x < x_ignite + width) {
            const double d = (x - x_ignite) / width;
            const double T = d * T_lo + (1.0 - d) * T_hi;
            Load(state, view, {i});
            equation_.SetReactorStateFromComplete(state);
            equation_.GetReactor().SetTemperature(
                std::max(state.temperature, T));
            equation_.CompleteFromReactor(state);
            Store(view, state, {i});
          }
        });
      });
      SetLastIgnitionTimePoint(level, next_ignition_time);
      fub::SeverityLogger log = GetInfoLogger();
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", options_.channel);
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", current_time.count());
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Level", level);
      BOOST_LOG(log) << "Detonation triggered at " << options_.ignite_position
                     << " [m], because equivalence ratio at " << x_fuel
                     << " [m] is " << equivalence_ratio << " [-]!";
    } else {
      fub::SeverityLogger debug = GetLogger(boost::log::trivial::debug);
      BOOST_LOG_SCOPED_LOGGER_TAG(debug, "Channel", options_.channel);
      BOOST_LOG_SCOPED_LOGGER_TAG(debug, "Time", current_time.count());
      BOOST_LOG_SCOPED_LOGGER_TAG(debug, "Level", level);
      BOOST_LOG(debug) << "Detonation did not trigger at "
                       << options_.ignite_position
                       << " [m], because equivalence ratio at " << x_fuel
                       << " [m] is " << equivalence_ratio << " [-]!";
    }
  }
  return boost::outcome_v2::success();
}

} // namespace fub::amrex
