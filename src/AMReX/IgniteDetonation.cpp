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
#include "fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp"

namespace fub::amrex {
namespace po = boost::program_options;

po::options_description IgniteDetonationOptions::GetCommandLineOptions() {
  IgniteDetonationOptions opts{};
  po::options_description desc{};
  desc.add_options()
  ("measurement_position", po::value<double>()->default_value(opts.measurement_position), "Position within the tube which measures the fuel equivalence ratio. This position effectively determines the fill fraction. [m]")
  ("equivalence_ratio_criterium", po::value<double>()->default_value(opts.equivalence_ratio_criterium), "Ignites the fuel if an equivalence ratio is measured which is greater than this value. [-]")
  ("ignite_position", po::value<double>()->default_value(opts.ignite_position), "The position within the tube where the temperature will be risen to ignite the tube. This position should be smaller than `measure_position`. [m]")
  ("ignite_interval", po::value<double>()->default_value(opts.ignite_interval.count()), "The time interval in which no second ignition can be triggered. [s]");
  return desc;
}

IgniteDetonationOptions::IgniteDetonationOptions(const po::variables_map& vm) {
  measurement_position = vm["measurement_position"].as<double>();
  equivalence_ratio_criterium = vm["equivalence_ratio_criterium"].as<double>();
  ignite_position = vm["ignite_position"].as<double>();
  ignite_interval = Duration(vm["ignite_interval"].as<double>());
}

IgniteDetonation::IgniteDetonation(
    const fub::IdealGasMix<1>& eq, std::shared_ptr<GriddingAlgorithm> grid,
    const fub::amrex::IgniteDetonationOptions& opts)
    : equation_(eq), gridding_(std::move(grid)), options_{opts} {}

Duration IgniteDetonation::ComputeStableDt() const noexcept {
  return Duration(std::numeric_limits<double>::infinity());
}

namespace {

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
} // namespace

Result<void, TimeStepTooLarge>
IgniteDetonation::AdvanceLevel(int level, fub::Duration /* dt */) {
  PatchHierarchy& hier = gridding_->GetPatchHierarchy();
  Duration current_time = hier.GetTimePoint(level);
  Duration next_ignition_time = last_ignition_ + options_.ignite_interval;
  if (next_ignition_time < current_time) {
    std::vector<double> moles =
        GatherMoles_(*gridding_, options_.measurement_position, equation_);
    const double equivalence_ratio = moles[Burke2012::sO2] ?
        0.5 * moles[Burke2012::sH2] / moles[Burke2012::sO2] : 0.0;
    if (equivalence_ratio > options_.equivalence_ratio_criterium) {
      ::amrex::MultiFab& data = hier.GetPatchLevel(level).data;
      const ::amrex::Geometry& geom = hier.GetGeometry(level);
      const double x_ignite = options_.ignite_position;
      const double T_lo = options_.temperature_low;
      const double T_hi = options_.temperature_high;
      const double width = options_.ramp_width;
      IdealGasMix<1>::Complete state(equation_);
      ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
        ::amrex::FArrayBox& fab = data[mfi];
        auto view =
            MakeView<Complete<IdealGasMix<1>>>(fab, equation_, mfi.tilebox());
        ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i) {
          const double x = geom.CellCenter(i, 0);
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
      last_ignition_ = current_time;
      ::amrex::Print() << fmt::format("[Info][t = {}] Detonation Ignited!\n", last_ignition_.count());
    }
  }
  return boost::outcome_v2::success();
}

} // namespace fub::amrex
