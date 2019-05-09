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

#ifndef FUB_AMREX_RUN_SIMULATION_HPP
#define FUB_AMREX_RUN_SIMULATION_HPP

#include "fub/Duration.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/core/assert.hpp"
#include "fub/ext/outcome.hpp"
#include "fub/RunSimulation.hpp"

#include <boost/optional.hpp>

#include <fmt/format.h>

namespace fub::amrex {

template <typename Plenum, typename Tube, typename Output, typename Print>
void RunCoupledSimulation(
    Plenum& plenum, Tube& tube, CoupledBoundary& boundary, RunOptions options,
    std::chrono::steady_clock::time_point wall_time_reference, Output output,
    Print print) {
  fub::Duration time_point = plenum.GetTimePoint();
  const fub::Duration eps = options.smallest_time_step_size;
  fub::Duration next_output_time = options.output_interval;
  std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
  std::chrono::steady_clock::duration wall_time = wall_time_reference - now;
  auto plenum_backup =
  std::make_shared<cutcell::GriddingAlgorithm>(*plenum.GetGriddingAlgorithm());
  auto tube_backup =
      std::make_shared<GriddingAlgorithm>(*tube.GetGriddingAlgorithm());
  boost::optional<Duration> failure_dt{};
  while (time_point + eps < options.final_time &&
         (options.max_cycles < 0 || plenum.GetCycles() < options.max_cycles)) {
    // We have a nested loop to exactly reach output time points.
    do {
      plenum.PreAdvanceHierarchy();
      tube.PreAdvanceHierarchy();

      // Compute the next time step size. If an estimate is available from a
      // prior failure we use that one.
      const fub::Duration stable_dt =
          failure_dt
              ? *failure_dt
              : std::min(tube.ComputeStableDt(), plenum.ComputeStableDt());
      FUB_ASSERT(stable_dt > eps);
      const fub::Duration limited_dt =
          std::min(next_output_time - time_point, options.cfl * stable_dt);

      // Advance the hierarchy in time!
      Result<void, TimeStepTooLarge> result =
          plenum.AdvanceHierarchy(limited_dt);
      if (result) {
        result = tube.AdvanceHierarchy(limited_dt);
      }

      if (result.has_error()) {
        // If the solver returned with an error, reduce the time step size with
        // the new estimate.
        failure_dt = result.error().coarse_dt;
        print(fmt::format("Pre-estimated coarse time step size (dt_old = {}s) "
                          "was too large.\nRestarting this time step with a "
                          "smaller coarse time step size (dt_new = {}s).\n",
                          limited_dt.count(),
                          options.cfl * failure_dt->count()));
        plenum.ResetHierarchyConfiguration(plenum_backup);
        tube.ResetHierarchyConfiguration(tube_backup);
      } else {
        plenum.PostAdvanceHierarchy();
        tube.PostAdvanceHierarchy();
        // If advancing the hierarchy was successfull print a successful time
        // step line and reset any failure indicators.
        now = std::chrono::steady_clock::now();
        time_point = plenum.GetTimePoint();
        std::chrono::steady_clock::duration wall_time_difference =
            (now - wall_time_reference) - wall_time;
        wall_time = now - wall_time_reference;
        const std::ptrdiff_t cycle = plenum.GetCycles();
        print(FormatTimeStepLine(cycle, time_point, limited_dt,
                                 options.final_time, wall_time,
                                 wall_time_difference));
        failure_dt.reset();
        plenum_backup =
        std::make_shared<cutcell::GriddingAlgorithm>(*plenum.GetGriddingAlgorithm());
        tube_backup =
            std::make_shared<GriddingAlgorithm>(*tube.GetGriddingAlgorithm());
        boundary.ComputeBoundaryData(plenum.GetPatchHierarchy(),
                                     tube.GetPatchHierarchy());
      }
    } while (time_point + eps < next_output_time &&
             (options.output_frequency <= 0 || failure_dt ||
              (plenum.GetCycles() % options.output_frequency) != 0));
    output(plenum.GetPatchHierarchy(), plenum.GetCycles(),
           plenum.GetTimePoint());
    next_output_time += options.output_interval;
  }
}

} // namespace fub::amrex

#endif
