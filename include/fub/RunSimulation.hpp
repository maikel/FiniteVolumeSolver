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

#ifndef FUB_SAMRAI_RUN_SIMULATION_HPP
#define FUB_SAMRAI_RUN_SIMULATION_HPP

#include "fub/Duration.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/core/assert.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include "fub/counter/CounterRegistry.hpp"
#include "fub/counter/Timer.hpp"

#include "fub/output/BasicOutput.hpp"

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>
#include <boost/outcome.hpp>

#include <fmt/format.h>

#include <optional>
#include <vector>

namespace fub {

struct RunOptions {
  RunOptions() = default;

  RunOptions(const ProgramOptions& vm);

  template <typename Logger> void Print(Logger& log) {
    BOOST_LOG(log) << " - final_time = " << final_time.count() << " [s]";
    BOOST_LOG(log) << " - max_cycles = " << max_cycles << " [-]";
    BOOST_LOG(log) << " - smallest_time_step_size = "
                   << smallest_time_step_size.count() << " [s]";
    BOOST_LOG(log) << " - cfl = " << cfl << " [-]";
    BOOST_LOG(log) << " - do_backup = " << do_backup << " [-]";
  }

  Duration final_time{1.0};
  std::ptrdiff_t max_cycles{-1};
  Duration smallest_time_step_size{1e-12};
  double cfl{0.8};
  int do_backup{1};
};

std::string
FormatTimeStepLine(std::ptrdiff_t cycle,
                   std::chrono::duration<double> time_point,
                   std::chrono::duration<double> time_step_size,
                   std::chrono::duration<double> final_time,
                   std::chrono::steady_clock::duration wall_time,
                   std::chrono::steady_clock::duration wall_time_difference);

template <typename Solver>
void MakeBackup(std::optional<Solver>& dest, const Solver& src,
                CounterRegistry& counter_database) {
  Timer counter = counter_database.get_timer("RunSimulation::MakeBackup");
  dest = src;
}

template <typename Solver,
          typename Grid = std::decay_t<
              decltype(*std::declval<Solver&>().GetGriddingAlgorithm())>>
void RunSimulation(Solver& solver, RunOptions options,
                   std::chrono::steady_clock::time_point wall_time_reference,
                   BasicOutput<Grid>& output) {
  namespace logger = boost::log;
  using namespace logger::trivial;
  logger::sources::severity_logger<severity_level> log(
      logger::keywords::severity = info);
  fub::Duration time_point = solver.GetTimePoint();
  const fub::Duration eps = options.smallest_time_step_size;
  std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
  std::chrono::steady_clock::duration wall_time = wall_time_reference - now;
  // Deeply copy the grid for a fallback scenario in case of time step errors
  std::optional<Solver> backup{};
  if (options.do_backup) {
    Timer counter = solver.GetCounterRegistry()->get_timer("RunSimulation::MakeBackup");
    backup = solver;
  }
  std::optional<Duration> failure_dt{};
  while (time_point + eps < options.final_time &&
         (options.max_cycles < 0 || solver.GetCycles() < options.max_cycles)) {
    // We have a nested loop to exactly reach output time points.
    do {
      const fub::Duration next_output_time = output.NextOutputTime(time_point);
      solver.PreAdvanceHierarchy();
      // Compute the next time step size. If an estimate is available from a
      // prior failure we use that one.
      const fub::Duration stable_dt =
          failure_dt ? *failure_dt : solver.ComputeStableDt();
      FUB_ASSERT(stable_dt > eps);
      const fub::Duration limited_dt =
          std::min({options.final_time - time_point,
                    next_output_time - time_point, options.cfl * stable_dt});

      // Advance the hierarchy in time!
      boost::outcome_v2::result<void, TimeStepTooLarge> result =
          solver.AdvanceHierarchy(limited_dt);

      if (result.has_error() && options.do_backup) {
        // If the solver returned with an error, reduce the time step size with
        // the new estimate.
        failure_dt = result.error().dt;
        BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", time_point.count());
        BOOST_LOG_SEV(log, warning) << fmt::format(
            "Pre-estimated coarse time step size (dt_old = {}s) "
            "was too large.\nRestarting this time step with a "
            "smaller coarse time step size (dt_new = {}s).\n",
            limited_dt.count(), options.cfl * failure_dt->count());
        solver = *backup;
      } else if (result.has_error() && !options.do_backup) {
        BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", time_point.count());
        BOOST_LOG_SEV(log, error) << fmt::format(
            "Pre-estimated coarse time step size (dt_old = {}s) was too large "
            "but no backup present. Aborting.",
            limited_dt.count());
        return;
      } else {
        solver.PostAdvanceHierarchy(limited_dt);
        // If advancing the hierarchy was successfull print a successful time
        // step line and reset any failure indicators.
        now = std::chrono::steady_clock::now();
        time_point = solver.GetTimePoint();
        std::chrono::steady_clock::duration wall_time_difference =
            (now - wall_time_reference) - wall_time;
        wall_time = now - wall_time_reference;
        const std::ptrdiff_t cycle = solver.GetCycles();
        BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", time_point.count());
        BOOST_LOG(log) << FormatTimeStepLine(cycle, time_point, limited_dt,
                                             options.final_time, wall_time,
                                             wall_time_difference);
        failure_dt.reset();
        if (options.do_backup) {
          Timer counter = solver.GetCounterRegistry()->get_timer("RunSimulation::MakeBackup");
          backup = solver;
        }
      }
    } while (
        time_point + eps < options.final_time &&
        (options.max_cycles < 0 || solver.GetCycles() < options.max_cycles) &&
        !output.ShallOutputNow(*solver.GetGriddingAlgorithm()));
    auto output_timer = solver.GetCounterRegistry()->get_timer("Output");
    output(*solver.GetGriddingAlgorithm());
  }
  // return solver;
}

} // namespace fub

#endif
