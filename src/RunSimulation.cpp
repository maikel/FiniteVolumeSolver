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

#include "fub/RunSimulation.hpp"

#include <fmt/format.h>
#include <fmt/chrono.h>

#include <cmath>

namespace fub {

boost::program_options::options_description GetCommandLineRunOptions() {
  namespace po = boost::program_options;
  RunOptions opts{};
  po::options_description desc{};
  // clang-format off
  desc.add_options()
    ("help", "Print help messages")
    ("max_cycles", po::value<int>()->default_value(opts.max_cycles), "Set maximal number of coarse time steps done.")
    ("cfl", po::value<double>()->default_value(opts.cfl), "Set the CFL condition")
    ("output_interval", po::value<double>()->default_value(0.0), "Sets the output interval")
    ("output_frequency", po::value<int>()->default_value(0), "Sets the output frequency")
    ("final_time", po::value<double>()->default_value(0.0), "Sets the final time point for this simulation");
  // clang-format on
  return desc;
}

RunOptions GetRunOptions(const boost::program_options::variables_map& vm) {
  RunOptions opts{};
  opts.cfl = vm["cfl"].as<double>();
  opts.final_time = Duration(vm["final_time"].as<double>());
  opts.output_interval = {Duration(vm["output_interval"].as<double>())};
  opts.output_frequency = {vm["output_frequency"].as<int>()};
  opts.max_cycles = vm["max_cycles"].as<int>();
  return opts;
}

void PrintRunOptions(std::ostream& out, const RunOptions& opts) {
  out << fmt::format(
      R"fmt([Info] ========================================================================
[Info] Run Options:
[Info] final_time = {}s
[Info] max_cycles = {}
[Info] output_interval = {}s
[Info] output_frequency = {}
[Info] cfl = {}
[Info] ========================================================================
)fmt",
      opts.final_time.count(), opts.max_cycles, opts.output_interval[0].count(),
      fmt::join(opts.output_frequency, ", "), opts.cfl);
}

std::string
FormatTimeStepLine(std::ptrdiff_t cycle,
                   std::chrono::duration<double> simulation_time,
                   std::chrono::duration<double> time_step_size,
                   std::chrono::duration<double> final_time,
                   std::chrono::steady_clock::duration wall_time,
                   std::chrono::steady_clock::duration wall_time_difference) {
  using namespace std::chrono;
  const int total_progress =
      static_cast<int>(simulation_time / final_time * 100.0);
  const auto wt_hours = duration_cast<hours>(wall_time).count();
  const auto wt_minutes = duration_cast<minutes>(wall_time).count() % 60;
  const double wt_seconds =
      std::fmod(duration_cast<duration<double>>(wall_time).count(), 60.);
  const double wt_diff_seconds =
      duration_cast<duration<double>>(wall_time_difference).count();
  return fmt::format(
      "[{:3}%] {:05}: T = {:>9.6}s dt = {:6e}s wall-time: {:02}h:"
      "{:02}m:{:6f}s (+{:6e}s)\n",
      total_progress, cycle, simulation_time.count(), time_step_size.count(),
      wt_hours, wt_minutes, wt_seconds, wt_diff_seconds);
}

std::optional<int> AnyOutputCondition(std::ptrdiff_t cycle, Duration time_point,
                                      const RunOptions& options) {
  const fub::Duration eps = options.smallest_time_step_size;
  if (!(time_point + eps < options.final_time &&
        (options.max_cycles < 0 || cycle < options.max_cycles))) {
    return -1;
  }
  for (std::size_t i = 0; i < options.output_interval.size(); ++i) {
    if (options.output_interval[i].count() > 0.0 &&
        std::fmod(time_point.count(), options.output_interval[i].count()) <
            2 * eps.count()) {
      return static_cast<int>(i);
    }
  }
  for (std::size_t i = 0; i < options.output_frequency.size(); ++i) {
    if (options.output_frequency[i] > 0 &&
        cycle % options.output_frequency[i] == 0) {
      return static_cast<int>(i);
    }
  }
  return {};
}

Duration NextOutputTime(Duration time_point, const RunOptions& options) {
  Duration next_output_time(std::numeric_limits<double>::infinity());
  for (std::size_t i = 0; i < options.output_interval.size(); ++i) {
    if (options.output_interval[i].count() > 0.0) {
      Duration dt(
          options.output_interval[i].count() -
          std::fmod(time_point.count(), options.output_interval[i].count()));
      dt += options.smallest_time_step_size;
      next_output_time = std::min(next_output_time, time_point + dt);
    }
  }
  return next_output_time;
}

} // namespace fub
