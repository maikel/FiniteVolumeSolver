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

#include <fmt/chrono.h>
#include <fmt/format.h>

#include <cmath>

namespace fub {

RunOptions::RunOptions(const boost::program_options::variables_map& vm) {
  auto GetOptionOr = [&](const char* opt, auto default_value) {
    if (vm.count(opt)) {
      return vm[opt].as<std::decay_t<decltype(default_value)>>();
    }
    return default_value;
  };
  final_time = Duration(GetOptionOr("final_time", 1.0));
  max_cycles = GetOptionOr("max_cycles", -1);
  std::vector<double> intervals = GetOptionOr("output.interval", std::vector{0.0});
  output_interval.resize(intervals.size());
  std::transform(intervals.begin(), intervals.end(), output_interval.begin(), [](double x) {
      return Duration(x);
  });
  output_frequency = GetOptionOr("output.frequency", std::vector{0});
  smallest_time_step_size = Duration(GetOptionOr("smallest_time_step_size", 1.0E-12));
  cfl = GetOptionOr("cfl", 0.8);
}

boost::program_options::options_description RunOptions::GetCommandLineOptions() {
  namespace po = boost::program_options;
  RunOptions opts{};
  po::options_description desc{"Run Options"};
  // clang-format off
  desc.add_options()
    ("help", "Print help messages")
    ("max_cycles", po::value<int>(), "Set maximal number of coarse time steps done.")
    ("cfl", po::value<double>(), "Set the CFL condition")
    ("output.interval", po::value<std::vector<double>>()->multitoken(), "Sets the output interval")
    ("output.frequency", po::value<std::vector<int>>()->multitoken(), "Sets the output frequency")
    ("final_time", po::value<double>(), "Sets the final time point for this simulation")
    ("smallest_time_step_size", po::value<double>(), "Set the minimal possible time step size");
  // clang-format on
  return desc;
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
      "[{:3}%] {:05}: dt = {:6e}s wall-time: {:02}h:"
      "{:02}m:{:6f}s (+{:6e}s)",
      total_progress, cycle, time_step_size.count(),
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
