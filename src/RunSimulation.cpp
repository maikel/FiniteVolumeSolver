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

RunOptions::RunOptions(const std::map<std::string, pybind11::object>& vm) {
  final_time = Duration(GetOptionOr(vm, "final_time", final_time.count()));
  max_cycles = GetOptionOr(vm, "max_cycles", max_cycles);
  smallest_time_step_size = Duration(GetOptionOr(
      vm, "smallest_time_step_size", smallest_time_step_size.count()));
  cfl = GetOptionOr(vm, "cfl", cfl);
  do_backup = GetOptionOr(vm, "do_backup", do_backup);
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
  return fmt::format("[{:3}%] {:05}: dt = {:6e}s wall-time: {:02}h:"
                     "{:02}m:{:6f}s (+{:6e}s)",
                     total_progress, cycle, time_step_size.count(), wt_hours,
                     wt_minutes, wt_seconds, wt_diff_seconds);
}

} // namespace fub
