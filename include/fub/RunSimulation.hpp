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
#include "fub/core/assert.hpp"

#include <iostream>

namespace fub {

struct RunOptions {
  std::chrono::duration<double> final_time;
  std::chrono::duration<double> output_interval{0.1};
  std::chrono::duration<double> smallest_time_step_size{1e-8};
  int regrid_frequency{2};
};

std::string
FormatTimeStepLine(std::ptrdiff_t cycle,
                   std::chrono::duration<double> time_point,
                   std::chrono::duration<double> time_step_size,
                   std::chrono::duration<double> final_time,
                   std::chrono::steady_clock::duration wall_time,
                   std::chrono::steady_clock::duration wall_time_difference);

template <typename Solver, typename Output, typename Print>
void RunSimulation(Solver& solver, RunOptions options,
                   std::chrono::steady_clock::time_point wall_time_reference,
                   Output output, Print print) {
  fub::Duration time_point = solver.GetTimePoint();
  const fub::Duration eps = options.smallest_time_step_size;
  fub::Duration next_output_time = options.output_interval;
  std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
  std::chrono::steady_clock::duration wall_time = wall_time_reference - now;
  while (time_point + eps < options.final_time) {
    while (time_point + eps < next_output_time) {
      const fub::Duration stable_dt = solver.ComputeStableDt();
      FUB_ASSERT(stable_dt > eps);
      const fub::Duration limited_dt =
          std::min(next_output_time - time_point, stable_dt);
      solver.AdvanceHierarchy(limited_dt);
      time_point = solver.GetTimePoint();
      now = std::chrono::steady_clock::now();
      std::chrono::steady_clock::duration wall_time_difference =
          (now - wall_time_reference) - wall_time;
      wall_time = now - wall_time_reference;
      std::ptrdiff_t cycle = solver.GetCycles();
      print(FormatTimeStepLine(cycle, time_point, limited_dt,
                               options.final_time, wall_time,
                               wall_time_difference));
    }
    output(solver.GetPatchHierarchy(), solver.GetCycles(), solver.GetTimePoint());
    next_output_time += options.output_interval;
  }
}

} // namespace fub

#endif