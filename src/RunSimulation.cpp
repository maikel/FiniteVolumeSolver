#include "fub/RunSimulation.hpp"

#include <fmt/format.h>

#include <cmath>

namespace fub {

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
  return fmt::format("[{:3}%] {:05}: T = {:>9.6}s dt = {:6e}s wall-time: {:02}h:"
                     "{:02}m:{:6f}s (+{:6e}s)\n",
                     total_progress, cycle, simulation_time.count(),
                     time_step_size.count(), wt_hours, wt_minutes, wt_seconds,
                     wt_diff_seconds);
}

} // namespace fub