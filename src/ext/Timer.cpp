//
// Created by Patrick Denzler on 10.10.19.
//

#include "fub/ext/Timer.hpp"

namespace fub {
Timer::Timer(String name) : name_(name) {}

void Timer::Start() {
  start_time_ = std::chrono::steady_clock::now();
}

void Timer::Stop() {
  stop_time_ = std::chrono::steady_clock::now();
  // Log
}

void Timer::Print() {
  if (start_time_.time_since_epoch() == std::chrono::steady_clock::duration(0)) {
    std::cout << "No measurement started!\n";
  } else if (stop_time_.time_since_epoch() == std::chrono::steady_clock::duration(0)) {
    std::cout << "No measurement completed!\n";
  } else {
    auto duration = stop_time_ - start_time_;
    auto nano = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
    std::cout << "Elapsed Time: " << nano.count() << " ns.\n";
  }
}

} // namespace fub