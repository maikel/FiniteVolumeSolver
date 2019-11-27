//
// Created by Patrick Denzler on 20.11.19.
//

#ifndef FUB_COUNTER_TIMER_HPP
#define FUB_COUNTER_TIMER_HPP

#include <string>
#include <chrono>
#include <functional>

#include "Counter.hpp"

namespace fub {

class Timer {
public:
  explicit Timer(Counter& counter);
  Timer(const Timer&) = delete;
  Timer& operator=(const Timer&) = delete;

  std::string const& get_name();

  ~Timer();

private:
  Counter* counter_;
  std::chrono::steady_clock::time_point start_time_;
  std::chrono::steady_clock::time_point stop_time_;

  void start();
  void stop();
  void submit();
};

} // namespace fub

#endif // FUB_COUNTER_TIMER_HPP