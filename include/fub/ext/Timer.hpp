//
// Created by Patrick Denzler on 10.10.19.
//

#ifndef FUB_EXT_TIMER_HPP
#define FUB_EXT_TIMER_HPP

#include <string>
#include <chrono>
#include <functional>

namespace fub {

class Timer {
public:
  Timer(String name);
  void Start();
  void Stop();
  void Measure(F&& function, Args&&... args);
  void Print();

private:
  String name_;
  std::chrono::time_point start_time_;
  std::chrono::time_point stop_time_;
};

template<typename F, typename... Args>
void Timer::Measure(F&& function, Args&&... args) {
  this->Start();
  std::invoke(function, args...);
  this->Stop();
}

} // namespace fub

#endif // FUB_EXT_TIMER_HPP
