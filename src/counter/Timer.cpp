//
// Created by Patrick Denzler on 20.11.19.
//

#include "fub/counter/Timer.hpp"

namespace fub {

Timer::Timer(Counter& counter) : counter_(&counter) { this->start(); }

Timer::~Timer() {
  this->stop();
  this->submit();
}

std::string const& Timer::get_name() { return counter_->get_name(); }

void Timer::start() { start_time_ = std::chrono::steady_clock::now(); }

void Timer::stop() { stop_time_ = std::chrono::steady_clock::now(); }

void Timer::submit() {
  auto duration = stop_time_ - start_time_;
  auto nano = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
  counter_->add(nano.count());
}

} // namespace fub