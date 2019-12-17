// Copyright (c) 2019 Patrick Denzler
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

#include "fub/counter/Timer.hpp"
#include "fub/core/assert.hpp"

namespace fub {

Timer::Timer(Counter& counter) : counter_(&counter) { this->start(); }

Timer::Timer(Timer&& other) noexcept
    : counter_(std::exchange(other.counter_, nullptr)),
      start_time_(std::exchange(other.start_time_,
                                std::chrono::steady_clock::time_point())),
      stop_time_(std::exchange(other.stop_time_,
                               std::chrono::steady_clock::time_point())) {}

Timer& Timer::operator=(Timer&& other) noexcept {
  counter_ = std::exchange(other.counter_, nullptr);
  start_time_ =
      std::exchange(other.start_time_, std::chrono::steady_clock::time_point());
  stop_time_ =
      std::exchange(other.stop_time_, std::chrono::steady_clock::time_point());
  return *this;
}

Timer::~Timer() {
  if (counter_) {
    this->stop();
    this->submit();
  }
}

std::string const& Timer::get_name() { return counter_->get_name(); }

void Timer::start() { start_time_ = std::chrono::steady_clock::now(); }

void Timer::stop() { stop_time_ = std::chrono::steady_clock::now(); }

void Timer::submit() {
  auto duration = stop_time_ - start_time_;
  auto nano = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
  static_assert(sizeof(long long) == sizeof(std::int64_t));
  static_assert(std::is_same_v<std::chrono::nanoseconds::rep, long long>);
  FUB_ASSERT(counter_);
  counter_->add(static_cast<long long>(nano.count()));
}

} // namespace fub