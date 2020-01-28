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

#ifndef FUB_COUNTER_TIMER_HPP
#define FUB_COUNTER_TIMER_HPP

#include <chrono>
#include <functional>
#include <string>

#include "Counter.hpp"

namespace fub {

class Timer {
public:
  Timer() = default;
  explicit Timer(Counter& counter);
  Timer(const Timer&) = delete;
  Timer& operator=(const Timer&) = delete;
  Timer(Timer&&) noexcept;
  Timer& operator=(Timer&&) noexcept;

  std::string const& get_name();

  ~Timer();

private:
  Counter* counter_{nullptr};
  std::chrono::steady_clock::time_point start_time_{};
  std::chrono::steady_clock::time_point stop_time_{};

  void start();
  void stop();
  void submit();
};

} // namespace fub

#endif // FUB_COUNTER_TIMER_HPP