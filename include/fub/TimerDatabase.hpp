// Copyright (c) 2019 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef FUB_TIMERDATABASE_HPP
#define FUB_TIMERDATABASE_HPP

#include "fub/core/span.hpp"

#include <chrono>
#include <string>
#include <vector>

namespace fub::profile {

class TimerDatabase {
public:
  using TimePoint = std::chrono::steady_clock::time_point;
  using Counter = std::chrono::steady_clock::duration;

  int RegisterTimer(std::string timer);

  void TraceStart(std::string id, std::chrono::steady_clock::duration count);
  void TraceStop(std::string id, std::chrono::steady_clock::duration count);

  span<const Counter> GetCounters(int id) const noexcept { return id_to_counters_[id]; }
  

private:
  std::vector<std::string> id_to_names_;
  std::vector<std::vector<Counter>> id_to_counters_;
  std::vector<int> event_to_id_;
  std::vector<TimePoint> start_events_;
  std::vector<TimePoint> stop_events_;
};

}

#endif

