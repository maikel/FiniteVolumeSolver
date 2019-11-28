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

#ifndef FUB_COUNTER_COUNTER_REGISTRY_HPP
#define FUB_COUNTER_COUNTER_REGISTRY_HPP

#include <fmt/format.h>
#include <map>
#include <string>
#include <unordered_map>

#include "Counter.hpp"
#include "Timer.hpp"

namespace fub {

class CounterRegistry {
public:
  CounterRegistry() = default;
  Counter* get(std::string const& name);
  Timer get_timer(std::string const& name);
  std::vector<CounterResult> gather_statistics();
  std::unordered_map<std::string, CounterResult> gather_statistics_map();

private:
  std::unordered_map<std::string, Counter> database_;
};

/*
 * Functions to print a table of all gathered counters (#call, min, max, avg,
 * stddev and total or percent)
 */
template <typename Duration> struct DurationName {
  static constexpr const char* name = "unknown";
};

template <> struct DurationName<std::chrono::nanoseconds> {
  static constexpr const char* name = "ns";
};

template <> struct DurationName<std::chrono::microseconds> {
  static constexpr const char* name = "us";
};

template <> struct DurationName<std::chrono::milliseconds> {
  static constexpr const char* name = "ms";
};

template <> struct DurationName<std::chrono::seconds> {
  static constexpr const char* name = "s";
};

template <typename Duration, typename T>
void print_statistics(T result, long long reference = -1) {
  const char* last_column_title = [](long long r) {
    return r == -1 ? "Total [{}]" : "Percent [%]";
  }(reference);

  fmt::print("{:-<155}\n", "Table of Counters ");
  fmt::print("{:>55} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n",
             "Counter Name", "#Calls",
             fmt::format("Min [{}]", DurationName<Duration>{}.name),
             fmt::format("Max [{}]", DurationName<Duration>{}.name),
             fmt::format("Avg [{}]", DurationName<Duration>{}.name),
             fmt::format("StdDev [{}]", DurationName<Duration>{}.name),
             fmt::format(last_column_title, DurationName<Duration>{}.name));
  if (result.empty()) {
    fmt::print("{}", "Registry does not countain counters!");
  }
  for (auto& c : result) {
    if constexpr (std::is_same_v<decltype(c), std::pair<const std::string,
                                                        CounterResult>&>) {
      double last_column_value =
          (reference == -1 ? std::chrono::duration_cast<Duration>(
                                 std::chrono::nanoseconds(std::llround(
                                     c.second.count * c.second.mean)))
                                 .count()
                           : ((c.second.count * c.second.mean) / reference)*100);
      fmt::print("{:>55} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15.2}\n",
                 c.second.name, c.second.count,
                 std::chrono::duration_cast<Duration>(
                     std::chrono::nanoseconds(std::llround(c.second.min)))
                     .count(),
                 std::chrono::duration_cast<Duration>(
                     std::chrono::nanoseconds(std::llround(c.second.max)))
                     .count(),
                 std::chrono::duration_cast<Duration>(
                     std::chrono::nanoseconds(std::llround(c.second.mean)))
                     .count(),
                 std::chrono::duration_cast<Duration>(
                     std::chrono::nanoseconds(std::llround(c.second.stddev)))
                     .count(),
                 last_column_value);
    } else if constexpr (std::is_same_v<decltype(c), CounterResult&>) {
      double last_column_value =
          (reference == -1
               ? std::chrono::duration_cast<Duration>(
                     std::chrono::nanoseconds(std::llround(c.count * c.mean)))
                     .count()
               : ((c.count * c.mean) / reference)*100);
      fmt::print("{:>55} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15.2}\n", c.name,
                 c.count,
                 std::chrono::duration_cast<Duration>(std::chrono::nanoseconds(std::llround(c.min))).count(),
                 std::chrono::duration_cast<Duration>(std::chrono::nanoseconds(std::llround(c.max))).count(),
                 std::chrono::duration_cast<Duration>(std::chrono::nanoseconds(std::llround(c.mean))).count(),
                 std::chrono::duration_cast<Duration>(std::chrono::nanoseconds(std::llround(c.stddev))).count(),
                 last_column_value);
    } else {
      fmt::print("{}", "Result type not supported!\n");
    }
    // fmt::print(format_string, c.name, c.count,
    // std::chrono::duration_cast<Duration>(std::chrono::nanoseconds(std::llround(c.min))).count(),
    // std::chrono::duration_cast<Duration>(std::chrono::nanoseconds(std::llround(c.mean))).count(),
    // std::chrono::duration_cast<Duration>(std::chrono::nanoseconds(std::llround(c.max))).count(),
    // std::chrono::duration_cast<Duration>(std::chrono::nanoseconds(std::llround(c.stddev))).count(),
    // std::chrono::duration_cast<Duration>(std::chrono::nanoseconds(std::llround(c.count*c.mean))).count());
  }
  fmt::print("{:=<155}\n", "");
}

} // namespace fub

#endif /* FUB_COUNTER_COUNTER_REGISTRY_HPP */