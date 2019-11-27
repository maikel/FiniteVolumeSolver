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

#include <string>
#include <unordered_map>
#include <fmt/format.h>

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
 * Functions to print a table of all gathered counters (#call, min, max, avg, stddev and total or percent)
 */
template <typename T> void print_statistics(const T& result) {
  std::string format_string =
      "{:>55} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n";
  fmt::print("{:-<155}\n", "Table of Counters ");
  fmt::print(format_string, "Counter Name", "#Calls", "Min [ns]", "Max [ns]",
             "Avg [ns]", "StdDev [ns]", "Total [ns]");
  for (auto& c : result) {
    if constexpr (std::is_same_v<
                      T, std::unordered_map<std::string, CounterResult>>) {
      fmt::print(format_string, c.second.name, c.second.count, c.second.min,
                 c.second.max, c.second.mean, c.second.stddev,
                 c.second.count * c.second.mean);
    } else {
      fmt::print(format_string, c.name, c.count, c.min, c.max, c.mean, c.stddev,
                 c.count * c.mean);
    }
    // fmt::print(format_string, c.name, c.count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.min))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.mean))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.max))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.stddev))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.count*c.mean))).count());
  }
  fmt::print("{:=<155}\n", "");
}

template <typename T> void print_statistics(T result, long long reference) {
  fmt::print("{:-<155}\n", "Table of Counters ");
  fmt::print("{:>55} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n",
             "Counter Name", "#Calls", "Min [ns]", "Max [ns]", "Avg [ns]",
             "StdDev [ns]", "Percent [%]");
  for (auto& c : result) {
    if constexpr (std::is_same_v<
                      T, std::unordered_map<std::string, CounterResult>>) {
      fmt::print("{:>55} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15.2}\n",
                 c.second.name, c.second.count, c.second.min, c.second.max,
                 c.second.mean, c.second.stddev,
                 (c.second.count * c.second.mean) / reference);
    } else {
      fmt::print("{:>55} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15.2}\n",
                 c.name, c.count, c.min, c.max, c.mean, c.stddev,
                 (c.count * c.mean) / reference);
    }
    // fmt::print(format_string, c.name, c.count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.min))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.mean))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.max))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.stddev))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.count*c.mean))).count());
  }
  fmt::print("{:=<155}\n", "");
}

} // namespace fub

#endif /* FUB_COUNTER_COUNTER_REGISTRY_HPP */