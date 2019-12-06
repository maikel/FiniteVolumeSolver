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

template <typename Range>
struct IsMapType
    : std::is_convertible<decltype(*std::declval<Range>().begin()),
                          std::pair<const std::string, CounterResult>> {};

template <typename ResultRange>
int LengthOfFirstColumn(const ResultRange& results) {
  if (results.size()) {
    if constexpr (IsMapType<ResultRange>()) {
      return std::max_element(results.begin(), results.end(),
                              [](auto&& p1, auto&& p2) {
                                const std::string& name1 = std::get<0>(p1);
                                const std::string& name2 = std::get<0>(p2);
                                return name1.size() < name2.size();
                              })
          ->first.size();
    } else {
      return std::max_element(
                 results.begin(), results.end(),
                 [](const CounterResult& r1, const CounterResult& r2) {
                   return r1.name.size() < r2.name.size();
                 })
          ->name.size();
    }
  }
  return std::string_view("      Counter Name").size();
}

template <typename Duration, typename ResultRange>
void print_statistics(const ResultRange& results, long long reference = -1) {
  const char* last_column_title = [](long long r) {
    return r == -1 ? "Total [{}]" : "Percent [%]";
  }(reference);

  std::string header_format = fmt::format(
      "{{:>{}}} {{:>15}} {{:>15}} {{:>15}} {{:>15}} {{:>15}} {{:>15}}\n",
      LengthOfFirstColumn(results));
  std::string row_format = fmt::format(
      "{{:>{}}} {{:>15}} {{:>15}} {{:>15}} {{:>15}} {{:>15}} {{:>15.2}}\n",
      LengthOfFirstColumn(results));

  fmt::print("{:-<155}\n", "Table of Counters ");
  fmt::print(header_format.c_str(), "Counter Name", "#Calls",
             fmt::format("Min [{}]", DurationName<Duration>::name),
             fmt::format("Max [{}]", DurationName<Duration>::name),
             fmt::format("Avg [{}]", DurationName<Duration>::name),
             fmt::format("StdDev [{}]", DurationName<Duration>::name),
             fmt::format(last_column_title, DurationName<Duration>::name));
  if (results.empty()) {
    fmt::print("Registry does not contain any counters!\n");
  }
  using namespace std::chrono;
  if constexpr (IsMapType<ResultRange>()) {
    for (auto&& [name, result] : results) {
      double last_column_value =
          (reference == -1
               ? duration_cast<Duration>(
                     nanoseconds(std::llround(result.count * result.mean)))
                     .count()
               : ((result.count * result.mean) / reference) * 100);
      fmt::print(
          row_format.c_str(), name, result.count,
          duration_cast<Duration>(nanoseconds(result.min)).count(),
          duration_cast<Duration>(nanoseconds(result.max)).count(),
          duration_cast<Duration>(nanoseconds(std::llround(result.mean)))
              .count(),
          duration_cast<Duration>(nanoseconds(std::llround(result.stddev)))
              .count(),
          last_column_value);
    }
  } else {
    ResultRange copy = results;
    std::sort(copy.begin(), copy.end(),
              [](const CounterResult& r1, const CounterResult& r2) {
                return r1.count * r1.mean > r2.count * r2.mean;
              });
    for (auto&& result : copy) {
      double last_column_value =
          (reference == -1
               ? duration_cast<Duration>(
                     nanoseconds(std::llround(result.count * result.mean)))
                     .count()
               : ((result.count * result.mean) / reference) * 100);
      fmt::print(
          row_format.c_str(), result.name, result.count,
          duration_cast<Duration>(nanoseconds(result.min)).count(),
          duration_cast<Duration>(nanoseconds(result.max)).count(),
          duration_cast<Duration>(nanoseconds(std::llround(result.mean)))
              .count(),
          duration_cast<Duration>(nanoseconds(std::llround(result.stddev)))
              .count(),
          last_column_value);
    }
  }
  fmt::print("{:=<155}\n", "");
}

} // namespace fub

#endif /* FUB_COUNTER_COUNTER_REGISTRY_HPP */