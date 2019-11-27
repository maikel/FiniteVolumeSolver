//
// Created by Patrick Denzler on 07.11.19.
//

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
template <typename T> void print_statistics(T result);
template <typename T> void print_statistics(T result, long long reference);

template <typename T> void print_statistics(T result) {
  std::string format_string =
      "{:>35} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n";
  fmt::print("{:-<135}\n", "Table of Counters ");
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
  fmt::print("{:=<135}\n", "");
}

template <typename T> void print_statistics(T result, long long reference) {
  fmt::print("{:-<135}\n", "Table of Counters ");
  fmt::print("{:>35} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n",
             "Counter Name", "#Calls", "Min [ns]", "Max [ns]", "Avg [ns]",
             "StdDev [ns]", "Percent [%]");
  for (auto& c : result) {
    if constexpr (std::is_same_v<
                      T, std::unordered_map<std::string, CounterResult>>) {
      fmt::print("{:>35} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15.2%}\n",
                 c.second.name, c.second.count, c.second.min, c.second.max,
                 c.second.mean, c.second.stddev,
                 (c.second.count * c.second.mean) / reference);
    } else {
      fmt::print("{:>35} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15.2%}\n",
                 c.name, c.count, c.min, c.max, c.mean, c.stddev,
                 (c.count * c.mean) / reference);
    }
    // fmt::print(format_string, c.name, c.count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.min))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.mean))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.max))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.stddev))).count(), std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(std::llround(c.count*c.mean))).count());
  }
  fmt::print("{:=<135}\n", "");
}

} // namespace fub

#endif /* FUB_COUNTER_COUNTER_REGISTRY_HPP */