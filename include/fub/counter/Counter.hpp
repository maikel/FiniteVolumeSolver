//
// Created by Patrick Denzler on 07.11.19.
//

#ifndef FUB_COUNTER_COUNTER_HPP
#define FUB_COUNTER_COUNTER_HPP

#include "CounterResult.hpp"

#include <string>
#include <vector>
#include <optional>
#include <numeric>

namespace fub {

class Counter {
public:
  explicit Counter(std::string name);
  void add(long long value);
  std::string const& get_name();
  long long get_value();
  unsigned long get_count();
  long long get_min();
  long long get_max();
  double get_mean();
  double get_variance();
  double get_stddev();
  void reset();

  std::optional<CounterResult> gather_statistics(int root = 0);

private:
  std::string name_;
  std::vector<long long> values_;
  double multiplier_{};

  template <typename T> double calculate_mean(std::vector<T> v);
  template <typename T> double calculate_variance(std::vector<T> v);
};

template <typename T> double Counter::calculate_mean(std::vector<T> v) {
  auto sum = std::accumulate(v.begin(), v.end(), static_cast<T>(0));

  return static_cast<double>(sum) / v.size();
}

template <typename T> double Counter::calculate_variance(std::vector<T> v) {
  std::vector<double> diff(v.size());
  double mean = calculate_mean(v);

  std::transform(v.begin(), v.end(), diff.begin(),
                 [mean](double x) { return x - mean; });

  double sq_sum =
      std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);

  return (sq_sum / v.size());
}

} // namespace fub

#endif /* FUB_COUNTER_COUNTER_HPP */