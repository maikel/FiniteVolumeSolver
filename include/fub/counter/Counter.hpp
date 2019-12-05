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

#ifndef FUB_COUNTER_COUNTER_HPP
#define FUB_COUNTER_COUNTER_HPP

#include "CounterResult.hpp"

#include <algorithm>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include <mpi.h>

namespace fub {

class Counter {
public:
  explicit Counter(std::string name);
  void add(long long value);
  std::string const& get_name() const noexcept;
  long long get_value() const noexcept;
  unsigned long get_count() const noexcept;
  long long get_min() const noexcept;
  long long get_max() const noexcept;
  double get_mean() const noexcept;
  double get_variance() const noexcept;
  double get_stddev() const noexcept;
  void reset();

  std::optional<CounterResult> gather_statistics(MPI_Comm comm = MPI_COMM_WORLD,
                                                 int root = 0) const;

private:
  std::string name_;
  std::vector<long long> values_;

  template <typename T> double calculate_mean(const std::vector<T>& v) const noexcept;
  template <typename T> double calculate_variance(const std::vector<T>& v) const noexcept;
};

template <typename T> double Counter::calculate_mean(const std::vector<T>& v) const noexcept {
  if (v.size()) {
    auto sum = std::accumulate(v.begin(), v.end(), static_cast<T>(0));
    return static_cast<double>(sum) / v.size();
  }
  return 0;
}

template <typename T>
double Counter::calculate_variance(const std::vector<T>& v) const noexcept {
  if (v.size()) {
    const double mean = calculate_mean(v);
    const double sq_sum = std::inner_product(
        v.begin(), v.end(), v.begin(), 0.0, std::plus<>(),
        [mean](double x, double y) { return (x - mean) * (y - mean); });
    return (sq_sum / v.size());
  }
  return 0;
}

} // namespace fub

#endif /* FUB_COUNTER_COUNTER_HPP */