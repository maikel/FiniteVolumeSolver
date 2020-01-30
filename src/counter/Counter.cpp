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

#include "fub/counter/Counter.hpp"

#include <cmath>
#include <mpi.h>

namespace fub {

Counter::Counter(std::string name) : name_(std::move(name)) {}

void Counter::add(long long value) { values_.push_back(value); }

std::string const& Counter::get_name() const noexcept { return name_; }

long long Counter::get_value() const noexcept {
  return std::accumulate(values_.begin(), values_.end(), 0ll);
}

unsigned long Counter::get_count() const noexcept { return values_.size(); }

long long Counter::get_min() const noexcept {
  if (values_.empty()) {
    return 0;
  }
  return *std::min_element(values_.begin(), values_.end());
}

long long Counter::get_max() const noexcept {
  if (values_.empty()) {
    return 0;
  }
  return *std::max_element(values_.begin(), values_.end());
}

double Counter::get_mean() const noexcept { return calculate_mean(values_); }

double Counter::get_variance() const noexcept {
  return calculate_variance(values_);
}

double Counter::get_stddev() const noexcept {
  return std::sqrt(get_variance());
}

void Counter::reset() { values_.clear(); }

std::optional<CounterResult> Counter::gather_statistics(MPI_Comm comm,
                                                        int root) const {
  int mpi_size = -1;
  int mpi_rank = -1;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  long long local_min = get_min();
  long long local_max = get_max();
  double local_mean = get_mean();

  const int send_size = 1;
  const int recv_count_per_process = 1;

  if (mpi_rank == root) {
    auto received_mins = std::vector<long long>(mpi_size);
    auto received_maxs = std::vector<long long>(mpi_size);
    auto received_means = std::vector<double>(mpi_size);

    MPI_Gather(&local_min, send_size, MPI_LONG_LONG, received_mins.data(),
               recv_count_per_process, MPI_LONG_LONG, root, comm);
    MPI_Gather(&local_max, send_size, MPI_LONG_LONG, received_maxs.data(),
               recv_count_per_process, MPI_LONG_LONG, root, comm);
    MPI_Gather(&local_mean, send_size, MPI_DOUBLE, received_means.data(),
               recv_count_per_process, MPI_DOUBLE, root, comm);

    double global_variance = calculate_variance(received_means);

    return CounterResult{
        get_name(),
        get_count(),
        *std::min_element(received_mins.begin(), received_mins.end()),
        *std::max_element(received_maxs.begin(), received_maxs.end()),
        calculate_mean(received_means),
        std::sqrt(global_variance),
        global_variance};

  } else {
    MPI_Gather(&local_min, send_size, MPI_LONG_LONG, nullptr,
               recv_count_per_process, MPI_LONG_LONG, root, comm);
    MPI_Gather(&local_max, send_size, MPI_LONG_LONG, nullptr,
               recv_count_per_process, MPI_LONG_LONG, root, comm);
    MPI_Gather(&local_mean, send_size, MPI_DOUBLE, nullptr,
               recv_count_per_process, MPI_DOUBLE, root, comm);

    return {};
  }
}

} // namespace fub