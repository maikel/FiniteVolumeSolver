//
// Created by Patrick Denzler on 07.11.19.
//

#include "fub/counter/Counter.hpp"

#include <cmath>
#include <mpi.h>

namespace fub {

Counter::Counter(std::string name) : name_(std::move(name)) {}

void Counter::add(long long value) { values_.push_back(value); }

std::string const& Counter::get_name() { return name_; }

long long Counter::get_value() {
  return std::accumulate(values_.begin(), values_.end(), 0ll);
}

unsigned long Counter::get_count() { return values_.size(); }

long long Counter::get_min() {
  return *std::min_element(values_.begin(), values_.end());
}

long long Counter::get_max() {
  return *std::max_element(values_.begin(), values_.end());
}

double Counter::get_mean() { return calculate_mean(values_); }

double Counter::get_variance() { return calculate_variance(values_); }

double Counter::get_stddev() { return std::sqrt(get_variance()); }

void Counter::reset() { values_.clear(); }

std::optional<CounterResult> Counter::gather_statistics(int root) {
  int mpi_size, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  auto local_min = get_min();
  auto local_max = get_max();
  auto local_mean = get_mean();

  if (mpi_rank == root) {
    auto received_mins = std::vector<long long>(mpi_size);
    auto received_maxs = std::vector<long long>(mpi_size);
    auto received_means = std::vector<double>(mpi_size);

    MPI_Gather(&local_min, 1, MPI_LONG_LONG, received_mins.data(), 1,
               MPI_LONG_LONG, root, MPI_COMM_WORLD);
    MPI_Gather(&local_max, 1, MPI_LONG_LONG, received_maxs.data(), 1,
               MPI_LONG_LONG, root, MPI_COMM_WORLD);
    MPI_Gather(&local_mean, 1, MPI_DOUBLE, received_means.data(), 1, MPI_DOUBLE,
               root, MPI_COMM_WORLD);

    auto global_variance = calculate_variance(received_means);

    return CounterResult{
        .name = get_name(),
        .count = get_count(),
        .min = *std::min_element(received_mins.begin(), received_mins.end()),
        .max = *std::max_element(received_maxs.begin(), received_maxs.end()),
        .mean = calculate_mean(received_means),
        .variance = global_variance,
        .stddev = std::sqrt(global_variance)};

  } else {
    MPI_Gather(&local_min, 1, MPI_LONG_LONG, nullptr, 1, MPI_LONG_LONG, root,
               MPI_COMM_WORLD);
    MPI_Gather(&local_max, 1, MPI_LONG_LONG, nullptr, 1, MPI_LONG_LONG, root,
               MPI_COMM_WORLD);
    MPI_Gather(&local_mean, 1, MPI_DOUBLE, nullptr, 1, MPI_DOUBLE, root,
               MPI_COMM_WORLD);

    return {};
  }
}

} // namespace fub