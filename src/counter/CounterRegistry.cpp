//
// Created by Patrick Denzler on 07.11.19.
//

#include "fub/counter/CounterRegistry.hpp"

namespace fub {

Counter* CounterRegistry::get(std::string const& name) {
  auto [iterator, _] = database_.emplace(name, name);

  return &iterator->second;
}

Timer CounterRegistry::get_timer(std::string const& name) {
  Counter* counter = this->get(name);

  return Timer(*counter);
}

std::vector<CounterResult> CounterRegistry::gather_statistics() {
  std::vector<CounterResult> statistics;
  for (auto& c : database_) {
    std::optional<CounterResult> result = c.second.gather_statistics();
    if (result.has_value()) {
      statistics.push_back(*result);
    }
  }
  return statistics;
}

std::unordered_map<std::string, CounterResult>
CounterRegistry::gather_statistics_map() {
  std::unordered_map<std::string, CounterResult> statistics;
  for (auto& c : database_) {
    std::optional<CounterResult> result = c.second.gather_statistics();
    if (result.has_value()) {
      statistics.emplace(result->name, *result);
    }
  }
  return statistics;
}

} // namespace fub