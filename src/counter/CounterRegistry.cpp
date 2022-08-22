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

#include "fub/counter/CounterRegistry.hpp"

namespace fub {

Counter* CounterRegistry::get(std::string const& name) {
  auto iterator = database_.emplace(name, name).first;
  return &iterator->second;
}

Timer CounterRegistry::get_timer(std::string const& name) {
  Counter* counter = this->get(name);
  return Timer(*counter);
}

std::vector<CounterResult> CounterRegistry::gather_statistics() {
  std::vector<CounterResult> statistics;
  for (const std::pair<const std::string, Counter>& c : database_) {
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
  for (const std::pair<const std::string, Counter>& c : database_) {
    std::optional<CounterResult> result = c.second.gather_statistics();
    if (result.has_value()) {
      statistics.emplace(result->name, *result);
    }
  }
  return statistics;
}

} // namespace fub