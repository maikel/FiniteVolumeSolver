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

#ifndef FUB_OUTPUT_COUNTER_OUTPUT_HPP
#define FUB_OUTPUT_COUNTER_OUTPUT_HPP

#include "fub/counter/CounterRegistry.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include <chrono>
#include <memory>

namespace fub {

template <typename Grid, typename PrintDuration = std::chrono::nanoseconds>
class CounterOutput : public OutputAtFrequencyOrInterval<Grid> {
public:
  CounterOutput(const ProgramOptions& po)
    : OutputAtFrequencyOrInterval<Grid>(po), reference_{std::chrono::steady_clock::now()} {}

  CounterOutput(std::chrono::steady_clock::time_point reference,
                std::vector<std::ptrdiff_t> frequencies,
                std::vector<Duration> intervals)
      : OutputAtFrequencyOrInterval<Grid>(std::move(frequencies),
                                          std::move(intervals)),
        reference_(reference) {}

  void operator()(const Grid& grid) override {
    std::chrono::steady_clock::time_point now =
        std::chrono::steady_clock::now();
    auto diff =
        std::chrono::duration_cast<std::chrono::nanoseconds>(now - reference_);
    std::vector<CounterResult> statistics = grid.GetPatchHierarchy().GetCounterRegistry()->gather_statistics();
    if (statistics.size()) {
      print_statistics<PrintDuration>(statistics, diff.count());
    }
  }

private:
  std::chrono::steady_clock::time_point reference_;
};

} // namespace fub

#endif