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

#ifndef FUB_OTPUT_OUTPUT_AT_FREQ_AND_INTERV_HPP
#define FUB_OTPUT_OUTPUT_AT_FREQ_AND_INTERV_HPP

#include "fub/ext/ProgramOptions.hpp"
#include "fub/output/BasicOutput.hpp"

#include <pybind11/stl.h>

namespace fub {

template <typename GriddingAlgorithm>
class OutputAtFrequencyOrInterval : public BasicOutput<GriddingAlgorithm> {
public:
  OutputAtFrequencyOrInterval() = default;

  /// Initialize this module by a program options map.
  OutputAtFrequencyOrInterval(
      const std::map<std::string, pybind11::object>& options) {
    frequencies_ = GetOptionOr(options, "frequencies", frequencies_);
    std::vector<double> intervals =
        GetOptionOr(options, "intervals", std::vector<double>{});
    std::transform(intervals.begin(), intervals.end(),
                   std::back_inserter(intervals_),
                   [](double dur) { return Duration(dur); });
    smallest_time_step_size_ = GetOptionOr(options, "smallest_time_step_size",
                                           smallest_time_step_size_);
  }

  OutputAtFrequencyOrInterval(
      std::vector<std::ptrdiff_t> frequencies,
      std::vector<fub::Duration> intervals,
      Duration smallest_time_step_size = Duration(1e-12))
      : frequencies_(frequencies), intervals_(intervals),
        smallest_time_step_size_(smallest_time_step_size) {}

  Duration NextOutputTime(Duration time_point) override {
    Duration next_output_time(std::numeric_limits<double>::max());
    for (std::size_t i = 0; i < intervals_.size(); ++i) {
      if (intervals_[i].count() > 0.0) {
        Duration dt(intervals_[i].count() -
                    std::fmod(time_point.count(), intervals_[i].count()));
        dt += smallest_time_step_size_;
        next_output_time = std::min<Duration>(next_output_time, time_point + dt);
      }
    }
    return next_output_time;
  }

  bool ShallOutputNow(const GriddingAlgorithm& grid) override {
    const std::ptrdiff_t cycle = grid.GetCycles();
    if (std::any_of(
            frequencies_.begin(), frequencies_.end(),
            [cycle](std::ptrdiff_t freq) { return cycle % freq == 0; })) {
      return true;
    }
    const Duration time_point = grid.GetTimePoint();
    if (time_point == Duration{} ||
        std::any_of(intervals_.begin(), intervals_.end(),
                    [this, time_point](Duration freq) {
                      return std::abs(
                                 std::fmod(time_point.count(), freq.count())) <
                             2 * smallest_time_step_size_.count();
                    })) {
      return true;
    }
    return false;
  }

  template <typename Log>
  void Print(Log& log) {
    BOOST_LOG(log) << fmt::format(" - frequencies = {{{}}} [-]", fmt::join(frequencies_, ", "));
    std::vector<double> interval_counts(intervals_.size());
    std::transform(intervals_.begin(), intervals_.end(), interval_counts.begin(), [](fub::Duration dur) {
      return dur.count();
    });
    BOOST_LOG(log) << fmt::format(" - intervals = {{{}}} [s]", fmt::join(interval_counts, ", "));
    BOOST_LOG(log) << fmt::format(" - smallest_time_step_size = {} [s]", smallest_time_step_size_.count());
  }

protected:
  std::vector<std::ptrdiff_t> frequencies_{};
  std::vector<Duration> intervals_{};
  Duration smallest_time_step_size_{1e-12};
};

} // namespace fub

#endif
