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

#ifndef FUB_OUTPUT_MULTIPLE_OUTPUTS_HPP
#define FUB_OUTPUT_MULTIPLE_OUTPUTS_HPP

#include "fub/ext/ProgramOptions.hpp"
#include "fub/output/OutputFactory.hpp"

#include <pybind11/stl.h>

#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>

namespace fub {

/// \ingroup Output
template <typename Grid> class MultipleOutputs : public BasicOutput<Grid> {
public:
  using ProgramOptions = std::map<std::string, pybind11::object>;

  MultipleOutputs() = default;

  MultipleOutputs(OutputFactory<Grid> factory, const ProgramOptions& opts)
      : factory_(std::move(factory)), outputs_{} {
    std::vector<int> frequencies =
        GetOptionOr(opts, "frequencies", std::vector<int>{});
    std::vector<double> intervals =
        GetOptionOr(opts, "intervals", std::vector<double>{});
    std::vector<pybind11::dict> outputs{};
    outputs = GetOptionOr(opts, "outputs", outputs);
    for (const pybind11::dict& output : outputs) {
      auto opts = ToMap(output);
      if (auto type = opts.find("type"); type != opts.end()) {
        using namespace std::literals;
        const std::string type_name = type->second.cast<std::string>();
        opts.emplace(std::pair{"frequencies", pybind11::cast(frequencies)});
        opts.emplace(std::pair{"intervals", pybind11::cast(intervals)});
        outputs_.push_back(factory_.MakeOutput(type_name, opts));
      }
    }
  }

  template <typename T> void AddOutput(std::unique_ptr<T>&& output) {
    outputs_.emplace_back(output.release());
  }

  /// Returns the time point at which the simulation shall stop to do some
  /// output.
  Duration NextOutputTime(Duration time_point) override {
    return std::accumulate(
        outputs_.begin(), outputs_.end(),
        Duration(std::numeric_limits<double>::max()),
        [time_point](Duration min,
                     const std::unique_ptr<BasicOutput<Grid>>& output) {
          return output ? std::min(min, output->NextOutputTime(time_point))
                        : min;
        });
  }

  /// Returns true if this output class shall be invoked at the specified time
  /// point.
  bool ShallOutputNow(const Grid& grid) override {
    return std::any_of(outputs_.begin(), outputs_.end(),
                       [&grid](const auto& output) {
                         return output && output->ShallOutputNow(grid);
                       });
  }

  /// Invoke the actual output logic.
  void operator()(const Grid& grid) override {
    for (const auto& output : outputs_) {
      if (output && output->ShallOutputNow(grid)) {
        (*output)(grid);
      }
    }
  }

private:
  OutputFactory<Grid> factory_{};
  std::vector<std::unique_ptr<BasicOutput<Grid>>> outputs_{};
};

} // namespace fub

#endif