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

#ifndef FUB_OUTPUT_INVOKE_FUNCTIONS_HPP
#define FUB_OUTPUT_INVOKE_FUNCTIONS_HPP

#include "fub/ext/ProgramOptions.hpp"
#include "fub/output/BasicOutput.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

namespace fub {

template <typename Grid, typename Fn>
class AsOutput : public OutputAtFrequencyOrInterval<Grid> {
public:
  AsOutput(std::vector<std::ptrdiff_t> frequencies,
           std::vector<Duration> intervals, Fn fn)
      : OutputAtFrequencyOrInterval<Grid>(std::move(frequencies),
                                          std::move(intervals)),
        fn_(std::move(fn)) {}

  AsOutput(const std::map<std::string, pybind11::object>& vm, Fn fn)
      : OutputAtFrequencyOrInterval<Grid>(vm), fn_(std::move(fn)) {}

  void operator()(const Grid& grid) override { std::invoke(fn_, grid); }

private:
  Fn fn_;
};

template <typename Grid, typename Fn>
std::unique_ptr<AsOutput<Grid, Fn>>
MakeOutput(std::vector<std::ptrdiff_t> frequencies,
           std::vector<Duration> intervals, Fn fn) {
  return std::make_unique<AsOutput<Grid, Fn>>(
      std::move(frequencies), std::move(intervals), std::move(fn));
}

template <typename Grid>
class AnyOutput : public OutputAtFrequencyOrInterval<Grid> {
public:
  AnyOutput(std::vector<std::ptrdiff_t> frequencies,
            std::vector<Duration> intervals,
            std::function<void(const Grid&)> fn)
      : OutputAtFrequencyOrInterval<Grid>(std::move(frequencies),
                                          std::move(intervals)),
        fn_(std::move(fn)) {}

  AnyOutput(const std::map<std::string, pybind11::object>& vm,
            std::function<void(const Grid&)> fn)
      : OutputAtFrequencyOrInterval<Grid>(vm), fn_(std::move(fn)) {}

  void operator()(const Grid& grid) override { fn_(grid); }

private:
  std::function<void(const Grid&)> fn_;
};
} // namespace fub

#endif
