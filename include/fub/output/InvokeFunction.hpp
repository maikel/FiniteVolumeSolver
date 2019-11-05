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

template <typename Grid>
class InvokeFunction : public OutputAtFrequencyOrInterval<Grid> {
public:
  InvokeFunction(std::vector<std::ptrdiff_t> frequencies,
      std::vector<fub::Duration> intervals, std::function<void(const Grid&)> fn)
    : OutputAtFrequencyOrInterval<Grid>(std::move(frequencies), std::move(intervals)), fn_(std::move(fn))
  {}

  InvokeFunction(const std::map<std::string, pybind11::object>& vm, std::function<void(const Grid&)> fn)
    : OutputAtFrequencyOrInterval<Grid>(vm), fn_(std::move(fn))
  {}

  void
  operator()(const Grid& grid) override {
    if (fn_) {
      fn_(grid);
    }
  }

private:
  std::function<void(const Grid&)> fn_;
};

}

#endif
