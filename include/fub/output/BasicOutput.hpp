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

#ifndef FUB_OUTPUT_BASIC_OUTPUT_HPP
#define FUB_OUTPUT_BASIC_OUTPUT_HPP

#include "fub/Duration.hpp"

#include <memory>

namespace fub {

/// This is a abstract base class for an output strategy.
/// Objects of this class are intended to be passt to the fub::RunSimulation function.
template <typename GriddingAlgorithm> struct BasicOutput {
  /// The destructor needs to be virtual to prevent leaking resources.
  virtual ~BasicOutput() = default;

  /// Returns the time point at which the simulation shall stop to do some output.
  virtual Duration NextOutputTime(Duration time_point) = 0;

  /// Returns true if this output class shall be invoked at the specified time point.
  virtual bool ShallOutputNow(const GriddingAlgorithm& grid) = 0;

  /// Invoke the actual output logic.
  virtual void
  operator()(const GriddingAlgorithm& grid) = 0;
};

} // namespace fub

#endif