// Copyright (c) 2018 Maikel Nadolski
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

#ifndef FUB_SPLITTING_METHOD_HPP
#define FUB_SPLITTING_METHOD_HPP

#include "fub/Duration.hpp"
#include "fub/core/function_ref.hpp"

namespace fub {
/// \ingroup Abstract
/// This is an abstract interface for implementing splitting methods.
//
/// Based on a strategy to evolve two operators this base class uses that
/// derived strategy recursively to generalize the derived splitting to N
/// operators.
struct SplittingMethod {
  using AdvanceFunction = function_ref<void(Duration)>;

  virtual ~SplittingMethod() = default;

  template <typename F> void Advance(Duration dt, F advance) const {
    advance(dt);
  }

  /// This method recurisvely applies the base case.
  template <typename F, typename... Fs>
  void Advance(Duration dt, AdvanceFunction advance1, AdvanceFunction advance2,
               F advance3, Fs... advances) const {
    Advance(dt, advance1, [&](Duration time_step_size) {
      Advance(time_step_size, advance2, advance3, advances...);
    });
  }

  /// This is the base case of applying the splitting method with two operators.
  virtual void Advance(Duration time_step_size, AdvanceFunction advance1,
                       AdvanceFunction advance2) const = 0;
};

} // namespace fub

#endif
