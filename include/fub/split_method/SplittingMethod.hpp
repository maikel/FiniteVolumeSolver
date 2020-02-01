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

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/core/function_ref.hpp"

#include <boost/outcome.hpp>

namespace fub {
/// \ingroup Abstract
/// This is an abstract interface for implementing splitting methods.
//
/// Based on a strategy to evolve two operators this base class uses that
/// derived strategy recursively to generalize the derived splitting to N
/// operators.
struct SplittingMethod {
  using AdvanceFunction =
      function_ref<boost::outcome_v2::result<void, TimeStepTooLarge>(Duration)>;

  virtual ~SplittingMethod() = default;

  template <typename F>
  boost::outcome_v2::result<void, TimeStepTooLarge> Advance(Duration dt,
                                                            F advance) const {
    return advance(dt);
  }

  /// This method recurisvely applies the base case.
  template <typename F, typename... Fs>
  boost::outcome_v2::result<void, TimeStepTooLarge>
  Advance(Duration dt, AdvanceFunction advance1, AdvanceFunction advance2,
          F advance3, Fs... advances) const {
    return Advance(dt, advance1, [&](Duration time_step_size) {
      return Advance(time_step_size, advance2, advance3, advances...);
    });
  }

  /// This is the base case of applying the splitting method with two operators.
  virtual boost::outcome_v2::result<void, TimeStepTooLarge>
  Advance(Duration time_step_size, AdvanceFunction advance1,
          AdvanceFunction advance2) const = 0;
};

template <int Rank>
constexpr std::array<Direction, static_cast<std::size_t>(Rank)>
MakeSplitDirections(std::pair<int, int> subcycle) noexcept {
  std::array<Direction, static_cast<std::size_t>(Rank)> directions{};
  int is_odd = subcycle.first % 2;
  int is_even = !is_odd;
  for (int i = 0; i < Rank; ++i) {
    int dir = (is_even * i + is_odd * (Rank - 1 - i));
    directions[static_cast<std::size_t>(i)] = static_cast<Direction>(dir);
  }
  return directions;
}

} // namespace fub

#endif
