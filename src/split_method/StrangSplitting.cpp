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

#include "fub/split_method/StrangSplitting.hpp"

namespace fub {

boost::outcome_v2::result<void, TimeStepTooLarge>
StrangSplitting::Advance(Duration time_step_size, AdvanceFunction advance1,
                         AdvanceFunction advance2) const {
  boost::outcome_v2::result<void, TimeStepTooLarge> result =
      advance1(0.5 * time_step_size);
  if (!result) {
    return result.as_failure();
  }
  result = advance2(0.5 * time_step_size);
  if (!result) {
    return result.as_failure();
  }
  result = advance2(0.5 * time_step_size);
  if (!result) {
    return result.as_failure();
  }
  return advance1(0.5 * time_step_size);
}

boost::outcome_v2::result<void, TimeStepTooLarge>
StrangSplitting::Advance(Duration time_step_size, AdvanceFunction advance1,
                         AdvanceFunction advance2, AdvanceFunction advance3) const {
  boost::outcome_v2::result<void, TimeStepTooLarge> result =
      advance1(0.5 * time_step_size);
  if (!result) {
    return result.as_failure();
  }
  result = advance2(0.5 * time_step_size);
  if (!result) {
    return result.as_failure();
  }
  result = advance3(0.5 * time_step_size);
  if (!result) {
    return result.as_failure();
  }
  result = advance3(0.5 * time_step_size);
  if (!result) {
    return result.as_failure();
  }
  result = advance2(0.5 * time_step_size);
  if (!result) {
    return result.as_failure();
  }
  return advance1(0.5 * time_step_size);
}

} // namespace fub
