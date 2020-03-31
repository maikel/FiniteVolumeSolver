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

#ifndef FUB_STRANG_SPLITTING_HPP
#define FUB_STRANG_SPLITTING_HPP

#include "fub/Duration.hpp"
#include "fub/core/function_ref.hpp"
#include "fub/split_method/SplittingMethod.hpp"

namespace fub {

/// \ingroup Solver
/// \brief This class implements a second order splitting method.
struct StrangSplitting : public SplittingMethod {
public:
  using AdvanceFunction = SplittingMethod::AdvanceFunction;

  using SplittingMethod::Advance;

  /// \brief Invokes two operators a1, a2 by the following scheme:
  /// a1(dt/2); a2(dt/2); a2(dt/2); a1(dt/2);
  ///
  /// \param[in] time_step_size  The time by which the hierarchy shall be
  ///                            advanced.
  /// \param[in] operator1  The first operator in the splitting.
  /// \param[in] operator2  The second operator in the spligging.
  boost::outcome_v2::result<void, TimeStepTooLarge>
  Advance(Duration time_step_size, AdvanceFunction operator1,
          AdvanceFunction operator2) const override;
  
  /// This is the base case of applying the splitting method with two operators.
  boost::outcome_v2::result<void, TimeStepTooLarge>
  Advance(Duration time_step_size, AdvanceFunction advance1,
          AdvanceFunction advance2, AdvanceFunction advance3) const override;
};

} // namespace fub

#endif
