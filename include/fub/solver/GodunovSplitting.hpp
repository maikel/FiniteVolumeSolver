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

#ifndef FUB_SOLVER_GODUNOV_SPLITTING_HPP
#define FUB_SOLVER_GODUNOV_SPLITTING_HPP

#include "fub/core/function_ref.hpp"
#include "fub/solver/SplittingMethod.hpp"

namespace fub {

/// \ingroup Solver
/// \brief This class implements a first order splitting method.
struct GodunovSplitting : public SplittingMethod {
public:
  using AdvanceFunction = SplittingMethod::AdvanceFunction;

  /// \brief Invokes two operators in sequence.
  ///
  /// \param[in, out] hierarchy  This is the hierarchy where operators are 
  ///                            acting on.
  /// \param[in] time_point  The current time point of the simulation.
  /// \param[in] time_step_size  The time by which the hierarchy shall be 
  ///                            advanced.
  /// \param[in] operator1  The first operator in the splitting.
  /// \param[in] operator2  The second operator in the spligging.
  void
  advanceTime(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
              double time_point, double time_step_size,
              AdvanceFunction operator1,
              AdvanceFunction operator2) const override;
};

} // namespace fub

#endif