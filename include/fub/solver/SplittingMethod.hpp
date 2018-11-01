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

#ifndef FUB_SOLVER_SPLITTING_METHOD_HPP
#define FUB_SOLVER_SPLITTING_METHOD_HPP

#include "fub/core/function_ref.hpp"

#include "SAMRAI/hier/PatchHierarchy.h"

#include <memory>

namespace fub {
/// \ingroup Abstract
/// This is an abstract interface for implementing splitting methods.
///
/// A splitting method takes two operators and executes them in some order such
/// that they solve a combined problem.
struct SplittingMethod {
public:
  using AdvanceFunction = function_ref<void(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&, double, double)>;

  virtual ~SplittingMethod() = default;

  /// This method recurisvely applies the base case.
  template <typename F, typename... Fs>
  void
  advanceTime(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
              double time_point, double time_step_size,
              AdvanceFunction advance1, AdvanceFunction advance2, F advance3,
              Fs... advances) const {
    advanceTime(
        hierarchy, time_point, time_step_size, advance1,
        [&](const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            double time_point, double time_step_size) {
          advanceTime(hierarchy, time_point, time_step_size, advance2, advance3,
                      advances...);
        });
  }

  /// This is the base case of applying the splitting method with two operators.
  virtual void
  advanceTime(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
              double time_point, double time_step_size,
              AdvanceFunction advance1, AdvanceFunction advance2) const = 0;
};

} // namespace fub

#endif