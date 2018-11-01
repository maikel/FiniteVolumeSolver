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

#ifndef FUB_SOLVER_INITIALIZE_STRATEGY_HPP
#define FUB_SOLVER_INITIALIZE_STRATEGY_HPP

#include "SAMRAI/hier/Patch.h"
#include <memory>

namespace fub {
/// \ingroup Solver
/// @{
/// \defgroup Abstract Abstract Interfaces
/// This submodule lists all abstract interfaces of the Solver module.
/// @}

/// \ingroup Abstract
/// \brief This is an abstract interface for setting initial data on a single
/// patch.
struct InitialCondition {
  virtual ~InitialCondition() = default;

  /// Initialize patch data values on the specified patch.
  ///
  /// A concrete implementation has to know which patch data ids have to be
  /// filled. This depends on problem basis.
  virtual void initializeDataOnPatch(const SAMRAI::hier::Patch&) const = 0;
};

} // namespace fub

#endif