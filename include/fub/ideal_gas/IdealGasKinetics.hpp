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

/// \defgroup Euler Euler Equations
/// This module defines classes and methods to be used to concrete solver
/// implementations wrt the euler equations.

#ifndef FUB_IDEAL_GAS_IDEAL_GAS_KINETICS_HPP
#define FUB_IDEAL_GAS_IDEAL_GAS_KINETICS_HPP

#include "fub/ideal_gas/IdealGasEquation.hpp"

#include <array>

namespace fub {
namespace ideal_gas {

/// \ingroup IdealGas
class IdealGasKinetics : public IdealGasEquation {
public:
  using IdealGasEquation::IdealGasEquation;

  /// \brief Computes a complete state from a given conservative state.
  virtual void AdvanceSourceTerm(const CompleteStatePatchData& complete,
                                 double time_step_size) const = 0;
};

} // namespace ideal_gas
} // namespace fub

#endif