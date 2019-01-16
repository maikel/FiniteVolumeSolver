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

#ifndef FUB_EULER_PERFECT_GAS_EQUATION_HPP
#define FUB_EULER_PERFECT_GAS_EQUATION_HPP

#include "fub/SAMRAI/ideal_gas/IdealGasEquation.hpp"

namespace fub {
namespace ideal_gas {

/// \ingroup Euler
class PerfectGasEquation : public IdealGasEquation {
public:
  PerfectGasEquation(std::string name, SAMRAI::tbox::Dimension dim)
      : IdealGasEquation(std::move(name), dim, 1) {}

  /// \brief Computes a complete state from a given conservative state.
  void FillFromCons(const CompletePatchData& q,
                    const ConsPatchData& u) const override;

  /// \brief Computes a complete state from a given conservative state.
  void FillFromPrim(const CompletePatchData& q,
                    const PrimPatchData& w) const override;
};

} // namespace ideal_gas
} // namespace fub

#endif