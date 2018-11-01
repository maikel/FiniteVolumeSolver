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

#ifndef FUB_SOLVER_EULER_KINETIC_SOURCE_TERM_HPP
#define FUB_SOLVER_EULER_KINETIC_SOURCE_TERM_HPP

#include "fub/solver/SourceTermIntegrator.hpp"
#include "fub/solver/euler/IdealGas.hpp"

#include <memory>

namespace fub {
namespace euler {
/// \ingroup Euler
class KineticSourceTerm : public SourceTermIntegrator {
public:
  KineticSourceTerm(std::shared_ptr<IdealGas> ideal_gas)
      : ideal_gas_{std::move(ideal_gas)} {}

  void advanceTimeOnPatch(const std::shared_ptr<SAMRAI::hier::Patch>& patch,
                          double time_point,
                          double time_step_size) const override;

private:
  std::shared_ptr<IdealGas> ideal_gas_;
};

} // namespace euler
} // namespace fub

#endif