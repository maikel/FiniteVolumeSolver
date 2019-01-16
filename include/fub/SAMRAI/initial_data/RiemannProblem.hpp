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

#ifndef FUB_INITIAL_DATA_RIEMANN_PROBLEM_HPP
#define FUB_INITIAL_DATA_RIEMANN_PROBLEM_HPP

#include "fub/core/type_traits.hpp"
#include "fub/geometry/PolymorphicGeometry.hpp"
#include "fub/SAMRAI/InitialCondition.hpp"

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Patch.h"

#include <memory>

namespace fub {

class RiemannProblem : public fub::InitialCondition {
public:
  RiemannProblem(PolymorphicGeometry geometry)
      : geometry_{std::move(geometry)} {}

  void InitializeDataOnPatch(const SAMRAI::hier::Patch& patch) const override;

  virtual void PostInitialize(const SAMRAI::hier::Patch& patch) const {}

  const PolymorphicGeometry& GetGeometry() const noexcept { return geometry_; }

private:
  virtual void FillLeftState(const SAMRAI::hier::Patch& patch,
                             const SAMRAI::hier::Index& index) const = 0;

  virtual void FillRightState(const SAMRAI::hier::Patch& patch,
                              const SAMRAI::hier::Index& index) const = 0;

  PolymorphicGeometry geometry_;
};

} // namespace fub

#endif