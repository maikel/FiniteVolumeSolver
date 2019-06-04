// Copyright (c) 2019 Maikel Nadolski
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

#ifndef FUB_AMREX_CUTCELL_FLUX_METHOD_HPP
#define FUB_AMREX_CUTCELL_FLUX_METHOD_HPP

#include <memory>
#include "fub/AMReX/cutcell/IntegratorContext.hpp"

namespace fub {
namespace amrex {
namespace cutcell {
namespace detail {
struct FluxMethodBase {
  virtual ~FluxMethodBase() = default;
  virtual std::unique_ptr<FluxMethodBase> Clone() const = 0;
  virtual void PreAdvanceHierarchy(IntegratorContext& context) = 0;
  virtual void ComputeNumericFluxes(IntegratorContext& context, int level, Duration dt, Direction dir) = 0;
};

template <typename T, typename... Args>
using PreAdvanceHierarchyT = decltype(std::declval<T>().PreAdvanceHierarchy(std::declval<Args>()...));

template <typename T, typename... Args>
using ComputeNumericFluxesT = decltype(std::declval<T>().ComputeNumericFluxes(std::declval<Args>()...));

template <typename FM>
struct FluxMethodWrapper : public FluxMethodBase {
  using Equation = std::decay_t<decltype(std::declval<FM&>().GetEquation())>;
  using Conservative = ::fub::Conservative<Equation>;
  using Complete = ::fub::Complete<Equation>;

  static constexpr int Rank = Equation::Rank();

  static constexpr bool HasPreAdvanceHierarchy = is_detected<PreAdvanceHierarchyT, FM&, IntegratorContext&>();
  static constexpr bool HasComputeNumericFluxesOnContext = is_detected<ComputeNumericFluxesT, FM&, IntegratorContext&, int, Duration, Direction>();


  FluxMethodWrapper(const FM& fm) : method_{fm} {}
  FluxMethodWrapper(FM&& fm) : method_{std::move(fm)} {}

  std::unique_ptr<FluxMethodBase> Clone() const override;
  void PreAdvanceHierarchy(IntegratorContext& context) override;
  void ComputeNumericFluxes(IntegratorContext& context, int level, Duration dt, Direction dir) override;

  FM method_;
};
}

class FluxMethod {
public:

private:
  std::unique_ptr<detail::FluxMethodBase> base_;
};

namespace detail {
template <typename FM>
std::unique_ptr<FluxMethodBase> FluxMethodWrapper<FM>::Clone() const override {
  return std::make_unique<FluxMethodWrapper<FM>>(method_);
}

template <typename FM>
void FluxMethodWrapper<FM>::PreAdvanceHierarchy(IntegratorContext& context) override {
  if constexpr (HasPreAdvanceHierarchy) {
    method_.PreAdvaceHierarchy(context);
  }
}

template <typename FM>
void FluxMethodWrapper<FM>::ComputeNumericFluxes(IntegratorContext& context, int level, Duration dt, Direction dir) override {
  if constexpr (HasComputeNumericFluxesOnContext) {
    method_.ComputeNumericFluxes(context, level, dt, dir);
  }
}
}

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
