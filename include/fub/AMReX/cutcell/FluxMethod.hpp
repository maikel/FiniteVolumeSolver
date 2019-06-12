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

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"

#include <memory>

namespace fub::amrex::cutcell {
class IntegratorContext;

struct FluxMethodBase {
  virtual ~FluxMethodBase() = default;

  virtual std::unique_ptr<FluxMethodBase> Clone() const = 0;

  virtual void PreAdvanceHierarchy(IntegratorContext& context) = 0;

  virtual void ComputeNumericFluxes(IntegratorContext& context, int level,
                                    Duration dt, Direction dir) = 0;

  virtual Duration ComputeStableDt(IntegratorContext& context, int level,
                                   Direction dir) = 0;

  virtual int GetStencilWidth() const = 0;
};

class FluxMethod {
public:
  template <typename FM> FluxMethod(FM&& fm);

  void PreAdvanceHierarchy(IntegratorContext& context);

  void ComputeNumericFluxes(IntegratorContext& context, int level, Duration dt,
                            Direction dir);

  Duration ComputeStableDt(IntegratorContext& context, int level,
                           Direction dir);

  int GetStencilWidth() const;

private:
  std::unique_ptr<FluxMethodBase> flux_method_;
};

template <typename FM>
FluxMethod::FluxMethod(FM&& flux_method)
    : flux_method_(
          std::make_unique<std::decay_t<FM>>(std::forward<FM>(flux_method_))) {}

inline void FluxMethod::PreAdvanceHierarchy(IntegratorContext& context) {
  if (!flux_method_) {
    throw std::runtime_error("Empty flux method.");
  }
  flux_method_->PreAdvanceHierarchy(context);
}

inline void FluxMethod::ComputeNumericFluxes(IntegratorContext& context,
                                             int level, Duration dt,
                                             Direction dir) {
  if (!flux_method_) {
    throw std::runtime_error("Empty flux method.");
  }
  flux_method_->ComputeNumericFluxes(context, level, dt, dir);
}

inline Duration FluxMethod::ComputeStableDt(IntegratorContext& context,
                                            int level, Direction dir) {
  if (!flux_method_) {
    throw std::runtime_error("Empty flux method.");
  }
  return flux_method_->ComputeStableDt(context, level, dir);
}

inline int FluxMethod::GetStencilWidth() const {
  if (!flux_method_) {
    throw std::runtime_error("Empty flux method.");
  }
  return flux_method_->GetStencilWidth();
}

} // namespace fub::amrex::cutcell

#endif
