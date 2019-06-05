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

#ifndef FUB_AMREX_FLUX_METHOD_HPP
#define FUB_AMREX_FLUX_METHOD_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"

#include <memory>

namespace fub::amrex {
class HyperbolicSplitIntegratorContext;

struct FluxMethodBase {
  virtual ~FluxMethodBase() = default;

  virtual std::unique_ptr<FluxMethodBase> Clone() const = 0;

  virtual Duration ComputeStableDt(HyperbolicSplitIntegratorContext& context,
                                   int level, Direction dir) = 0;

  virtual void ComputeNumericFluxes(HyperbolicSplitIntegratorContext& context,
                                    int level, Duration dt, Direction dir) = 0;

  virtual int GetStencilWidth() const = 0;
};


template <typename FM> struct FluxMethodWrapper : FluxMethodBase {
  FluxMethodWrapper() = default;

  FluxMethodWrapper(const FM& fm);
  FluxMethodWrapper(FM&& fm) noexcept;

  std::unique_ptr<FluxMethodBase> Clone() const override;

  double ComputeStableDt(const ::amrex::FArrayBox& fab, const ::amrex::Box& box,
                         const ::amrex::Geometry& geom, Direction dir) override;

  void ComputeNumericFluxes(::amrex::FArrayBox& fluxes, const ::amrex::Box& box,
                            const ::amrex::FArrayBox& states,
                            const ::amrex::Geometry& geom, Duration dt,
                            Direction dir) override;

  int GetStencilWidth() const override;

  FM flux_method_{};
};

class FluxMethod {
public:
  FluxMethod() = delete;

template <typename FM, typename = std::enable_if_t<std::is_base_of<
                           FluxMethodBase, std::decay_t<FM>>::value>>
  FluxMethod(FM&& flux_method);

  FluxMethod(const FluxMethod& other);
  FluxMethod& operator=(const FluxMethod&);

  FluxMethod(FluxMethod&&) = default;
  FluxMethod& operator=(FluxMethod&&) = default;

  Duration ComputeStableDt(HyperbolicSplitIntegratorContext& context, int level,
                           Direction dir);

  void ComputeNumericFluxes(HyperbolicSplitIntegratorContext& context,
                            int level, Duration dt, Direction dir);

  int GetStencilWidth() const;

private:
  std::unique_ptr<FluxMethodBase> flux_method_;
};

template <typename FM, typename>
FluxMethod::FluxMethod(FM&& flux_method)
    : flux_method_{
          std::make_unique<std::decay_t<FM>>(std::forward<FM>(flux_method))} {}

inline FluxMethod::FluxMethod(const FluxMethod& other)
  : flux_method_{other.flux_method_->Clone()}
{}

inline FluxMethod& FluxMethod::operator=(const FluxMethod& other) {
  flux_method_ = other.flux_method_->Clone();
  return *this;
}

inline Duration
FluxMethod::ComputeStableDt(HyperbolicSplitIntegratorContext& context,
                            int level, Direction dir) {
  if (!flux_method_) {
    throw std::runtime_error("Empty flux method.");
  }
  return flux_method_->ComputeStableDt(context, level, dir);
}

inline void
FluxMethod::ComputeNumericFluxes(HyperbolicSplitIntegratorContext& context,
                                 int level, Duration dt, Direction dir)

{
  if (!flux_method_) {
    throw std::runtime_error("Empty flux method.");
  }
  return flux_method_->ComputeNumericFluxes(context, level, dt, dir);
}

inline int FluxMethod::GetStencilWidth() const
{
  if (!flux_method_) {
    throw std::runtime_error("Empty flux method.");
  }
  return flux_method_->GetStencilWidth();
}

} // namespace fub::amrex

#endif
