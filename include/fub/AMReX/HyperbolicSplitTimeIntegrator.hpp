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

#ifndef FUB_AMREX_HYPERBOLIC_SPLIT_TIME_INTEGRATOR_HPP
#define FUB_AMREX_HYPERBOLIC_SPLIT_TIME_INTEGRATOR_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/State.hpp"

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <memory>

namespace fub::amrex {

struct HyperbolicSplitTimeIntegratorBase {
  virtual ~HyperbolicSplitTimeIntegratorBase() = default;
  virtual void UpdateConservatively(::amrex::MultiFab& dest,
                                    const ::amrex::MultiFab& src,
                                    const ::amrex::MultiFab& fluxes,
                                    const ::amrex::Geometry& geom, Duration dt,
                                    Direction dir) = 0;
  virtual std::unique_ptr<HyperbolicSplitTimeIntegratorBase> Clone() const = 0;
};

template <typename T>
struct HyperbolicSplitTimeIntegratorWrapper
    : HyperbolicSplitTimeIntegratorBase {
  HyperbolicSplitTimeIntegratorWrapper(const T& x) : time_integrator_{x} {}
  HyperbolicSplitTimeIntegratorWrapper(T&& x)
      : time_integrator_{std::move(x)} {}

  std::unique_ptr<HyperbolicSplitTimeIntegratorBase> Clone() const override {
    return std::make_unique<HyperbolicSplitTimeIntegratorWrapper<T>>(
        time_integrator_);
  }

  void UpdateConservatively(::amrex::MultiFab& dest,
                            const ::amrex::MultiFab& src,
                            const ::amrex::MultiFab& fluxes,
                            const ::amrex::Geometry& geom, Duration dt,
                            Direction dir) override {
    time_integrator_.UpdateConservatively(dest, src, fluxes, geom, dt, dir);
  }

  T time_integrator_;
};

class HyperbolicSplitTimeIntegrator {
public:
  HyperbolicSplitTimeIntegrator() = delete;

  HyperbolicSplitTimeIntegrator(const HyperbolicSplitTimeIntegrator& other);
  HyperbolicSplitTimeIntegrator&
  operator=(const HyperbolicSplitTimeIntegrator& other);

  HyperbolicSplitTimeIntegrator(HyperbolicSplitTimeIntegrator&&) = default;
  HyperbolicSplitTimeIntegrator&
  operator=(HyperbolicSplitTimeIntegrator&&) = default;

  template <typename I>
  HyperbolicSplitTimeIntegrator(I&& integrator)
      : time_integrator_{std::make_unique<
            HyperbolicSplitTimeIntegratorWrapper<std::decay_t<I>>>(
            std::forward<I>(integrator))} {}

  void UpdateConservatively(::amrex::MultiFab& dest,
                            const ::amrex::MultiFab& src,
                            const ::amrex::MultiFab& fluxes,
                            const ::amrex::Geometry& geom, Duration dt,
                            Direction dir);

private:
  std::unique_ptr<HyperbolicSplitTimeIntegratorBase> time_integrator_;
};

struct ForwardIntegrator {
  void UpdateConservatively(::amrex::MultiFab& dest,
                            const ::amrex::MultiFab& src,
                            const ::amrex::MultiFab& fluxes,
                            const ::amrex::Geometry& geom, Duration dt,
                            Direction dir);
};

inline HyperbolicSplitTimeIntegrator::HyperbolicSplitTimeIntegrator(
    const HyperbolicSplitTimeIntegrator& other)
    : time_integrator_{other.time_integrator_->Clone()} {}

inline HyperbolicSplitTimeIntegrator& HyperbolicSplitTimeIntegrator::
operator=(const HyperbolicSplitTimeIntegrator& other) {
  time_integrator_ = other.time_integrator_->Clone();
  return *this;
}

inline void HyperbolicSplitTimeIntegrator::UpdateConservatively(
    ::amrex::MultiFab& dest, const ::amrex::MultiFab& src,
    const ::amrex::MultiFab& fluxes, const ::amrex::Geometry& geom, Duration dt,
    Direction dir) {
  time_integrator_->UpdateConservatively(dest, src, fluxes, geom, dt, dir);
}

} // namespace fub::amrex

#endif
