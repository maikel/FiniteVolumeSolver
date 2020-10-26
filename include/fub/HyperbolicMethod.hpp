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

#ifndef FUB_HYPERBOLIC_METHOD_HPP
#define FUB_HYPERBOLIC_METHOD_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/core/type_traits.hpp"

#include <memory>
#include <type_traits>

namespace fub {

/// \defgroup PolymorphicValueType Polymorphic Value Types
/// \brief This group summarizes all types which have a polymorphic behaviour.

///////////////////////////////////////////////////////////////////////////////
//                                                          Strategy Interfaces

namespace detail {
template <typename IntegratorContext> struct ReconstructionBase {
  virtual ~ReconstructionBase() = default;
  virtual std::unique_ptr<ReconstructionBase> Clone() const = 0;
  virtual void CompleteFromCons(IntegratorContext& context, int level,
                                Duration dt) = 0;
};

template <typename IntegratorContext> struct TimeIntegratorBase {
  virtual ~TimeIntegratorBase() = default;
  virtual std::unique_ptr<TimeIntegratorBase> Clone() const = 0;
  virtual void UpdateConservatively(IntegratorContext& context, int level,
                                    Duration dt, Direction dir) = 0;
};

template <typename IntegratorContext> struct FluxMethodBase {
  virtual ~FluxMethodBase() = default;
  virtual std::unique_ptr<FluxMethodBase> Clone() const = 0;
  virtual void PreAdvanceHierarchy(IntegratorContext& context) = 0;
  virtual Duration ComputeStableDt(IntegratorContext& context, int level,
                                   Direction dir) = 0;
  virtual void ComputeNumericFluxes(IntegratorContext& context, int level,
                                    Duration dt, Direction dir) = 0;
  virtual int GetStencilWidth() const = 0;
};
} // namespace detail

///////////////////////////////////////////////////////////////////////////////
//                                                          Polymorphic Classes

/// \ingroup PolymorphicValueType
/// This is a polymorphic wrapper class for reconstruction strategies used by
/// \a IntegratorContext.
template <typename IntegratorContext> class AnyReconstruction {
public:
  AnyReconstruction() = delete;

  AnyReconstruction(const AnyReconstruction& other);
  AnyReconstruction& operator=(const AnyReconstruction& other);

  AnyReconstruction(AnyReconstruction&& other) noexcept = default;
  AnyReconstruction& operator=(AnyReconstruction&& other) noexcept = default;

  template <typename R,
            typename = std::enable_if_t<!decays_to<R, AnyReconstruction>()>>
  AnyReconstruction(R&& r); // NOLINT

  void CompleteFromCons(IntegratorContext& context, int level, Duration dt);

private:
  std::unique_ptr<detail::ReconstructionBase<IntegratorContext>> reconstruct_;
};

/// \ingroup PolymorphicValueType
/// This is a polymorphic wrapper class for TimeIntegator strategies used by
/// \a IntegratorContext.
template <typename IntegratorContext> class AnyTimeIntegrator {
public:
  AnyTimeIntegrator() = delete;

  AnyTimeIntegrator(const AnyTimeIntegrator& other);
  AnyTimeIntegrator& operator=(const AnyTimeIntegrator& other);

  AnyTimeIntegrator(AnyTimeIntegrator&& other) noexcept = default;
  AnyTimeIntegrator& operator=(AnyTimeIntegrator&& other) noexcept = default;

  template <typename R,
            typename = std::enable_if_t<!decays_to<R, AnyTimeIntegrator>()>>
  AnyTimeIntegrator(R&& r); // NOLINT

  void UpdateConservatively(IntegratorContext& context, int level, Duration dt,
                            Direction dir);

private:
  std::unique_ptr<detail::TimeIntegratorBase<IntegratorContext>> integrator_;
};

/// \ingroup PolymorphicValueType FluxMethod
/// This is a polymorphic wrapper class for FluxMethod strategies used by
/// \a IntegratorContext.
template <typename IntegratorContext> class AnyFluxMethod {
public:
  AnyFluxMethod() = delete;

  AnyFluxMethod(const AnyFluxMethod& other);
  AnyFluxMethod& operator=(const AnyFluxMethod& other);

  AnyFluxMethod(AnyFluxMethod&& other) noexcept = default;
  AnyFluxMethod& operator=(AnyFluxMethod&& other) noexcept = default;

  template <typename R,
            typename = std::enable_if_t<!decays_to<R, AnyFluxMethod>()>>
  AnyFluxMethod(R&& r); // NOLINT

  Duration ComputeStableDt(IntegratorContext& context, int level,
                           Direction dir);

  void PreAdvanceHierarchy(IntegratorContext& context);

  void ComputeNumericFluxes(IntegratorContext& context, int level, Duration dt,
                            Direction dir);

  int GetStencilWidth() const;

private:
  std::unique_ptr<detail::FluxMethodBase<IntegratorContext>> flux_method_;
};

///////////////////////////////////////////////////////////////////////////////
//                                                     Collection of Strategies

template <typename IntegratorContext> struct HyperbolicMethod {
  AnyFluxMethod<IntegratorContext> flux_method;
  AnyTimeIntegrator<IntegratorContext> time_integrator;
  AnyReconstruction<IntegratorContext> reconstruction;
};

///////////////////////////////////////////////////////////////////////////////
//                                                               IMPLEMENTATION

///////////////////////////////////////////////////////////////////////////////
//                                                               Reconstruction

template <typename IntegratorContext>
AnyReconstruction<IntegratorContext>::AnyReconstruction(
    const AnyReconstruction& other)
    : reconstruct_{other.reconstruct_ ? other.reconstruct_->Clone() : nullptr} {
}

template <typename IntegratorContext>
AnyReconstruction<IntegratorContext>&
AnyReconstruction<IntegratorContext>::operator=(
    const AnyReconstruction& other) {
  if (other.reconstruct_) {
    reconstruct_ = other.reconstruct_->Clone();
  } else {
    reconstruct_ = nullptr;
  }
  return *this;
}

template <typename IntegratorContext>
void AnyReconstruction<IntegratorContext>::CompleteFromCons(
    IntegratorContext& context, int level, Duration dt) {
  reconstruct_->CompleteFromCons(context, level, dt);
}

namespace detail {
template <typename IntegratorContext, typename R>
struct ReconstructionWrapper : ReconstructionBase<IntegratorContext> {
  ReconstructionWrapper(const R& rec) : rec_{rec} {}
  ReconstructionWrapper(R&& rec) : rec_{std::move(rec)} {}
  std::unique_ptr<ReconstructionBase<IntegratorContext>>
  Clone() const override {
    return std::make_unique<ReconstructionWrapper>(rec_);
  }
  void CompleteFromCons(IntegratorContext& context, int level,
                        Duration dt) override {
    rec_.CompleteFromCons(context, level, dt);
  }
  R rec_;
};
} // namespace detail

template <typename IntegratorContext>
template <typename R, typename>
AnyReconstruction<IntegratorContext>::AnyReconstruction(R&& rec)
    : reconstruct_(
          std::make_unique<detail::ReconstructionWrapper<IntegratorContext,
                                                         std::decay_t<R>>>(
              std::forward<R>(rec))) {}

///////////////////////////////////////////////////////////////////////////////
//                                                               TimeIntegrator

template <typename IntegratorContext>
AnyTimeIntegrator<IntegratorContext>::AnyTimeIntegrator(
    const AnyTimeIntegrator& other)
    : integrator_{other.integrator_ ? other.integrator_->Clone() : nullptr} {}

template <typename IntegratorContext>
AnyTimeIntegrator<IntegratorContext>&
AnyTimeIntegrator<IntegratorContext>::operator=(
    const AnyTimeIntegrator& other) {
  if (other.integrator_) {
    integrator_ = other.integrator_->Clone();
  } else {
    integrator_ = nullptr;
  }
  return *this;
}

template <typename IntegratorContext>
void AnyTimeIntegrator<IntegratorContext>::UpdateConservatively(
    IntegratorContext& context, int level, Duration dt, Direction dir) {
  integrator_->UpdateConservatively(context, level, dt, dir);
}

namespace detail {
template <typename IntegratorContext, typename I>
struct TimeIntegratorWrapper : TimeIntegratorBase<IntegratorContext> {
  TimeIntegratorWrapper(const I& integrator) : integrator_{integrator} {}
  TimeIntegratorWrapper(I&& integrator) : integrator_{std::move(integrator)} {}
  std::unique_ptr<TimeIntegratorBase<IntegratorContext>>
  Clone() const override {
    return std::make_unique<TimeIntegratorWrapper>(integrator_);
  }
  void UpdateConservatively(IntegratorContext& context, int level, Duration dt,
                            Direction dir) override {
    integrator_.UpdateConservatively(context, level, dt, dir);
  }
  I integrator_;
};
} // namespace detail

template <typename IntegratorContext>
template <typename I, typename>
AnyTimeIntegrator<IntegratorContext>::AnyTimeIntegrator(I&& integrator)
    : integrator_(
          std::make_unique<detail::TimeIntegratorWrapper<IntegratorContext,
                                                         std::decay_t<I>>>(
              std::forward<I>(integrator))) {}

///////////////////////////////////////////////////////////////////////////////
//                                                                   FluxMethod

template <typename IntegratorContext>
AnyFluxMethod<IntegratorContext>::AnyFluxMethod(const AnyFluxMethod& other)
    : flux_method_{other.flux_method_ ? other.flux_method_->Clone() : nullptr} {
}

template <typename IntegratorContext>
AnyFluxMethod<IntegratorContext>&
AnyFluxMethod<IntegratorContext>::operator=(const AnyFluxMethod& other) {
  if (other.flux_method_) {
    flux_method_ = other.flux_method_->Clone();
  } else {
    flux_method_ = nullptr;
  }
  return *this;
}

template <typename IntegratorContext>
void AnyFluxMethod<IntegratorContext>::PreAdvanceHierarchy(
    IntegratorContext& context) {
  flux_method_->PreAdvanceHierarchy(context);
}

template <typename IntegratorContext>
Duration
AnyFluxMethod<IntegratorContext>::ComputeStableDt(IntegratorContext& context,
                                                  int level, Direction dir) {
  return flux_method_->ComputeStableDt(context, level, dir);
}

template <typename IntegratorContext>
void AnyFluxMethod<IntegratorContext>::ComputeNumericFluxes(
    IntegratorContext& context, int level, Duration dt, Direction dir) {
  flux_method_->ComputeNumericFluxes(context, level, dt, dir);
}

template <typename IntegratorContext>
int AnyFluxMethod<IntegratorContext>::GetStencilWidth() const {
  return flux_method_->GetStencilWidth();
}

namespace detail {
template <typename I, typename... Args>
using PreAdvanceHierarchyT =
    decltype(std::declval<I>().PreAdvanceHierarchy(std::declval<Args>()...));

template <typename IntegratorContext, typename I>
struct FluxMethodWrapper : FluxMethodBase<IntegratorContext> {
  FluxMethodWrapper(const I& flux_method) : flux_method_{flux_method} {}
  FluxMethodWrapper(I&& flux_method) : flux_method_{std::move(flux_method)} {}
  std::unique_ptr<FluxMethodBase<IntegratorContext>> Clone() const override {
    return std::make_unique<FluxMethodWrapper>(flux_method_);
  }
  Duration ComputeStableDt(IntegratorContext& context, int level,
                           Direction dir) override {
    return flux_method_.ComputeStableDt(context, level, dir);
  }
  void ComputeNumericFluxes(IntegratorContext& context, int level, Duration dt,
                            Direction dir) override {
    flux_method_.ComputeNumericFluxes(context, level, dt, dir);
  }
  int GetStencilWidth() const override {
    return flux_method_.GetStencilWidth();
  }
  void PreAdvanceHierarchy(IntegratorContext& context) override {
    if constexpr (is_detected<PreAdvanceHierarchyT, I&,
                              IntegratorContext&>::value) {
      flux_method_.PreAdvanceHierarchy(context);
    }
  }
  I flux_method_;
};
} // namespace detail

template <typename IntegratorContext>
template <typename I, typename>
AnyFluxMethod<IntegratorContext>::AnyFluxMethod(I&& flux_method)
    : flux_method_(
          std::make_unique<
              detail::FluxMethodWrapper<IntegratorContext, std::decay_t<I>>>(
              std::forward<I>(flux_method))) {}
} // namespace fub
#endif
