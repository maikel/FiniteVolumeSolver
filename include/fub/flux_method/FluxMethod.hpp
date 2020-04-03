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

#ifndef FUB_FLUX_METHOD_HPP
#define FUB_FLUX_METHOD_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/Execution.hpp"
#include "fub/ForEach.hpp"
#include "fub/Meta.hpp"
#include "fub/State.hpp"
#include "fub/StateArray.hpp"
#include "fub/StateRow.hpp"

#include <array>

namespace fub {

/// \defgroup FluxMethod Flux Methods
/// \ingroup IntegratorContext
/// \brief This module collects all types and functions relevant to compute
/// numerical fluxes


/// \ingroup FluxMethod
/// \brief This class applies a base flux nethod on a view of states.
///
/// The base class only needs to define its logic on single states or state
/// arrays instead of on a multi dimensional arrray. By using this class the
/// developer does not need to repeat the logic which might involve iterating
/// through all indices of a patch, loading the states and applying the base
/// method.
///
/// There are in total two strategies implemented, a cell-wise and a simdified
/// one. The sequential strategy is easier to debug but less performant. The
/// simd version does use the spatial grid index. This class also assumes that
/// base class does not use the local coordinates.
template <typename BaseMethod> class FluxMethod : public BaseMethod {
public:
  using Equation = meta::Equation<const BaseMethod&>;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;

  static constexpr int Rank = Equation::Rank();

  /// \brief Returns the stencil width of this flux method.
  static constexpr int GetStencilWidth() noexcept;

  /// \brief Returns the number of elements in a stencil of this flux method.
  static constexpr int GetStencilSize() noexcept;

  /// \brief Returns the number of elements in a stencil of this flux method
  /// (std::size_t version).
  static constexpr std::size_t StencilSize =
      static_cast<std::size_t>(GetStencilSize());

  FluxMethod() = default;
  ~FluxMethod() = default;

  FluxMethod(const FluxMethod&) = default;
  FluxMethod(FluxMethod&&) = default;

  FluxMethod& operator=(const FluxMethod&) = default;
  FluxMethod& operator=(FluxMethod&&) = default;

  template <typename... Args> FluxMethod(Args&&... args);

  /// @{
  /// \brief Returns the Implementation class which will be used to compute
  /// single fluxes.
  const BaseMethod& Base() const noexcept;
  BaseMethod& Base() noexcept;
  /// @}

  using BaseMethod::GetEquation;

  /// @{
  /// \brief This function computes numerical fluxes.
  ///
  /// \param[out] fluxes  The destination view where numeric fluxes will be
  ///                     stored at.
  ///
  /// \param[in] states  A view over cell states which will be used to fill the
  ///                    stencil for this method.
  ///
  /// \param[in] dir  The direction parameter specifies which directional flux
  ///                 method to use.
  ///
  /// \param[in] dt  The time step size which will be used to compute the
  ///                numerical flux.
  ///
  /// \param[in] dx  The cell width size which will be used to compute the
  ///                numerical flux.
  void ComputeNumericFluxes(const View<Conservative>& fluxes,
                            const View<const Complete>& states, Duration dt,
                            double dx, Direction dir);

  void ComputeNumericFluxes(execution::SequentialTag,
                            const View<Conservative>& fluxes,
                            const View<const Complete>& states, Duration dt,
                            double dx, Direction dir);

  void ComputeNumericFluxes(execution::SimdTag,
                            const View<Conservative>& fluxes,
                            const View<const Complete>& states, Duration dt,
                            double dx, Direction dir);

  void ComputeNumericFluxes(execution::OpenMpTag,
                            const View<Conservative>& fluxes,
                            const View<const Complete>& states, Duration dt,
                            double dx, Direction dir);

  void ComputeNumericFluxes(execution::OpenMpSimdTag,
                            const View<Conservative>& fluxes,
                            const View<const Complete>& states, Duration dt,
                            double dx, Direction dir);
  /// @}

  using BaseMethod::ComputeNumericFlux;

  /// @{
  /// \brief This function computes a time step size such that no signal will
  /// leave any cell covered by this view.
  ///
  /// \param[in] states  A view over cell states which will be used to fill the
  ///                    stencil for this method.
  ///
  /// \param[in] dir  The direction parameter specifies which directional flux
  ///                 method to use.
  ///
  /// \param[in] dx  The cell width size which will be used to compute the
  ///                numerical flux.
  double ComputeStableDt(const View<const Complete>& states, double dx,
                         Direction dir);

  double ComputeStableDt(execution::SequentialTag,
                         const View<const Complete>& states, double dx,
                         Direction dir);

  double ComputeStableDt(execution::SimdTag, const View<const Complete>& states,
                         double dx, Direction dir);

  double ComputeStableDt(execution::OpenMpTag,
                         const View<const Complete>& states, double dx,
                         Direction dir);

  double ComputeStableDt(execution::OpenMpSimdTag,
                         const View<const Complete>& states, double dx,
                         Direction dir);
  /// @}

  using BaseMethod::ComputeStableDt;

  std::array<Complete, StencilSize> stencil_{};
  Conservative numeric_flux_{GetEquation()};

  std::array<CompleteArray, StencilSize> stencil_array_{};
  ConservativeArray numeric_flux_array_{GetEquation()};
};

template <typename BaseMethod>
constexpr int FluxMethod<BaseMethod>::GetStencilWidth() noexcept {
  return BaseMethod::GetStencilWidth();
}

template <typename BaseMethod>
constexpr int FluxMethod<BaseMethod>::GetStencilSize() noexcept {
  return 2 * GetStencilWidth();
}

template <typename BaseMethod>
template <typename... Args>
FluxMethod<BaseMethod>::FluxMethod(Args&&... args)
    : BaseMethod(std::forward<Args>(args)...) {
  stencil_.fill(Complete{BaseMethod::GetEquation()});
  stencil_array_.fill(CompleteArray{BaseMethod::GetEquation()});
}

template <typename BaseMethod>
void FluxMethod<BaseMethod>::ComputeNumericFluxes(
    execution::SequentialTag, const View<Conservative>& fluxes,
    const View<const Complete>& states, Duration dt, double dx, Direction dir) {
  ForEachIndex(Box<0>(fluxes), [&](const auto... is) {
    using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
    const Index face{is...};
    const Index cell0 = Shift(face, dir, -GetStencilWidth());
    for (std::size_t i = 0; i < stencil_.size(); ++i) {
      const Index cell = Shift(cell0, dir, static_cast<int>(i));
      Load(stencil_[i], states, cell);
    }
    BaseMethod::ComputeNumericFlux(numeric_flux_, stencil_, dt, dx, dir);
    Store(fluxes, numeric_flux_, face);
  });
}

template <typename BaseMethod>
void FluxMethod<BaseMethod>::ComputeNumericFluxes(
    const View<Conservative>& fluxes, const View<const Complete>& states,
    Duration dt, double dx, Direction dir) {
  ComputeNumericFluxes(execution::seq, fluxes, states, dt, dx, dir);
}

template <typename BaseMethod>
void FluxMethod<BaseMethod>::ComputeNumericFluxes(
    execution::OpenMpTag, const View<Conservative>& fluxes,
    const View<const Complete>& states, Duration dt, double dx, Direction dir) {
  ComputeNumericFluxes(execution::seq, fluxes, states, dt, dx, dir);
}

template <typename BaseMethod>
void FluxMethod<BaseMethod>::ComputeNumericFluxes(
    execution::OpenMpSimdTag, const View<Conservative>& fluxes,
    const View<const Complete>& states, Duration dt, double dx, Direction dir) {
  ComputeNumericFluxes(execution::simd, fluxes, states, dt, dx, dir);
}

template <typename FM, typename... Args>
using HasComputeNumericFluxType =
    decltype(std::declval<FM>().ComputeNumericFlux(std::declval<Args>()...));

template <typename FM, typename... Args>
using HasComputeNumericFlux =
    is_detected<HasComputeNumericFluxType, FM, Args...>;

template <typename BaseMethod>
void FluxMethod<BaseMethod>::ComputeNumericFluxes(
    execution::SimdTag, const View<Conservative>& fluxes,
    const View<const Complete>& states, Duration dt, double dx, Direction dir) {
  static_assert(HasComputeNumericFlux<BaseMethod&, ConservativeArray&,
                                      std::array<CompleteArray, StencilSize>,
                                      Duration, double, Direction>());
  IndexBox<Rank> fluxbox = Box<0>(fluxes);
  static constexpr int kWidth = GetStencilWidth();
  IndexBox<Rank> cellbox = Grow(fluxbox, dir, {kWidth, kWidth - 1});
  BasicView base = Subview(states, cellbox);
  std::array<View<const Complete>, StencilSize> stencil_views;
  for (std::size_t i = 0; i < StencilSize; ++i) {
    stencil_views[i] =
        Shrink(base, dir,
               {static_cast<std::ptrdiff_t>(i),
                static_cast<std::ptrdiff_t>(StencilSize - i) - 1});
  }
  std::tuple views = std::apply(
      [&fluxes](const auto&... vs) {
        return std::tuple{fluxes, vs...};
      },
      stencil_views);
  ForEachRow(views, [this, dt, dx, dir](const Row<Conservative>& fluxes,
                                        auto... rows) {
    ViewPointer fit = Begin(fluxes);
    ViewPointer fend = End(fluxes);
    std::array<ViewPointer<const Complete>, StencilSize> states{Begin(rows)...};
    int n = static_cast<int>(get<0>(fend) - get<0>(fit));
    while (n >= kDefaultChunkSize) {
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Load(stencil_array_[i], states[i]);
      }
      ComputeNumericFlux(numeric_flux_array_, stencil_array_, dt, dx, dir);
      Store(fit, numeric_flux_array_);
      Advance(fit, kDefaultChunkSize);
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Advance(states[i], kDefaultChunkSize);
      }
      n = static_cast<int>(get<0>(fend) - get<0>(fit));
    }
    for (std::size_t i = 0; i < StencilSize; ++i) {
      LoadN(stencil_array_[i], states[i], n);
    }
    ComputeNumericFlux(numeric_flux_array_, stencil_array_, dt, dx, dir);
    StoreN(fit, numeric_flux_array_, n);
  });
}

template <typename BaseMethod>
double
FluxMethod<BaseMethod>::ComputeStableDt(execution::SequentialTag,
                                        const View<const Complete>& states,
                                        double dx, Direction dir) {
  double min_dt = std::numeric_limits<double>::infinity();
  ForEachIndex(Shrink(Box<0>(states), dir, {0, GetStencilSize()}),
               [&](const auto... is) {
                 using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
                 const Index face{is...};
                 for (std::size_t i = 0; i < stencil_.size(); ++i) {
                   const Index cell = Shift(face, dir, static_cast<int>(i));
                   Load(stencil_[i], states, cell);
                 }
                 const double dt =
                     BaseMethod::ComputeStableDt(stencil_, dx, dir);
                 min_dt = std::min(dt, min_dt);
               });
  return min_dt;
}

template <typename BaseMethod>
double
FluxMethod<BaseMethod>::ComputeStableDt(const View<const Complete>& states,
                                        double dx, Direction dir) {
  return ComputeStableDt(execution::seq, states, dx, dir);
}

template <typename BaseMethod>
double
FluxMethod<BaseMethod>::ComputeStableDt(execution::OpenMpTag,
                                        const View<const Complete>& states,
                                        double dx, Direction dir) {
  return ComputeStableDt(execution::seq, states, dx, dir);
}

template <typename BaseMethod>
double
FluxMethod<BaseMethod>::ComputeStableDt(execution::OpenMpSimdTag,
                                        const View<const Complete>& states,
                                        double dx, Direction dir) {
  return ComputeStableDt(execution::simd, states, dx, dir);
}

template <typename T, typename... Ts> T Head(T&& head, Ts&&...) { return head; }

template <typename BaseMethod>
double
FluxMethod<BaseMethod>::ComputeStableDt(execution::SimdTag,
                                        const View<const Complete>& states,
                                        double dx, Direction dir) {
  Array1d min_dt(std::numeric_limits<double>::infinity());
  std::array<View<const Complete>, StencilSize> stencil_views;
  for (std::size_t i = 0; i < StencilSize; ++i) {
    stencil_views[i] =
        Shrink(states, dir,
               {static_cast<std::ptrdiff_t>(i),
                static_cast<std::ptrdiff_t>(StencilSize - i) - 1});
  }
  std::tuple views = std::apply(
      [](const auto&... vs) { return std::tuple{vs...}; }, stencil_views);
  ForEachRow(views, [this, dx, dir, &min_dt](auto... rows) {
    std::array<ViewPointer<const Complete>, StencilSize> states{Begin(rows)...};
    ViewPointer<const Complete> send = End(Head(rows...));
    int n = static_cast<int>(get<0>(send) - get<0>(states[0]));
    while (n >= kDefaultChunkSize) {
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Load(stencil_array_[i], states[i]);
      }
      min_dt = min_dt.min(ComputeStableDt(stencil_array_, dx, dir));
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Advance(states[i], kDefaultChunkSize);
      }
      n = static_cast<int>(get<0>(send) - get<0>(states[0]));
    }
    for (std::size_t i = 0; i < StencilSize; ++i) {
      LoadN(stencil_array_[i], states[i], n);
    }
    Array1d mask_;
    for (int i = 0; i < mask_.rows(); ++i) {
      mask_[i] = static_cast<double>(i);
    }
    auto mask = (mask_ < n).eval();
    min_dt = mask.select(min_dt.min(ComputeStableDt(stencil_array_, dx, dir)),
                         min_dt);
  });
  return min_dt.minCoeff();
}

} // namespace fub

#endif
