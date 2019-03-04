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
#include "fub/State.hpp"

#include <array>

namespace fub {

template <typename Method> struct FluxMethodTraits {
  using Equation = typename Method::Equation;
  using Complete = typename Method::Complete;
  using Conservative = typename Method::Conservative;
  static constexpr int StencilWidth() noexcept {
    return Method::GetStencilWidth();
  }
};

/// This class applies a Base Method on a an array of states.
///
/// The Base class only needs to define its logic on single states instead
/// of a multi dimensional arrray. By using this class the developer does
/// not need to repeat that logic.
template <typename BaseMethod> class FluxMethod {
public:
  using Equation = typename FluxMethodTraits<BaseMethod>::Equation;
  using Complete = typename FluxMethodTraits<BaseMethod>::Complete;
  using Conservative = typename FluxMethodTraits<BaseMethod>::Conservative;

  /// Allocates memory space for properly sized stencil and numeric flux
  /// buffers.
  template <typename... Args>
  FluxMethod(Args&&... args) : base_(std::forward<Args>(args)...) {}

  const BaseMethod& Base() { return base_; }

  const Equation& GetEquation() const noexcept { return base_.GetEquation(); }

  static constexpr int GetStencilWidth() noexcept {
    return FluxMethodTraits<BaseMethod>::StencilWidth();
  }

  void ComputeNumericFluxes(View<Conservative> fluxes,
                            View<const Complete> states, Duration dt, double dx,
                            Direction dir) {
    ForEachIndex(Mapping<0>(fluxes), [&](const auto... is) {
      using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
      const Index face{is...};
      for (std::size_t i = 0; i < stencil_.size(); ++i) {
        const Index cell = Shift(face, dir, i);
        Load(stencil_[i], states, cell);
      }
      base_.ComputeNumericFlux(numeric_flux_, stencil_, dt, dx, dir);
      Store(fluxes, numeric_flux_, face);
    });
  }

  void ComputeNumericFlux(Conservative& numeric_flux,
                          span<const Complete, 2> states, Duration dt,
                          double dx, Direction dir) {
    base_.ComputeNumericFlux(numeric_flux, states, dt, dx, dir);
  }

  double ComputeStableDt(span<const Complete, 2> states, double dx,
                         Direction dir) {
    return base_.ComputeStableDt(states, dx, dir);
  }

  double ComputeStableDt(View<const Complete> states, double dx,
                         Direction dir) {
    double min_dt = std::numeric_limits<double>::infinity();
    constexpr int stencil = FluxMethodTraits<BaseMethod>::StencilWidth();
    ForEachIndex(Shrink(Mapping<0>(states), dir, 2 * stencil),
                 [&](const auto... is) {
                   using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
                   const Index face{is...};
                   for (std::size_t i = 0; i < stencil_.size(); ++i) {
                     const Index cell = Shift(face, dir, i);
                     Load(stencil_[i], states, cell);
                   }
                   double dt = base_.ComputeStableDt(stencil_, dx, dir);
                   min_dt = std::min(dt, min_dt);
                 });
    return min_dt;
  }

  BaseMethod base_;
  std::array<Complete, 2 * FluxMethodTraits<BaseMethod>::StencilWidth()>
      stencil_{};
  Conservative numeric_flux_{};
};

} // namespace fub

#endif