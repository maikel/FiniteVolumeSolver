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
#include "fub/State.hpp"

#include <array>

namespace fub {

template <typename Method> struct FluxMethodTraits {
  using Equation = typename Method::Equation;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
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

  using CompleteArray = typename FluxMethodTraits<BaseMethod>::CompleteArray;
  using ConservativeArray =
      typename FluxMethodTraits<BaseMethod>::ConservativeArray;

  static constexpr int ChunkSize = BaseMethod::ChunkSize;

  /// Allocates memory space for properly sized stencil and numeric flux
  /// buffers.
  template <typename... Args>
  FluxMethod(Args&&... args) : base_(std::forward<Args>(args)...) {}

  const BaseMethod& Base() { return base_; }

  const Equation& GetEquation() const noexcept { return base_.GetEquation(); }

  static constexpr int GetStencilWidth() noexcept {
    return FluxMethodTraits<BaseMethod>::StencilWidth();
  }

  void ComputeNumericFluxes(const View<Conservative>& fluxes,
                            const View<const Complete>& states, Direction dir,
                            Duration dt, double dx) {
    ForEachIndex(Box<0>(fluxes), [&](const auto... is) {
      using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
      const Index face{is...};
      const Index cell0 = Shift(face, dir, -GetStencilWidth());
      for (std::size_t i = 0; i < stencil_.size(); ++i) {
        const Index cell = Shift(cell0, dir, static_cast<int>(i));
        Load(stencil_[i], states, cell);
      }
      base_.ComputeNumericFlux(numeric_flux_, stencil_, dt, dx, dir);
      Store(fluxes, numeric_flux_, face);
    });
  }

  void ComputeNumericFluxes(execution::SimdTag,
                            const View<Conservative>& fluxes,
                            const View<const Complete>& states, Direction dir,
                            Duration dt, double dx) {
    static constexpr int stencil = 2 * GetStencilWidth();
    [[maybe_unused]] const int d = static_cast<int>(dir);
    FUB_ASSERT(Extents<0>(states).extent(d) ==
               Extents<0>(fluxes).extent(d) + stencil - 1);
    if (dir == Direction::X) {
      ForEachRow(
          [&](auto fluxes_row, auto states_row) {
            std::ptrdiff_t i = 0;
            const std::ptrdiff_t ssize = Extents<0>(states_row).extent(0);
            const std::ptrdiff_t fsize = Extents<0>(fluxes_row).extent(0);
            while (i + stencil - 1 + ChunkSize <= ssize) {
              for (std::size_t k = 0; k < stencil; ++k) {
                const std::ptrdiff_t index = i + static_cast<std::ptrdiff_t>(k);
                Load(stencil_array_[k], states_row, {index});
              }
              base_.ComputeNumericFlux(numeric_flux_array_, stencil_array_, dt,
                                       dx, dir);
              FUB_ASSERT(i + ChunkSize <= fsize);
              Store(fluxes_row, numeric_flux_array_, {i});
              i += ChunkSize;
            }
            const int rest = ssize - (i + stencil);
            FUB_ASSERT(rest < ChunkSize);
            if (rest > 0) {
              for (std::size_t k = 0; k < stencil; ++k) {
                const std::ptrdiff_t index = i + static_cast<std::ptrdiff_t>(k);
                LoadN(stencil_array_[k], states_row, rest, {index});
              }
              base_.ComputeNumericFlux(numeric_flux_array_, stencil_array_, dt,
                                       dx, dir);
              StoreN(fluxes_row, numeric_flux_array_, rest, {i});
            }
          },
          fluxes, states);
    } else if (dir == Direction::Y) {
      if constexpr (Equation::Rank() == 2) {
        const int n = Extents<0>(states).extent(1);
        std::array<View<const Complete, layout_stride>, stencil>
            stencil_views{};
        for (int i = 0; i < stencil; ++i) {
          stencil_views[i] =
              Slice<Direction::Y>(states, std::pair{i, n - stencil + i});
        }
        std::apply(
            [&](auto... rows) {
              ForEachRow(
                  [&](auto fluxes_row, auto first_row, auto... other_rows) {
                    std::array<decltype(first_row), stencil> rows{
                        first_row, other_rows...};
                    std::ptrdiff_t i = 0;
                    const std::ptrdiff_t ssize =
                        Extents<0>(first_row).extent(0);
                    while (i + ChunkSize <= ssize) {
                      for (std::size_t k = 0; k < stencil; ++k) {
                        Load(stencil_array_[k], rows[k], {i});
                      }
                      base_.ComputeNumericFlux(numeric_flux_array_,
                                               stencil_array_, dt, dx, dir);
                      Store(fluxes_row, numeric_flux_array_, {i});
                      i += ChunkSize;
                    }
                    const int rest = ssize - i;
                    FUB_ASSERT(rest < ChunkSize);
                    if (rest > 0) {
                      for (std::size_t k = 0; k < stencil; ++k) {
                        LoadN(stencil_array_[k], rows[k], rest, {i});
                      }
                      base_.ComputeNumericFlux(numeric_flux_array_,
                                               stencil_array_, dt, dx, dir);
                      StoreN(fluxes_row, numeric_flux_array_, rest, {i});
                    }
                  },
                  fluxes, rows...);
            },
            stencil_views);
      }
    } else {
      if constexpr (Equation::Rank() == 3) {
        FUB_ASSERT(dir == Direction::Z);
        const int n = Extents<0>(states).extent(1);
        std::array<View<const Complete, layout_stride>, stencil>
            stencil_views{};
        for (int i = 0; i < stencil; ++i) {
          stencil_views[i] =
              Slice<Direction::Z>(states, std::pair{i, n - stencil + i});
        }
        std::apply(
            [&](auto... rows) {
              ForEachRow(
                  [&](auto fluxes_row, auto first_row, auto... other_rows) {
                    std::array<decltype(first_row), stencil> rows{
                        first_row, other_rows...};
                    std::ptrdiff_t i = 0;
                    const std::ptrdiff_t size = Extents<0>(first_row).extent(0);
                    while (i + ChunkSize + stencil - 1 <= size) {
                      for (std::size_t k = 0; k < stencil; ++k) {
                        Load(stencil_array_[k], rows[k], {i});
                      }
                      base_.ComputeNumericFlux(numeric_flux_array_,
                                               stencil_array_, dt, dx, dir);
                      Store(fluxes_row, numeric_flux_array_, {i});
                      i += ChunkSize;
                    }
                    const int rest = size - i;
                    FUB_ASSERT(rest < ChunkSize);
                    if (rest > 0) {
                      for (std::size_t k = 0; k < stencil; ++k) {
                        LoadN(stencil_array_[k], rows[k], rest, {i});
                      }
                      base_.ComputeNumericFlux(numeric_flux_array_,
                                               stencil_array_, dt, dx, dir);
                      StoreN(fluxes_row, numeric_flux_array_, rest, {i});
                    }
                  },
                  fluxes, rows...);
            },
            stencil_views);
      }
    }
  }

  void ComputeNumericFlux(Conservative& numeric_flux,
                          span<const Complete, 2 * GetStencilWidth()> states,
                          Duration dt, double dx, Direction dir) {
    base_.ComputeNumericFlux(numeric_flux, states, dt, dx, dir);
  }

  void
  ComputeNumericFlux(ConservativeArray& numeric_flux,
                     span<const CompleteArray, 2 * GetStencilWidth()> states,
                     Duration dt, double dx, Direction dir) {
    base_.ComputeNumericFlux(numeric_flux, states, dt, dx, dir);
  }

  double ComputeStableDt(span<const Complete, 2 * GetStencilWidth()> states,
                         double dx, Direction dir) {
    return base_.ComputeStableDt(states, dx, dir);
  }

  double ComputeStableDt(View<const Complete> states, double dx,
                         Direction dir) {
    double min_dt = std::numeric_limits<double>::infinity();
    constexpr int stencil = FluxMethodTraits<BaseMethod>::StencilWidth();
    ForEachIndex(Shrink(Box<0>(states), dir, {0, 2 * stencil}),
                 [&](const auto... is) {
                   using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
                   const Index face{is...};
                   for (std::size_t i = 0; i < stencil_.size(); ++i) {
                     const Index cell = Shift(face, dir, static_cast<int>(i));
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

  std::array<CompleteArray, 2 * FluxMethodTraits<BaseMethod>::StencilWidth()>
      stencil_array_{};
  ConservativeArray numeric_flux_array_{};
};

} // namespace fub

#endif
