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

#ifndef FUB_EQUATION_HPP
#define FUB_EQUATION_HPP

#include "fub/Direction.hpp"
#include "fub/State.hpp"
#include "fub/StateArray.hpp"
#include "fub/core/type_traits.hpp"

#include "fub/ext/Eigen.hpp"

namespace fub {

template <typename Eq, typename... Args>
using FluxT = decltype(std::declval<Eq>().Flux(std::declval<Args>()...));

template <typename Eq>
using ScalarFluxT =
    decltype(std::declval<const Eq&>().Flux(std::declval<Conservative<Eq>&>(),
                                            std::declval<const Complete<Eq>&>(),
                                            Direction::X));

template <typename Eq, typename N = int_constant<kDefaultChunkSize>>
using VectorizedFluxT = decltype(std::declval<const Eq&>().Flux(
    std::declval<ConservativeArray<Eq, N::value>&>(),
    std::declval<const CompleteArray<Eq, N::value>&>(), Direction::X));

template <typename Eq>
using ScalarReconstructT = decltype(std::declval<const Eq&>().Reconstruct(
    std::declval<Complete<Eq>&>(), std::declval<const Conservative<Eq>&>()));

template <typename Eq, typename N = int_constant<kDefaultChunkSize>>
using VectorizedReconstructT = decltype(std::declval<const Eq&>().Reconstruct(
    std::declval<CompleteArray<Eq, N::value>&>(),
    std::declval<const ConservativeArray<Eq, N::value>&>()));

template <typename Equation>
struct HasScalarFlux : is_detected<ScalarFluxT, Equation> {};

template <typename Equation, int N = kDefaultChunkSize>
struct HasVectorizedFlux
    : is_detected<VectorizedFluxT, Equation, int_constant<N>> {};

template <typename Equation>
struct HasScalarReconstruction : is_detected<ScalarReconstructT, Equation> {};

template <typename Equation, int N = kDefaultChunkSize>
struct HasVectorizedReconstruction
    : is_detected<VectorizedReconstructT, Equation, int_constant<N>> {};

template <typename Equation>
struct HasReconstruction : disjunction<HasScalarReconstruction<Equation>,
                                       HasVectorizedReconstruction<Equation>> {
};

template <typename Extent>
auto Shrink(const layout_left::mapping<Extent>& layout, Direction dir,
            std::ptrdiff_t n = 1) {
  const std::array<std::ptrdiff_t, Extent::rank()> extents =
      AsArray(layout.extents());
  return layout_left::mapping<Extent>(Extent(Shift(extents, dir, -n)));
}

template <typename State, typename Layout, int Rank>
View<State> Subview(const BasicView<State, Layout, Rank>& state,
                    const IndexBox<Rank>& box) {
  auto subview = [](const auto& pdv, const IndexBox<Rank>& box) {
    using PatchDataView = std::decay_t<decltype(pdv)>;
    if constexpr (PatchDataView::Rank() == Rank) {
      return pdv.Subview(box);
    } else if constexpr (PatchDataView::Rank() == Rank + 1) {
      const std::ptrdiff_t lower = pdv.Box().lower[Rank];
      const std::ptrdiff_t upper = pdv.Box().upper[Rank];
      const IndexBox<Rank + 1> embedded_box =
          Embed<Rank + 1>(box, {lower, upper});
      return pdv.Subview(embedded_box);
    }
  };
  View<State> strided{};
  ForEachVariable(
      [&](auto& dest, const auto& src) { dest = subview(src, box); }, strided,
      state);
  return strided;
}

template <typename Eq, typename Equation = std::decay_t<Eq>>
void Flux(Eq&& equation, Conservative<Equation>& flux,
          const Complete<Equation>& state, Direction dir,
          [[maybe_unused]] double x = 0.0) {
  if constexpr (is_detected<FluxT, Eq, Conservative<Equation>&,
                            const Complete<Equation>&, Direction,
                            double>::value) {
    equation.Flux(flux, state, dir, x);
  } else if constexpr (is_detected_exact<Conservative<Equation>, FluxT, Eq,
                                         const Complete<Equation>&, Direction,
                                         double>::value) {
    flux = equation.Flux(state, dir, x);
  } else if constexpr (is_detected_exact<Conservative<Equation>, FluxT, Eq,
                                         const Complete<Equation>&,
                                         Direction>::value) {
    flux = equation.Flux(state, dir);
  } else {
    static_assert(is_detected<FluxT, Eq, Conservative<Equation>&,
                              const Complete<Equation>&, Direction>::value);
    equation.Flux(flux, state, dir);
  }
}

template <typename Eq, typename Equation = std::decay_t<Eq>>
void Flux(Eq&& equation, ConservativeArray<Equation>& flux,
          const CompleteArray<Equation>& state, Direction dir,
          [[maybe_unused]] double x = 0.0) {
  if constexpr (is_detected<FluxT, Eq, ConservativeArray<Equation>&,
                            const CompleteArray<Equation>&, Direction,
                            double>::value) {
    equation.Flux(flux, state, dir, x);
  } else if constexpr (is_detected_exact<ConservativeArray<Equation>, FluxT, Eq,
                                         const CompleteArray<Equation>&,
                                         Direction, double>::value) {
    flux = equation.Flux(state, dir, x);
  } else if constexpr (is_detected_exact<ConservativeArray<Equation>, FluxT, Eq,
                                         const CompleteArray<Equation>&,
                                         Direction>::value) {
    flux = equation.Flux(state, dir);
  } else if constexpr (is_detected<FluxT, Eq, ConservativeArray<Equation>&,
                                   const CompleteArray<Equation>&,
                                   Direction>::value) {
    equation.Flux(flux, state, dir);
  }
}

template <typename Eq, typename Equation = std::decay_t<Eq>>
void Flux(Eq&& equation, ConservativeArray<Equation>& flux,
          const CompleteArray<Equation>& state, MaskArray mask, Direction dir,
          [[maybe_unused]] double x = 0.0) {
  if constexpr (is_detected<FluxT, Eq, ConservativeArray<Equation>&,
                            const CompleteArray<Equation>&, MaskArray, Direction,
                            double>::value) {
    equation.Flux(flux, state, mask, dir, x);
  } else if constexpr (is_detected_exact<ConservativeArray<Equation>, FluxT, Eq,
                                         const CompleteArray<Equation>&, MaskArray,
                                         Direction, double>::value) {
    flux = equation.Flux(state, mask, dir, x);
  } else if constexpr (is_detected_exact<ConservativeArray<Equation>, FluxT, Eq,
                                         const CompleteArray<Equation>&, MaskArray,
                                         Direction>::value) {
    flux = equation.Flux(state, mask, dir);
  } else if constexpr (is_detected<FluxT, Eq, ConservativeArray<Equation>&,
                                   const CompleteArray<Equation>&, MaskArray,
                                   Direction>::value) {
    equation.Flux(flux, state, mask, dir);
  }
}

} // namespace fub

#endif
