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
#include "fub/core/type_traits.hpp"

namespace fub {

template <typename Eq>
using ScalarFluxT =
    decltype(std::declval<const Eq&>().Flux(std::declval<Conservative<Eq>&>(),
                                            std::declval<const Complete<Eq>&>(),
                                            Direction::X));

template <typename Eq, typename N = constant<kDefaultChunkSize>>
using VectorizedFluxT = decltype(std::declval<const Eq&>().Flux(
    std::declval<ConservativeArray<Eq, N::value>&>(),
    std::declval<const CompleteArray<Eq, N::value>&>(), Direction::X));

template <typename Eq>
using ScalarReconstructT = decltype(std::declval<const Eq&>().Reconstruct(
    std::declval<Complete<Eq>&>(), std::declval<const Conservative<Eq>&>()));

template <typename Eq, typename N = constant<kDefaultChunkSize>>
using VectorizedReconstructT = decltype(std::declval<const Eq&>().Reconstruct(
    std::declval<CompleteArray<Eq, N::value>&>(),
    std::declval<const ConservativeArray<Eq, N::value>&>()));

template <typename Equation>
struct HasScalarFlux : is_detected<ScalarFluxT, Equation> {};

template <typename Equation, int N = kDefaultChunkSize>
struct HasVectorizedFlux : is_detected<VectorizedFluxT, Equation, constant<N>> {
};

template <typename Equation>
struct HasScalarReconstruction : is_detected<ScalarReconstructT, Equation> {};

template <typename Equation, int N = kDefaultChunkSize>
struct HasVectorizedReconstruction
    : is_detected<VectorizedReconstructT, Equation, constant<N>> {};

template <typename Equation>
struct HasReconstruction : disjunction<HasScalarReconstruction<Equation>,
                                       HasVectorizedReconstruction<Equation>> {
};

template <typename Extents>
constexpr std::array<std::ptrdiff_t, Extents::rank()>
AsArray(Extents e) noexcept {
  std::array<std::ptrdiff_t, Extents::rank()> array{};
  for (std::size_t r = 0; r < Extents::rank(); ++r) {
    array[r] = e.extent(r);
  }
  return array;
}

template <typename Extent>
auto Shrink(const layout_left::mapping<Extent>& layout, Direction dir,
            std::ptrdiff_t n = 1) {
  const std::array<std::ptrdiff_t, Extent::rank()> extents =
      AsArray(layout.extents());
  return layout_left::mapping<Extent>(Extent(Shift(extents, dir, -n)));
}

template <typename State, typename Layout, int Rank>
StridedView<State> Subview(const View<State, Layout, Rank>& state, const IndexBox<Rank>& box) {
  return boost::hana::unpack(state.Members(), [&](const auto&... pdview) {
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
    return StridedView<State>{subview(pdview, box)...};
  });
}

} // namespace fub

#endif
