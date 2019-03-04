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

template <typename Eq, typename N = constant<kChunkSize>>
using VectorizedFluxT = decltype(std::declval<const Eq&>().Flux(
    std::declval<ConsArray<Eq, N::value>&>(),
    std::declval<const CompleteArray<Eq, N::value>&>(), Direction::X));

template <typename Eq>
using ScalarReconstructT = decltype(std::declval<const Eq&>().Reconstruct(
    std::declval<Complete<Eq>&>(), std::declval<const Conservative<Eq>&>()));

template <typename Eq, typename N = constant<kChunkSize>>
using VectorizedReconstructT = decltype(std::declval<const Eq&>().Reconstruct(
    std::declval<CompleteArray<Eq, N::value>&>(),
    std::declval<const ConsArray<Eq, N::value>&>()));

template <typename Equation>
struct HasScalarFlux : is_detected<ScalarFluxT, Equation> {};

template <typename Equation, int N = kChunkSize>
struct HasVectorizedFlux : is_detected<VectorizedFluxT, Equation, constant<N>> {
};

template <typename Equation>
struct HasScalarReconstruction : is_detected<ScalarReconstructT, Equation> {};

template <typename Equation, int N = kChunkSize>
struct HasVectorizedReconstruction
    : is_detected<VectorizedReconstructT, Equation, constant<N>> {};

template <typename Equation>
struct HasReconstruction : disjunction<HasScalarReconstruction<Equation>,
                                       HasVectorizedReconstruction<Equation>> {
};

template <std::size_t N>
std::array<std::ptrdiff_t, N> Shift(const std::array<std::ptrdiff_t, N>& idx,
                                    Direction dir, std::ptrdiff_t shift) {
  auto shifted(idx);
  shifted[int(dir)] += shift;
  return shifted;
}

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

template <typename T, typename E, typename L, typename A>
strided_mdspan<T, E, A> ViewInnerRegion(basic_mdspan<T, E, L, A> mdspan,
                                        Direction dir, int gcw) {
  constexpr int Rank = E::rank();
  std::array<std::pair<std::ptrdiff_t, std::ptrdiff_t>, Rank> slices;
  for (int r = 0; r < Rank; ++r) {
    slices[r] = std::pair{0L, mdspan.extent(r)};
  }
  const int dir_v = int(dir);
  FUB_ASSERT(dir_v < Rank);
  slices[dir_v].first += gcw;
  slices[dir_v].second -= gcw;
  return std::apply([&](auto... slices) { return subspan(mdspan, slices...); },
                    slices);
}

template <typename State>
StridedView<State> ViewInnerRegion(View<State> state, Direction dir, int gcw) {
  return boost::hana::unpack(state.Members(), [&](auto... mdspan) {
    return StridedView<State>{ViewInnerRegion(mdspan, dir, gcw)...};
  });
}

} // namespace fub

#endif