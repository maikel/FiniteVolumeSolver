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
    decltype(std::declval<const Eq&>().Flux(std::declval<Cons<Eq>&>(),
                                            std::declval<const Complete<Eq>&>(),
                                            Direction::X));

template <typename Eq, typename N = constant<kChunkSize>>
using VectorizedFluxT = decltype(std::declval<const Eq&>().Flux(
    std::declval<ConsArray<Eq, N::value>&>(),
    std::declval<const CompleteArray<Eq, N::value>&>(), Direction::X));

template <typename Eq>
using ScalarReconstructT = decltype(std::declval<const Eq&>().Reconstruct(
    std::declval<Complete<Eq>&>(), std::declval<const Cons<Eq>&>()));

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

template <typename E, typename Function>
Function ForEachIndex(const layout_left::mapping<E>& mapping,
                      Function function) {
  if constexpr (E::rank() == 1) {
    for (int i = 0; i < mapping.extents().extent(0); ++i) {
      function(i);
    }
  } else if constexpr (E::rank() == 2) {
    for (int i = 0; i < mapping.extents().extent(1); ++i) {
      for (int j = 0; j < mapping.extents().extent(0); ++j) {
        function(j, i);
      }
    }
  } else if constexpr (E::rank() == 3) {
    for (int i = 0; i < mapping.extents().extent(2); ++i) {
      for (int j = 0; j < mapping.extents().extent(1); ++j) {
        for (int k = 0; k < mapping.extents().extent(0); ++k) {
          function(k, j, i);
        }
      }
    }
  }
  return function;
}

template <typename E, typename Function>
Function ForEachIndex(const layout_stride::mapping<E>& mapping,
                      Function function) {
  if constexpr (E::rank() == 1) {
    for (int i = 0; i < mapping.extents().extent(0); ++i) {
      function(i);
    }
  } else if constexpr (E::rank() == 2) {
    for (int i = 0; i < mapping.extents().extent(1); ++i) {
      for (int j = 0; j < mapping.extents().extent(0); ++j) {
        function(j, i);
      }
    }
  } else if constexpr (E::rank() == 3) {
    for (int i = 0; i < mapping.extents().extent(2); ++i) {
      for (int j = 0; j < mapping.extents().extent(1); ++j) {
        for (int k = 0; k < mapping.extents().extent(0); ++k) {
          function(k, j, i);
        }
      }
    }
  }
  return function;
}

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
  std::array<std::ptrdiff_t, Extents::rank()> array;
  for (int r = 0; r < Extents::rank(); ++r) {
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

template <typename Equation, typename CompleteView, typename ConsView>
void ReconstructStatesFromCons(const Equation& eq, CompleteView sview,
                               ConsView cview) {
  FUB_ASSERT(Extents(sview) == Extents(cview));
  Complete<Equation> complete;
  Cons<Equation> cons;
  ForEachIndex(Mapping(sview), [&](auto... is) {
    Load(cons, cview, {is...});
    if constexpr (HasReconstruction<Equation>()) {
      eq.Reconstruct(complete, cons);
    } else if (sizeof(Complete<Equation>) == sizeof(Cons<Equation>)) {
      ::std::memcpy(&complete, &cons, sizeof(Cons<Equation>));
    }
    Store(sview, complete, {is...});
  });
}

} // namespace fub

#endif