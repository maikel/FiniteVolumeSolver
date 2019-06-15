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

#ifndef FUB_PATCH_DATA_VIEW_HPP
#define FUB_PATCH_DATA_VIEW_HPP

#include "fub/Direction.hpp"
#include "fub/core/mdspan.hpp"

#include <array>
#include <functional>
#include <tuple>
#include <utility>

namespace fub {
template <int Rank>
using Index = std::array<std::ptrdiff_t, static_cast<std::size_t>(Rank)>;

template <int Rank> struct IndexBox {
  static_assert(Rank >= 0);
  Index<Rank> lower;
  Index<Rank> upper;
};

template <int Rank>
bool operator==(const IndexBox<Rank>& b1, const IndexBox<Rank>& b2) {
  return b1.lower == b2.lower && b1.upper == b2.upper;
}

template <int Rank>
bool operator!=(const IndexBox<Rank>& b1, const IndexBox<Rank>& b2) {
  return !(b1 == b2);
}

template <int Rank>
bool Contains(const IndexBox<Rank>& b1, const IndexBox<Rank>& b2) {
  return std::equal(b1.lower.begin(), b1.lower.end(), b2.lower.begin(),
                    std::less_equal<>{}) &&
         std::equal(b2.upper.begin(), b2.upper.end(), b1.upper.begin(),
                    std::less_equal<>{});
}

template <int Rank>
bool Contains(const IndexBox<Rank>& box, const Index<Rank>& index) {
  Index<Rank> upper;
  std::transform(index.begin(), index.end(), upper.begin(),
                 [](std::ptrdiff_t i) { return i + 1; });
  return Contains(box, IndexBox<Rank>{index, upper});
}

template <int Rank>
IndexBox<Rank> Intersect(const IndexBox<Rank>& b1, const IndexBox<Rank>& b2) {
  constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  std::array<std::ptrdiff_t, sRank> lower;
  std::transform(
      b1.lower.begin(), b1.lower.end(), b2.lower.begin(), lower.begin(),
      [](std::ptrdiff_t l1, std::ptrdiff_t l2) { return std::max(l1, l2); });
  std::array<std::ptrdiff_t, sRank> upper;
  std::transform(
      b1.upper.begin(), b1.upper.end(), b2.upper.begin(), upper.begin(),
      [](std::ptrdiff_t u1, std::ptrdiff_t u2) { return std::min(u1, u2); });
  std::transform(
      lower.begin(), lower.end(), upper.begin(), upper.begin(),
      [](std::ptrdiff_t l, std::ptrdiff_t u) { return std::max(l, u); });
  return IndexBox<Rank>{lower, upper};
}

template <int Rank>
IndexBox<Rank> Grow(const IndexBox<Rank>& box, Direction dir,
                    const std::array<std::ptrdiff_t, 2>& shifts) {
  constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  std::array<std::ptrdiff_t, sRank> lower = box.lower;
  lower[static_cast<std::size_t>(dir)] -= shifts[0];
  std::array<std::ptrdiff_t, sRank> upper = box.upper;
  upper[static_cast<std::size_t>(dir)] += shifts[1];
  return IndexBox<Rank>{lower, upper};
}

template <int Rank>
IndexBox<Rank> Shrink(const IndexBox<Rank>& box, Direction dir,
                      const std::array<std::ptrdiff_t, 2>& shifts) {
  constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  std::array<std::ptrdiff_t, sRank> lower = box.lower;
  lower[static_cast<std::size_t>(dir)] += shifts[0];
  std::array<std::ptrdiff_t, sRank> upper = box.upper;
  upper[static_cast<std::size_t>(dir)] -= shifts[1];
  return IndexBox<Rank>{lower, upper};
}

template <int Rank>
IndexBox<Rank> Embed(const IndexBox<Rank - 1>& box,
                     const std::array<std::ptrdiff_t, 2>& limits) {
  constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  std::array<std::ptrdiff_t, sRank> lower;
  std::copy_n(box.lower.begin(), Rank - 1, lower.begin());
  lower[Rank - 1] = limits[0];
  std::array<std::ptrdiff_t, sRank> upper;
  std::copy_n(box.upper.begin(), Rank - 1, upper.begin());
  upper[Rank - 1] = limits[1];
  return IndexBox<Rank>{lower, upper};
}

template <int Rank, int OtherRank>
IndexBox<Rank> Project(const IndexBox<OtherRank>& box) {
  constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  std::array<std::ptrdiff_t, sRank> lower;
  std::copy_n(box.lower.begin(), Rank, lower.begin());
  std::array<std::ptrdiff_t, sRank> upper;
  std::copy_n(box.upper.begin(), Rank, upper.begin());
  return IndexBox<Rank>{lower, upper};
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

template <typename T, int Rank, typename Layout = layout_left>
struct PatchDataView;

template <typename T, int Rank, typename Layout> struct PatchDataViewBase {
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);

  PatchDataViewBase() = default;
  PatchDataViewBase(const PatchDataViewBase&) = default;
  PatchDataViewBase& operator=(const PatchDataViewBase&) = default;

  PatchDataViewBase(const mdspan<T, sRank, Layout>& data,
                    const std::array<std::ptrdiff_t, sRank>& origin)
      : mdspan_{data}, origin_{origin} {}

  mdspan<T, sRank, Layout> mdspan_;
  std::array<std::ptrdiff_t, sRank> origin_;
};

template <typename T, int Rank, typename Layout>
struct PatchDataViewBase<const T, Rank, Layout> {
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);

  PatchDataViewBase() = default;

  PatchDataViewBase(const PatchDataViewBase&) = default;
  PatchDataViewBase& operator=(const PatchDataViewBase&) = default;

  PatchDataViewBase(const PatchDataView<T, Rank, Layout>& other)
      : mdspan_{other.MdSpan()}, origin_{other.Origin()} {}

  PatchDataViewBase(const mdspan<const T, sRank, Layout>& data,
                    const std::array<std::ptrdiff_t, sRank>& origin)
      : mdspan_{data}, origin_{origin} {}

  mdspan<const T, sRank, Layout> mdspan_;
  std::array<std::ptrdiff_t, sRank> origin_;
};

template <typename T, int R, typename Layout>
struct PatchDataView : public PatchDataViewBase<T, R, Layout> {
  static constexpr std::size_t sRank = static_cast<std::size_t>(R);

  using PatchDataViewBase<T, R, Layout>::PatchDataViewBase;

  const mdspan<T, sRank, Layout>& MdSpan() const noexcept {
    return this->mdspan_;
  }

  const std::array<std::ptrdiff_t, sRank>& Origin() const noexcept {
    return this->origin_;
  }

  dynamic_extents<sRank> Extents() const noexcept {
    return this->mdspan_.extents();
  }

  typename Layout::template mapping<dynamic_extents<sRank>> Mapping() const
      noexcept {
    return this->mdspan_.mapping();
  }

  static constexpr int Rank() noexcept { return R; }

  span<T> Span() const noexcept { return this->mdspan_.get_span(); }

  std::ptrdiff_t Stride(std::size_t n) const { return this->mdspan_.stride(n); }

  std::ptrdiff_t Extent(std::size_t n) const { return this->mdspan_.extent(n); }

  IndexBox<R> Box() const noexcept {
    IndexBox<R> box;
    box.lower = this->origin_;
    std::array<std::ptrdiff_t, sRank> extents = AsArray(Extents());
    std::transform(extents.begin(), extents.end(), box.lower.begin(),
                   box.upper.begin(),
                   [](std::ptrdiff_t extents, std::ptrdiff_t lower) {
                     FUB_ASSERT(extents >= 0);
                     return lower + extents;
                   });
    return box;
  }

  PatchDataView<T, R, layout_stride> Subview(const IndexBox<R>& box) const {
    std::array<std::ptrdiff_t, sRank> offset;
    std::transform(box.lower.begin(), box.lower.end(), this->origin_.begin(),
                   offset.begin(),
                   [](std::ptrdiff_t lower, std::ptrdiff_t origin) {
                     FUB_ASSERT(origin <= lower);
                     return lower - origin;
                   });
    std::array<std::ptrdiff_t, sRank> extents;
    std::transform(box.lower.begin(), box.lower.end(), box.upper.begin(),
                   extents.begin(),
                   [](std::ptrdiff_t lower, std::ptrdiff_t upper) {
                     FUB_ASSERT(lower <= upper);
                     return upper - lower;
                   });
    std::array<std::pair<std::ptrdiff_t, std::ptrdiff_t>, sRank> slice_array;
    std::transform(offset.begin(), offset.end(), extents.begin(),
                   slice_array.begin(),
                   [](std::ptrdiff_t offset, std::ptrdiff_t extents) {
                     return std::make_pair(offset, offset + extents);
                   });
    return std::apply(
        [&](const auto&... slices) {
          return PatchDataView<T, R, layout_stride>(
              subspan(MdSpan(), slices...), box.lower);
        },
        slice_array);
  }

  template <typename... IndexType,
            typename = std::enable_if_t<conjunction<
                std::is_convertible<IndexType, std::ptrdiff_t>...>::value>>
  auto& operator()(IndexType... indices) const {
    std::array<std::ptrdiff_t, sRank> index{
        static_cast<std::ptrdiff_t>(indices)...};
    return this->operator()(index);
  }

  template <typename... IndexType,
            typename = std::enable_if_t<conjunction<
                std::is_convertible<IndexType, std::ptrdiff_t>...>::value>>
  auto& operator()(const std::array<std::ptrdiff_t, sRank>& indices) const {
    std::array<std::ptrdiff_t, sRank> local_index;
    std::transform(indices.begin(), indices.end(), Origin().begin(),
                   local_index.begin(),
                   [](std::ptrdiff_t i, std::ptrdiff_t o) { return i - o; });
    return this->mdspan_(local_index);
  }
};

template <typename T, int Rank, typename Layout>
PatchDataView<T, Rank - 1, Layout>
SliceLast(const PatchDataView<T, Rank, Layout>& pdv, int component = 0) {
  constexpr std::size_t sRank = static_cast<std::size_t>(Rank - 1);
  std::array<std::ptrdiff_t, sRank + 1> index{};
  index[sRank] = component;
  std::array<std::ptrdiff_t, sRank> extents;
  for (std::size_t r = 0; r < sRank; ++r) {
    extents[r] = pdv.Extent(r);
  }
  std::array<std::ptrdiff_t, sRank> strides;
  for (std::size_t r = 0; r < sRank; ++r) {
    strides[r] = pdv.Stride(r);
  }
  std::array<std::ptrdiff_t, sRank> origin;
  std::copy_n(pdv.Origin().begin(), Rank - 1, origin.begin());
  if constexpr (std::is_same_v<Layout, layout_stride>) {
    layout_stride::mapping<dynamic_extents<sRank>> mapping{
        dynamic_extents<sRank>(extents), strides};
    mdspan<T, sRank, Layout> mds(&pdv.MdSpan()(index), mapping);
    return PatchDataView<T, Rank - 1, Layout>(mds, origin);
  } else {
    mdspan<T, sRank, Layout> mds(&pdv.MdSpan()(index), extents);
    return PatchDataView<T, Rank - 1, Layout>(mds, origin);
  }
}

} // namespace fub

#endif
