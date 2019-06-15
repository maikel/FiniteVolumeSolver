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

#ifndef FUB_STATE_ROW_HPP
#define FUB_STATE_ROW_HPP

#include "fub/State.hpp"

namespace fub {

template <typename T, typename Depth> struct DepthToRowType;

template <typename T> struct DepthToRowType<T, ScalarDepth> {
  using type = span<T>;
};

template <typename T, int Rank> struct DepthToRowType<T, VectorDepth<Rank>> {
  using type = mdspan<T, 2, layout_stride>;
};

template <typename State> struct RowBaseImpl {
  template <typename Depth>
  using fn = typename DepthToRowType<double, Depth>::type;
  using type = boost::mp11::mp_transform<fn, typename State::Depths>;
};

template <typename State> struct RowBaseImpl<const State> {
  template <typename Depth>
  using fn = typename DepthToRowType<const double, Depth>::type;
  using type = boost::mp11::mp_transform<fn, typename State::Depths>;
};

template <typename State> using RowBase = typename RowBaseImpl<State>::type;

template <typename State> struct Row : RowBase<State> {
  static constexpr int Rank = 1;
  using Traits = StateTraits<RowBase<State>>;
  using ElementType =
      std::conditional_t<std::is_const_v<State>, const double, double>;

  Row() = default;

  Row(const ViewPointer<State>& pointer, std::ptrdiff_t extent) {
    ForEachVariable(
        overloaded{[extent](span<ElementType>& s, auto p) {
                     s = span<ElementType>(p, extent);
                   },
                   [extent](mdspan<ElementType, 2, layout_stride>& mds,
                            const auto& p) {
                     layout_stride::mapping<dynamic_extents<2>> mapping{
                         dynamic_extents<2>{extent, 2}, {1, p.second}};
                     mds = mdspan<ElementType, 2, layout_stride>(p.first,
                                                                 mapping);
                   }},
        *this, pointer);
  }
};

template <typename State>
struct StateTraits<Row<State>> : StateTraits<RowBase<State>> {};

template <typename State> ViewPointer<State> Begin(const Row<State>& row) {
  ViewPointer<State> pointer;
  ForEachVariable(
      overloaded{
          [&](double*& p, span<double> s) { p = s.begin(); },
          [&](const double*& p, span<const double> s) { p = s.begin(); },
          [&](auto& p, const mdspan<double, 2, layout_stride>& mds) {
            p = std::pair{mds.get_span().begin(), mds.stride(1)};
          },
          [&](auto& p, const mdspan<const double, 2, layout_stride>& mds) {
            p = std::pair{mds.get_span().begin(), mds.stride(1)};
          }},
      pointer, row);
  return pointer;
}

template <typename State> ViewPointer<State> End(const Row<State>& row) {
  ViewPointer<State> pointer;
  ForEachVariable(
      overloaded{
          [&](double*& p, span<double> s) { p = s.end(); },
          [&](const double*& p, span<const double> s) { p = s.end(); },
          [&](auto& p, const mdspan<double, 2, layout_stride>& mds) {
            p = std::pair{mds.get_span().end(), mds.stride(1)};
          },
          [&](auto& p, const mdspan<const double, 2, layout_stride>& mds) {
            p = std::pair{mds.get_span().end(), mds.stride(1)};
          }},
      pointer, row);
  return pointer;
}

template <typename Tuple, typename Function>
auto Transform(Tuple&& tuple, Function f) {
  return std::apply([&](auto&&... xs) { return std::tuple{f(xs)...}; }, tuple);
}

template <Direction Dir> struct ToStride {
  template <typename T, int R, typename L>
  std::ptrdiff_t operator()(const PatchDataView<T, R, L>& pdv) const {
    if constexpr (R == 1) {
      return pdv.Extent(0);
    } else {
      return pdv.Stride(static_cast<int>(Dir));
    }
  }

  template <typename T, typename L, int R>
  std::ptrdiff_t operator()(const BasicView<T, L, R>& view) const {
    return get<0>(view).Stride(static_cast<int>(Dir));
  }
};

template <Direction Dir, typename T, int R, typename L>
std::ptrdiff_t Extent(const PatchDataView<T, R, L>& pdv) {
  return pdv.Extent(static_cast<std::size_t>(Dir));
}

template <Direction Dir, typename T, typename L, int R>
std::ptrdiff_t Extent(const BasicView<T, L, R>& view) {
  return Extent<Dir>(get<0>(view));
}

struct ToRow {
  std::ptrdiff_t extent;

  template <typename T> span<T> operator()(T* pointer) {
    return {pointer, extent};
  }

  template <typename T> Row<T> operator()(const ViewPointer<T>& pointer) {
    return {pointer, extent};
  }
};

template <typename T> void Advance(T*& pointer, std::ptrdiff_t n) {
  pointer += n;
}

template <typename T, int R, typename L>
T* Begin(const PatchDataView<T, R, L>& pdv) {
  return pdv.Span().begin();
}

template <typename T, int R, typename L>
T* End(const PatchDataView<T, R, L>& pdv) {
  return pdv.Span().end();
}

template <int N, typename T> constexpr T* GetOrForward(T* pointer) noexcept {
  return pointer;
}

template <int N, typename T>
constexpr decltype(auto) GetOrForward(const ViewPointer<T>& pointer) noexcept {
  return get<N>(pointer);
}

template <typename Tuple, typename Function>
void ForEachRow(const Tuple& views, Function f) {
  std::tuple firsts = Transform(views, [](const auto& v) { return Begin(v); });
  std::tuple lasts = Transform(views, [](const auto& v) { return End(v); });
  std::tuple strides = Transform(views, ToStride<Direction::Y>());
  const std::ptrdiff_t row_extent = Extent<Direction::X>(std::get<0>(views));
  const auto& first = std::get<0>(firsts);
  const auto& last = std::get<0>(lasts);
  while (GetOrForward<0>(last) - GetOrForward<0>(first) >=
         std::get<0>(strides)) {
    std::tuple rows = Transform(firsts, ToRow{row_extent});
    std::apply(f, rows);
    std::tuple firsts_and_strides = Zip(firsts, strides);
    std::apply(
        [](auto&... ps) { (Advance(std::get<0>(ps), std::get<1>(ps)), ...); },
        firsts_and_strides);
    firsts =
        std::apply([](auto&&... fs) { return std::tuple{std::get<0>(fs)...}; },
                   firsts_and_strides);
  }
}

} // namespace fub

#endif
