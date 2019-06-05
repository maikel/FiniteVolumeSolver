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

#ifndef FUB_TAGGING_GRADIENT_DETECTOR_HPP
#define FUB_TAGGING_GRADIENT_DETECTOR_HPP

#include "fub/ForEach.hpp"
#include "fub/StateArray.hpp"
#include "fub/StateRow.hpp"
#include "fub/core/mdspan.hpp"
#include "fub/ext/Eigen.hpp"

#include <utility>

namespace fub {

template <std::size_t DestRank, std::size_t SrcRank>
std::array<std::ptrdiff_t, DestRank>
EmbedInSpace(const std::array<std::ptrdiff_t, SrcRank>& index) {
  if constexpr (DestRank == SrcRank) {
    return index;
  } else {
    std::array<std::ptrdiff_t, DestRank> dest{};
    std::copy_n(index.begin(), static_cast<std::ptrdiff_t>(SrcRank),
                dest.begin());
    return dest;
  }
}

template <typename Proj, typename State>
using Projection_t =
    decltype(std::invoke(std::declval<Proj>(), std::declval<State>()));

template <typename Proj, typename State>
using IsProjection = is_detected<Projection_t, Proj, State>;

template <typename Equation, typename... Projections>
class GradientDetectorBase {
public:
  GradientDetectorBase(const Equation& equation,
                       const std::pair<Projections, double>&... conds)
      : equation_{equation}, conditions_{conds...} {}

  Equation& GetEquation() noexcept { return equation_; }
  const Equation& GetEquation() const noexcept { return equation_; }

  const std::tuple<std::pair<Projections, double>...>& GetConditions() const
      noexcept {
    return conditions_;
  }

protected:
  Equation equation_;
  std::tuple<std::pair<Projections, double>...> conditions_;
};

template <typename Equation, typename... Projections>
class ScalarGradientDetector
    : public GradientDetectorBase<Equation, Projections...> {
public:
  ScalarGradientDetector(const Equation& equation,
                         const std::pair<Projections, double>&... conds);

  template <int TagRank>
  void
  TagCellsForRefinement(const PatchDataView<char, TagRank, layout_stride>& tags,
                        const View<const Complete<Equation>>& states);

private:
  using Base = GradientDetectorBase<Equation, Projections...>;
  Complete<Equation> sL{Base::equation_};
  Complete<Equation> sM{Base::equation_};
  Complete<Equation> sR{Base::equation_};
};

template <typename Equation, typename... Projections>
class ArrayGradientDetector
    : public GradientDetectorBase<Equation, Projections...> {
public:
  static constexpr int Rank = Equation::Rank();

  ArrayGradientDetector(const Equation& equation,
                        const std::pair<Projections, double>&... conds);

  void
  TagCellsForRefinement(const PatchDataView<char, Rank, layout_stride>& tags,
                        const View<const Complete<Equation>>& states);

private:
  using Base = GradientDetectorBase<Equation, Projections...>;
  CompleteArray<Equation> sL_{Base::equation_};
  CompleteArray<Equation> sM_{Base::equation_};
  CompleteArray<Equation> sR_{Base::equation_};
};

template <bool IsArray, typename Eq, typename... Ps>
struct GradientDetectorImpl {
  using type = ArrayGradientDetector<Eq, Ps...>;
};

template <typename Eq, typename... Ps>
struct GradientDetectorImpl<false, Eq, Ps...> {
  using type = ScalarGradientDetector<Eq, Ps...>;
};

template <typename Eq, typename... Ps>
struct GradientDetector
    : public GradientDetectorImpl<
          (IsProjection<Ps&, const CompleteArray<Eq>&>::value && ...), Eq,
          Ps...>::type {
public:
  static constexpr bool is_array_based =
      (IsProjection<Ps&, const CompleteArray<Eq>&>::value && ...);
  using Base = typename GradientDetectorImpl<is_array_based, Eq, Ps...>::type;

  GradientDetector(const Eq& equation, const std::pair<Ps, double>&... projs)
      : Base(equation, projs...) {}
};

template <typename Eq, typename... Ps>
GradientDetector(const Eq&, const std::pair<Ps, double>&...)
    ->GradientDetector<Eq, Ps...>;

template <typename Equation, typename... Projections>
ScalarGradientDetector<Equation, Projections...>::ScalarGradientDetector(
    const Equation& equation, const std::pair<Projections, double>&... conds)
    : Base(equation, conds...) {}

template <typename Equation, typename... Projections>
template <int TagRank>
void ScalarGradientDetector<Equation, Projections...>::TagCellsForRefinement(
    const PatchDataView<char, TagRank, layout_stride>& tags,
    const View<const Complete<Equation>>& states) {
  constexpr std::size_t sTagRank = static_cast<std::size_t>(TagRank);
  const auto tagbox = tags.Box();
  for (std::size_t dir = 0; dir < Extents<0>(states).rank(); ++dir) {
    ForEachIndex(
        Shrink(Box<0>(states), Direction(dir), {1, 1}), [&](auto... is) {
          using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
          const Index mid{is...};
          if (Contains(tagbox, EmbedInSpace<sTagRank>(mid))) {
            const Index left = Shift(mid, Direction(dir), -1);
            const Index right = Shift(mid, Direction(dir), 1);
            boost::mp11::tuple_for_each(Base::conditions_, [&](auto cond) {
              auto&& [proj, tolerance] = cond;
              Load(sL, states, left);
              Load(sM, states, mid);
              Load(sR, states, right);
              auto&& xL = std::invoke(proj, sL);
              auto&& xM = std::invoke(proj, sM);
              auto&& xR = std::invoke(proj, sR);
              if (xM != xR || xM != xL) {
                const double left =
                    std::abs(xM - xL) / (std::abs(xM) + std::abs(xL));
                const double right =
                    std::abs(xM - xR) / (std::abs(xM) + std::abs(xR));
                tags(EmbedInSpace<sTagRank>(mid)) |=
                    static_cast<char>(left > tolerance || right > tolerance);
              }
            });
          }
        });
  }
}

template <typename Equation, typename... Projections>
ArrayGradientDetector<Equation, Projections...>::ArrayGradientDetector(
    const Equation& equation, const std::pair<Projections, double>&... conds)
    : Base(equation, conds...) {}

template <typename Equation, typename... Projections>
void ArrayGradientDetector<Equation, Projections...>::TagCellsForRefinement(
    const PatchDataView<char, Rank, layout_stride>& tags,
    const View<const Complete<Equation>>& states) {
  Direction dir = Direction::X;
  BasicView shifted = Subview(states, Grow(tags.Box(), dir, {1, 1}));
  std::tuple views{tags, Shrink(shifted, dir, {0, 2}),
                   Shrink(shifted, dir, {1, 1}), Shrink(shifted, dir, {2, 0})};
  FUB_ASSERT(Box<0>(std::get<2>(views)) == tags.Box());
  ForEachRow(views, [&](span<char> tags_row,
                        const Row<const Complete<Equation>>& rL,
                        const Row<const Complete<Equation>>& rM,
                        const Row<const Complete<Equation>>& rR) {
    char* tags_pointer = tags_row.begin();
    ViewPointer<const Complete<Equation>> pL = Begin(rL);
    ViewPointer<const Complete<Equation>> pM = Begin(rM);
    ViewPointer<const Complete<Equation>> pR = Begin(rR);
    ViewPointer<const Complete<Equation>> pEnd = End(rR);
    int n = static_cast<int>(get<0>(pEnd) - get<0>(pR));
    if (n >= kDefaultChunkSize) {
      Load(sM_, pL);
      Load(sR_, pM);
    }
    const Array<char, 1> ones = Array<char, 1>::Constant(char(1));
    while (n >= kDefaultChunkSize) {
      sL_ = sM_;
      sM_ = sR_;
      Load(sR_, pR);
      auto tags_array = Array<char, 1>::Map(tags_pointer);
      boost::mp11::tuple_for_each(Base::conditions_, [&](auto&& condition) {
        auto&& [proj, tolerance] = condition;
        const Array1d xL = std::invoke(proj, sL_);
        const Array1d xM = std::invoke(proj, sM_);
        const Array1d xR = std::invoke(proj, sR_);
        const Array1d left = (xM - xL).abs() / (xM.abs() + xL.abs());
        const Array1d right = (xM - xR).abs() / (xM.abs() + xR.abs());
        tags_array = (left > tolerance).select(tags_array, ones);
        tags_array = (right > tolerance).select(tags_array, ones);
        Advance(pR, kDefaultChunkSize);
        tags_pointer += kDefaultChunkSize;
        n = static_cast<int>(get<0>(pEnd) - get<0>(pR));
      });
    }
    LoadN(sL_, pL, n);
    LoadN(sM_, pM, n);
    LoadN(sR_, pR, n);
    Array<char, 1> tags_array = Array<char, 1>::Zero();
    for (int i = 0; i < n; ++i) {
      tags_array[i] = tags_pointer[i];
    }
    boost::mp11::tuple_for_each(Base::conditions_, [&](auto&& condition) {
      auto&& [proj, tolerance] = condition;
      const Array1d xL = std::invoke(proj, sL_);
      const Array1d xM = std::invoke(proj, sM_);
      const Array1d xR = std::invoke(proj, sR_);
      const Array1d left = (xM - xL).abs() / (xM.abs() + xL.abs());
      const Array1d right = (xM - xR).abs() / (xM.abs() + xR.abs());
      tags_array = (left > tolerance).select(tags_array, ones);
      tags_array = (right > tolerance).select(tags_array, ones);
    });
    for (int i = 0; i < n; ++i) {
      tags_pointer[i] = tags_array[i];
    }
  });
}

} // namespace fub

#endif
