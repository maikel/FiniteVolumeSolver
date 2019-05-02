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

#include "fub/CartesianCoordinates.hpp"
#include "fub/ForEach.hpp"
#include "fub/core/mdspan.hpp"

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

template <typename Equation, typename... Projections> struct GradientDetector {
  Equation equation_;
  std::tuple<std::pair<Projections, double>...> conditions_;
  Complete<Equation> sL{equation_};
  Complete<Equation> sM{equation_};
  Complete<Equation> sR{equation_};

  GradientDetector(const Equation& equation,
                   const std::pair<Projections, double>&... conds)
      : equation_{equation}, conditions_{conds...} {}

  template <int TagRank>
  void TagCellsForRefinement(const PatchDataView<char, TagRank>& tags,
                             const BasicView<const Complete<Equation>>& states,
                             const CartesianCoordinates&) {
    constexpr std::size_t sTagRank = static_cast<std::size_t>(TagRank);
    const auto tagbox = tags.Box();
    for (std::size_t dir = 0; dir < Extents<0>(states).rank(); ++dir) {
      ForEachIndex(Shrink(Box<0>(states), Direction(dir), {1, 1}), [&](auto... is) {
        using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
        const Index mid{is...};
        if (Contains(tagbox, EmbedInSpace<sTagRank>(mid))) {
            const Index left = Shift(mid, Direction(dir), -1);
            const Index right = Shift(mid, Direction(dir), 1);
            boost::mp11::tuple_for_each(conditions_, [&](auto cond) {
              auto&& [proj, tolerance] = cond;
              Load(sL, states, left);
              Load(sM, states, mid);
              Load(sR, states, right);
              auto&& xL = std::invoke(proj, sL);
              auto&& xM = std::invoke(proj, sM);
              auto&& xR = std::invoke(proj, sR);
              //  Eigen::Vector3d dx = coords.dx();
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

  template <int TagRank, typename CutCellData>
  void TagCellsForRefinement(const PatchDataView<char, TagRank>& tags,
                             const BasicView<const Complete<Equation>>& states,
                             const CutCellData& cutcell_data,
                             const CartesianCoordinates&) {
    const auto& flags = cutcell_data.flags;
    FUB_ASSERT(Contains(flags.Box(), Box<0>(states)));
    const auto tagbox = tags.Box();
    for (std::size_t dir = 0; dir < Extents<0>(states).rank(); ++dir) {
      ForEachIndex(
          Shrink(Box<0>(states), Direction(dir), {0, 2}), [&](auto... is) {
            using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
            if (Contains(tagbox, Index{is...})) {
            std::array<std::ptrdiff_t, sizeof...(is)> left{is...};
            std::array<std::ptrdiff_t, sizeof...(is)> mid =
                Shift(left, Direction(dir), 1);
            std::array<std::ptrdiff_t, sizeof...(is)> right =
                Shift(mid, Direction(dir), 1);
            if (flags(mid).isRegular()) {
              boost::mp11::tuple_for_each(conditions_, [&](auto cond) {
                auto&& [proj, tolerance] = cond;
                Load(sL, states, left);
                Load(sM, states, mid);
                Load(sR, states, right);
                auto&& xL = std::invoke(proj, sL);
                auto&& xM = std::invoke(proj, sM);
                auto&& xR = std::invoke(proj, sR);
                // Eigen::Vector3d dx = coords.dx();
                if (xM != xR || xM != xL) {
                  const double left =
                      std::abs(xM - xL) / (std::abs(xM) + std::abs(xL));
                  const double right =
                      std::abs(xM - xR) / (std::abs(xM) + std::abs(xR));
                  tags(mid) |=
                      static_cast<char>(left > tolerance || right > tolerance);
                }
              });
            }
            }
          });
    }
  }
};

template <typename Eq, typename... Ps>
GradientDetector(const Eq&, const std::pair<Ps, double>&...)
    ->GradientDetector<Eq, Ps...>;
} // namespace fub

#endif
