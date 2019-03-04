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
#include "fub/core/mdspan.hpp"
#include "fub/ForEach.hpp"

#include <boost/hana/tuple.hpp>

#include <utility>

namespace fub {
template <typename... Projections> struct GradientDetector {
  boost::hana::tuple<std::pair<Projections, double>...> conditions;

  GradientDetector(const std::pair<Projections, double>&... conds)
      : conditions{conds...} {}

  template <typename Tags, typename StateView>
  void TagCellsForRefinement(Tags tags, StateView states,
                             const CartesianCoordinates& coords) {
    using Equation = typename StateView::EquationType;
    using Complete = typename Equation::Complete;
    Complete sL;
    Complete sM;
    Complete sR;
    for (std::size_t dir = 0; dir < Extents<0>(states).rank(); ++dir) {
      ForEachIndex(Shrink(Mapping<0>(states), Direction(dir), 2), [&](auto... is) {
        boost::hana::for_each(conditions, [&](auto cond) {
          auto&& [proj, tolerance] = cond;
          std::array<std::ptrdiff_t, sizeof...(is)> index{is...};
          Load(sL, states, index);
          Load(sM, states, Shift(index, Direction(dir), 1));
          Load(sR, states, Shift(index, Direction(dir), 2));
          auto&& xL = std::invoke(proj, sL);
          auto&& xM = std::invoke(proj, sM);
          auto&& xR = std::invoke(proj, sR);
          Eigen::Vector3d dx = coords.dx();
          tags(is...) |= ((std::abs(xM - xL) + std::abs(xR - xM)) /
                          (2 * dx[dir])) > tolerance;
        });
      });
    }
  }
};

template <typename... Ps>
GradientDetector(const std::pair<Ps, double>&...)->GradientDetector<Ps...>;
} // namespace fub

#endif