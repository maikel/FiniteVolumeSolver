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

#ifndef FUB_TAGGING_TAG_BUFFER_HPP
#define FUB_TAGGING_TAG_BUFFER_HPP

#include "fub/CartesianCoordinates.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"
#include "fub/core/mdspan.hpp"

#include <boost/hana/tuple.hpp>

#include <utility>

namespace fub {
struct TagBuffer {
  int buffer_width_;

  explicit TagBuffer(int width) : buffer_width_{width} {}

  template <typename Tags, typename StateView>
  void TagCellsForRefinement(const Tags& tags, const StateView& states,
                             const CartesianCoordinates& /* coords */) {
    for (std::size_t dir = 0; dir < Extents<0>(states).rank(); ++dir) {
      ForEachIndex(Shrink(Box<0>(states), Direction(dir), {2, 2}), [&](auto... is) {
        if (tags(is...) == 1) {
          std::array<std::ptrdiff_t, sizeof...(is)> index{is...};
          for (int width = 1; width <= buffer_width_; ++width) {
            if (index[dir] + width < Extents<0>(states).extent(dir) &&
                !tags(Shift(index, Direction(dir), width))) {
              tags(Shift(index, Direction(dir), width)) = 2;
            }
            if (0 <= index[dir] - width &&
                !tags(Shift(index, Direction(dir), -width))) {
              tags(Shift(index, Direction(dir), -width)) = 2;
            }
          }
        }
      });
    }
  }

  template <typename Tags, typename StateView, typename CutCellData>
  void TagCellsForRefinement(const Tags& tags, const StateView& states,
                             const CutCellData&,
                             const CartesianCoordinates& /* coords */) {
    for (std::size_t dir = 0; dir < Extents<0>(states).rank(); ++dir) {
      ForEachIndex(Shrink(Box<0>(states), Direction(dir), {2, 2}), [&](auto... is) {
        if (tags(is...) == 1) {
          std::array<std::ptrdiff_t, sizeof...(is)> index{is...};
          for (int width = 1; width <= buffer_width_; ++width) {
            if (index[dir] + width < Extents<0>(states).extent(dir) &&
                !tags(Shift(index, Direction(dir), width))) {
              tags(Shift(index, Direction(dir), width)) = 2;
            }
            if (0 <= index[dir] - width &&
                !tags(Shift(index, Direction(dir), -width))) {
              tags(Shift(index, Direction(dir), -width)) = 2;
            }
          }
        }
      });
    }
  }
};

} // namespace fub

#endif
