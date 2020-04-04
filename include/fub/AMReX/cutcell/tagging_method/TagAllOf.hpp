// Copyright (c) 2020 Maikel Nadolski
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

#ifndef FUB_AMREX_CUTCELL_TAGGING_METHOD_TAG_ALL_OF_HPP
#define FUB_AMREX_CUTCELL_TAGGING_METHOD_TAG_ALL_OF_HPP

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/Duration.hpp"

#include <tuple>

namespace fub::amrex::cutcell {

template <typename... Tagging> struct TagAllOf {
  TagAllOf(Tagging... tagging) : tagging_{std::move(tagging)...} {}

  void TagCellsForRefinement(::amrex::TagBoxArray& tags,
                             GriddingAlgorithm& gridding, int level,
                             Duration time_point) {
    std::apply(
        [&tags, time_point, level, &gridding](Tagging&... tagging) {
          (tagging.TagCellsForRefinement(tags, gridding, level, time_point),
           ...);
        },
        tagging_);
  }

  std::tuple<Tagging...> tagging_;
};

} // namespace fub::amrex::cutcell

#endif