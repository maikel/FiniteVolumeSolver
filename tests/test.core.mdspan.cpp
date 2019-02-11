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

#include "fub/core/mdspan.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

TEST_CASE("Construct extents") {
  SECTION("static") {
    fub::extents<1, 2> e;
    REQUIRE(e.static_extent(0) == 1);
    REQUIRE(e.static_extent(1) == 2);
    REQUIRE(e.extent(0) == 1);
    REQUIRE(e.extent(1) == 2);
    REQUIRE(e.rank() == 2);
    REQUIRE(e.rank_dynamic() == 0);
    REQUIRE(sizeof(e) == 1);
  }
  SECTION("mixed") {
    fub::extents<1, fub::dynamic_extent> e(2);
    REQUIRE(e.static_extent(0) == 1);
    REQUIRE(e.static_extent(1) == fub::dynamic_extent);
    REQUIRE(e.extent(0) == 1);
    REQUIRE(e.extent(1) == 2);
    REQUIRE(e.rank() == 2);
    REQUIRE(e.rank_dynamic() == 1);
    REQUIRE(sizeof(e) == sizeof(std::ptrdiff_t));
  }
  SECTION("full") {
    fub::extents<fub::dynamic_extent, fub::dynamic_extent> e(1, 2);
    REQUIRE(e.static_extent(0) == fub::dynamic_extent);
    REQUIRE(e.static_extent(1) == fub::dynamic_extent);
    REQUIRE(e.extent(0) == 1);
    REQUIRE(e.extent(1) == 2);
    REQUIRE(e.rank() == 2);
    REQUIRE(e.rank_dynamic() == 2);
    REQUIRE(sizeof(e) == 2 * sizeof(std::ptrdiff_t));
  }
}

TEST_CASE("Construct mdspan") {
  std::array<int, 8> arr{0, 1, 2, 3, 4, 5, 6, 7};
  fub::DynamicMdSpan<int, 2> view(arr.data(), fub::DynamicExtents<2>(2, 4));
  REQUIRE(view(0, 0) == 0);
  REQUIRE(view(0, 1) == 1);
  REQUIRE(view(0, 2) == 2);
  REQUIRE(view(0, 3) == 3);
  REQUIRE(view(1, 0) == 4);
  REQUIRE(view(1, 1) == 5);
  REQUIRE(view(1, 2) == 6);
  REQUIRE(view(1, 3) == 7);
  view(1, 1) = 55;
  REQUIRE(arr[5] == 55);
}

TEST_CASE("subspan of mdspan") {
  std::array<int, 100> array;
  std::iota(array.begin(), array.end(), 0);
  fub::mdspan<int, 10, 10> ghost_view(array.data());
  SECTION("Keep Dimension") {
    auto inner =
        fub::subspan(ghost_view, std::make_pair(1, 9), std::make_pair(1, 9));
    REQUIRE(inner.rank() == 2);
    REQUIRE(inner.extent(0) == 8);
    REQUIRE(inner.extent(1) == 8);
    for (int i = 0; i < inner.extent(0); ++i) {
      for (int j = 0; j < inner.extent(1); ++j) {
        REQUIRE(inner(i, j) == ghost_view(1 + i, 1 + j));
      }
    }
  }
  SECTION("Reduce Dimension") {
    auto inner = fub::subspan(ghost_view, 8, std::make_pair(1, 9));
    REQUIRE(inner.rank() == 1);
    REQUIRE(inner.extent(0) == 8);
    for (int i = 0; i < inner.extent(0); ++i) {
      REQUIRE(inner(i) == ghost_view(8, 1 + i));
    }

    auto left = fub::subspan(inner, std::make_pair(0, 7));
    auto right = fub::subspan(inner, std::make_pair(1, 8));
    REQUIRE(left.rank() == 1);
    REQUIRE(left.extent(0) == 7);
    REQUIRE(right.rank() == 1);
    REQUIRE(right.extent(0) == 7);
    for (int i = 0; i < 6; ++i) {
      REQUIRE(left(i + 1) == right(i));
    }
  }
}