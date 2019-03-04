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

#include "fub/ext/Eigen.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

using namespace fub;

TEST_CASE("Load from linear memory") {
  std::array<double, 10> xs{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  mdspan<const double, 1> mdspan(xs.data(), 10);
  for (int i = 0; i < 10; i += 2) {
    Eigen::Array<double, 2, 1> x = Load(constant<2>(), mdspan, i);
    REQUIRE(x[0] == i);
    REQUIRE(x[1] == i + 1);
  }
  Eigen::Array<double, 2, 1> x = Load(constant<2>(), mdspan, 9);
  REQUIRE(x[0] == 9);
  x[0] = 42;

  mdspan<double, 1> span(xs.data(), 10);
  Store(span, x, 9);
  std::array<double, 10> ys{0, 1, 2, 3, 4, 5, 6, 7, 8, 42};
  REQUIRE(xs == ys);
}

TEST_CASE("Load from strided memory") {
  std::array<double, 10> xs{0, 1, 2, 3, 4, 
                            5, 6, 7, 8, 9};
  {
    mdspan<const double, 2> mdspan(xs.data(), 5, 2);
    auto inner = subspan(mdspan, 3, all);
    Eigen::Array<double, 2, 1> x = Load(constant<2>(), inner, 0);
    REQUIRE(x[0] == 3);
    REQUIRE(x[1] == 8);
  }
  {
    mdspan<double, 2> mdspan(xs.data(), 5, 2);
    Eigen::Array<double, 2, 1> x{42, 24};
    auto inner = subspan(mdspan, 3, all);
    Store(inner, x, 0);
  }
  std::array<double, 10> ys{0, 1, 2, 42, 4, 
                            5, 6, 7, 24, 9};
  REQUIRE(xs == ys);
}
