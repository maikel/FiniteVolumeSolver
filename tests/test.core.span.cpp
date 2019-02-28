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

#include "fub/core/span.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <boost/hana/tuple.hpp>

#include <vector>

TEST_CASE("Default construct spans") {
  SECTION("Implicit constructor zeroes its elements.") {
    fub::span<int> span;
    REQUIRE(span.size() == 0);
    REQUIRE(span.size_bytes() == 0);
    REQUIRE(span.data() == nullptr);
    REQUIRE(span.empty());
  }
  SECTION("Braced default constructor creates an empty span") {
    fub::span<int> span{};
    REQUIRE(span.size() == 0);
    REQUIRE(span.size_bytes() == 0);
    REQUIRE(span.data() == nullptr);
    REQUIRE(span.empty());
  }
  SECTION("Statically sized spans with size 0 can be default constructed.") {
    fub::span<int, 0> span;
    REQUIRE(span.size() == 0);
    REQUIRE(span.size_bytes() == 0);
    REQUIRE(span.data() == nullptr);
    REQUIRE(span.empty());
  }
  SECTION("Statically sized spans with size 0 can be default constructed, 2.") {
    fub::span<int, 0> span{};
    REQUIRE(span.size() == 0);
    REQUIRE(span.size_bytes() == 0);
    REQUIRE(span.data() == nullptr);
    REQUIRE(span.empty());
  }
  REQUIRE(not std::is_default_constructible<fub::span<int, 1>>::value);
}

TEST_CASE("SFINAE Constructor Test") {
  SECTION("std::array compile time requirements") {
    REQUIRE(std::is_convertible<std::array<int, 3>&, fub::span<int, 3>>());
    REQUIRE(
        std::is_convertible<std::array<int, 3>&, fub::span<const int, 3>>());
    REQUIRE(std::is_convertible<const std::array<int, 3>&,
                                fub::span<const int, 3>>());
  }
  SECTION("std::vector runtime requirements") {
    REQUIRE(std::is_convertible<const std::array<int, 3>&,
                                fub::span<const int, 3>>());
    REQUIRE(std::is_convertible<std::vector<int>&, fub::span<int, 3>>());
  }
}

TEST_CASE("span<T, N> to boost::hana::tuple") {
  std::vector<int> v{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  fub::span<int, 10> span = v;
  auto tuple = boost::hana::to<boost::hana::tuple_tag>(span);
  REQUIRE(boost::hana::at_c<0>(tuple) == span[0]);
  REQUIRE(boost::hana::at_c<1>(tuple) == span[1]);
  REQUIRE(boost::hana::at_c<2>(tuple) == span[2]);
  REQUIRE(boost::hana::at_c<3>(tuple) == span[3]);
  REQUIRE(boost::hana::at_c<4>(tuple) == span[4]);
  REQUIRE(boost::hana::at_c<5>(tuple) == span[5]);
  REQUIRE(boost::hana::at_c<6>(tuple) == span[6]);
  REQUIRE(boost::hana::at_c<7>(tuple) == span[7]);
  REQUIRE(boost::hana::at_c<8>(tuple) == span[8]);
  REQUIRE(boost::hana::at_c<9>(tuple) == span[9]);
}
