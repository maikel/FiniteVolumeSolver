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

#include "fub/simple_grid/SimpleGrid.hpp"
#include "fub/equations/Advection.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

using namespace fub;

TEST_CASE("Test Construction Size") {
  Advection equation({1.0});
  const int gcw = 1;
  const std::ptrdiff_t n_cells = 200;
  const double dx = 1.0 / n_cells;

  SimpleGrid<Advection> grid(equation, n_cells, gcw, dx);
  
  span<double> mass = grid.states().mass.get_span();
  std::iota(mass.begin(), mass.end(), 0);

  SECTION("states have correct size") {
    REQUIRE(grid.states().mass.rank() == 3);
    REQUIRE(grid.states().mass.stride(0) == 1);
    REQUIRE(grid.states().mass.extent(0) == n_cells);
    REQUIRE(grid.states().mass.extent(1) == 1);
    REQUIRE(grid.states().mass.extent(2) == 1);

    REQUIRE(grid.fluxes().mass.rank() == 3);
    REQUIRE(grid.fluxes().mass.stride(0) == 1);
    REQUIRE(grid.fluxes().mass.extent(0) == n_cells + 1);
    REQUIRE(grid.fluxes().mass.extent(1) == 1);
    REQUIRE(grid.fluxes().mass.extent(2) == 1);

    REQUIRE(grid.scratch().mass.rank() == 3);
    REQUIRE(grid.scratch().mass.stride(0) == 1);
    REQUIRE(grid.scratch().mass.extent(0) == n_cells + 2);
    REQUIRE(grid.scratch().mass.extent(1) == 1);
    REQUIRE(grid.scratch().mass.extent(2) == 1);
  }

  SECTION("view inner of scratch") {
    auto scratch = grid.scratch();
    auto states = grid.states();
    CopyToInnerRegion(scratch, states, gcw);
    auto inner = ViewInnerRegion(scratch, gcw);
    REQUIRE(inner.mass.stride(0) == 1);
    REQUIRE(states.mass.stride(0) == 1);
    REQUIRE(Extents(inner).extent(0) == Extents(states).extent(0));
    auto extents = Extents(inner);
    for (int i = 0; i < extents.extent(0); ++i) {
      REQUIRE(inner.mass(i, 0, 0) == states.mass(i, 0, 0));
      REQUIRE(scratch.mass(i + 1, 0, 0) == states.mass(i, 0, 0));
    }
  }
}