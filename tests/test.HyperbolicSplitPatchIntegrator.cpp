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

#include "fub/HyperbolicPatchIntegrator.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <vector>
#include <limits>

template <typename T, typename... Is>
fub::PatchDataView<T, sizeof...(Is)> MakePdv(fub::span<T> data, Is... extents) {
  fub::mdspan<T, sizeof...(Is)> mdspan(data.data(), extents...);
  return fub::PatchDataView<T, sizeof...(Is)>{mdspan, {}};
}

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

TEST_CASE("sequential is correct") {
  SECTION ("1D: 1 cell") {
    std::vector<double> next_data(1, NaN);
    const std::vector<double> prev_data{1.0};
    const std::vector<double> flux_data{0.3, 0.4};
    fub::IndexBox<2> cells{{0, 0}, {1, 1}};
    fub::IndexBox<2> faces{{0, 0}, {2, 1}};
    auto next = MakePdv(fub::span(next_data), 1, 1).Subview(cells);
    auto prev = MakePdv(fub::span(prev_data), 1, 1).Subview(cells);
    auto flux = MakePdv(fub::span(flux_data), 2, 1).Subview(faces);
    auto integrator = fub::HyperbolicPatchIntegrator(fub::execution::seq);
    integrator.UpdateConservatively(next, prev, flux, fub::Duration(1.0), 1.0,
                                    fub::Direction::X);
    REQUIRE(next_data[0] == Approx(0.9));
  }
  SECTION ("1D: 8 cells") {
    std::vector<double> next_data(8, NaN);
    const std::vector<double> prev_data(8, 1.0);
    const std::vector<double> flux_data{0.3, 0.4, -0.2, 0.0, 0.1, 0.2, 0.0, 0.0, 0.0};
    fub::IndexBox<2> cells{{0, 0}, {8, 1}};
    fub::IndexBox<2> faces{{0, 0}, {9, 1}};
    auto next = MakePdv(fub::span(next_data), 8, 1).Subview(cells);
    auto prev = MakePdv(fub::span(prev_data), 8, 1).Subview(cells);
    auto flux = MakePdv(fub::span(flux_data), 9, 1).Subview(faces);
    auto integrator = fub::HyperbolicPatchIntegrator(fub::execution::seq);
    integrator.UpdateConservatively(next, prev, flux, fub::Duration(1.0), 1.0,
                                    fub::Direction::X);
    for (int i = 0; i < next_data.size(); ++i) {
        REQUIRE(next_data[i] == Approx(prev_data[i] + flux_data[i] - flux_data[i + 1]));
    }
  }
}

TEST_CASE("simd is correct") {
    SECTION ("1D: 1 cell") {
  std::vector<double> next_data(1);
  const std::vector<double> prev_data{1.0};
  const std::vector<double> flux_data{0.3, 0.4};
  fub::IndexBox<2> cells{{0, 0}, {1, 1}};
  fub::IndexBox<2> faces{{0, 0}, {2, 1}};
  auto next = MakePdv(fub::span(next_data), 1, 1).Subview(cells);
  auto prev = MakePdv(fub::span(prev_data), 1, 1).Subview(cells);
  auto flux = MakePdv(fub::span(flux_data), 2, 1).Subview(faces);
  auto integrator = fub::HyperbolicPatchIntegrator(fub::execution::simd);
  integrator.UpdateConservatively(next, prev, flux, fub::Duration(1.0), 1.0,
                                  fub::Direction::X);
  REQUIRE(next_data[0] == Approx(0.9));
    }

      SECTION ("1D: 8 cells") {
    std::vector<double> next_data(8, NaN);
    const std::vector<double> prev_data(8, 1.0);
    const std::vector<double> flux_data{0.3, 0.4, -0.2, 0.0, 0.1, 0.2, 0.0, 0.0, 0.0};
    fub::IndexBox<2> cells{{0, 0}, {8, 1}};
    fub::IndexBox<2> faces{{0, 0}, {9, 1}};
    auto next = MakePdv(fub::span(next_data), 8, 1).Subview(cells);
    auto prev = MakePdv(fub::span(prev_data), 8, 1).Subview(cells);
    auto flux = MakePdv(fub::span(flux_data), 9, 1).Subview(faces);
    auto integrator = fub::HyperbolicPatchIntegrator(fub::execution::simd);
    integrator.UpdateConservatively(next, prev, flux, fub::Duration(1.0), 1.0,
                                    fub::Direction::X);
    for (int i = 0; i < next_data.size(); ++i) {
        REQUIRE(next_data[i] == Approx(prev_data[i] + flux_data[i] - flux_data[i + 1]));
    }
  }
}