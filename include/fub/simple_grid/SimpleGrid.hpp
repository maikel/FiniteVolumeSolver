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

#ifndef FUB_SIMPLE_GRID_HPP
#define FUB_SIMPLE_GRID_HPP

#include "fub/Equation.hpp"

#include <vector>

namespace fub {
template <typename Equation> class SimpleGrid {
public:
  using State = typename Equation::State;
  using Cons = typename Equation::Cons;

  template <typename T> using StateSpan = StateView<T, Equation>;
  template <typename T> using ConsSpan = ConsView<T, Equation>;

  SimpleGrid(Equation eq, std::ptrdiff_t n_cells, int ghost_width,
             double cell_width)
      : equation_{eq}, cell_width_{cell_width},
        ghost_width_{ghost_width}, states_{MakeStateArray(equation_, n_cells)},
        scratch_{MakeStateArray(equation_, n_cells + 2 * ghost_width_)},
        fluxes_{MakeConsArray(equation_, n_cells + 1)} {}

  double dx() const noexcept { return cell_width_; }
  int ghost_cell_width() const noexcept { return ghost_width_; }
  double time_point() const noexcept {
    return time_point_;
  }

  void time_point(double t) {
    time_point_ = t;
  }

  StateSpan<double> states() noexcept {
    auto sizes = GetSizes<State>(equation_);
    return Transform(
        overloaded{
            [&](std::vector<double>& state, int n) {
              return basic_mdspan<double, DynamicExtents<4>, layout_left>(
                  state.data(), state.size() / n, 1, 1, n);
            },
            [&](std::vector<double>& state, std::integral_constant<int, 1>) {
              return basic_mdspan<double, DynamicExtents<3>, layout_left>(
                  state.data(), state.size(), 1, 1);
            }},
        states_, sizes);
  }

  StateSpan<const double> states() const noexcept {
    auto sizes = GetSizes<State>(equation_);
    return Transform(
        overloaded{
            [&](const std::vector<double>& state, int n) {
              return basic_mdspan<const double, DynamicExtents<4>, layout_left>(
                  state.data(), state.size() / n, 1, 1, n);
            },
            [&](const std::vector<double>& state, std::integral_constant<int, 1>) {
              return basic_mdspan<const double, DynamicExtents<3>, layout_left>(
                  state.data(), state.size(), 1, 1);
            }},
        states_, sizes);
  }

  StateSpan<double> scratch() noexcept {
    auto sizes = GetSizes<State>(equation_);
    return Transform(
        overloaded{
            [&](std::vector<double>& state, int n) {
              return basic_mdspan<double, DynamicExtents<4>, layout_left>(
                  state.data(), state.size() / n, 1, 1, n);
            },
            [&](std::vector<double>& state, std::integral_constant<int, 1>) {
              return basic_mdspan<double, DynamicExtents<3>, layout_left>(
                  state.data(), state.size(), 1, 1);
            }},
        scratch_, sizes);
  }

  StateSpan<const double> scratch() const noexcept {
    auto sizes = GetSizes<State>(equation_);
    return Transform(
        overloaded{
            [&](const std::vector<double>& state, int n) {
              return basic_mdspan<const double, DynamicExtents<4>, layout_left>(
                  state.data(), state.size() / n, 1, 1, n);
            },
            [&](const std::vector<double>& state, std::integral_constant<int, 1>) {
              return basic_mdspan<const double, DynamicExtents<3>, layout_left>(
                  state.data(), state.size(), 1, 1);
            }},
        scratch_, sizes);
  }

  ConsSpan<double> fluxes() noexcept {
    auto sizes = GetSizes<State>(equation_);
    return Transform(
        overloaded{
            [&](std::vector<double>& state, int n) {
              return basic_mdspan<double, DynamicExtents<4>, layout_left>(
                  state.data(), state.size() / n, 1, 1, n);
            },
            [&](std::vector<double>& state, std::integral_constant<int, 1>) {
              return basic_mdspan<double, DynamicExtents<3>, layout_left>(
                  state.data(), state.size(), 1, 1);
            }},
        fluxes_, sizes);
  }

  ConsSpan<const double> fluxes() const noexcept {
    auto sizes = GetSizes<State>(equation_);
    return Transform(
        overloaded{
            [&](const std::vector<double>& state, int n) {
              return basic_mdspan<const double, DynamicExtents<4>, layout_left>(
                  state.data(), state.size() / n, 1, 1, n);
            },
            [&](const std::vector<double>& state, std::integral_constant<int, 1>) {
              return basic_mdspan<const double, DynamicExtents<3>, layout_left>(
                  state.data(), state.size(), 1, 1);
            }},
        fluxes_, sizes);
  }

private:
  Equation equation_;
  double cell_width_;
  int ghost_width_;
  StateArray<Equation> states_;
  StateArray<Equation> scratch_;
  ConsArray<Equation> fluxes_;
  double time_point_{};
};

template <typename State> auto ViewInnerRegion(const State& state, int gcw) {
  return Transform([gcw](auto mdspan) {
    if constexpr (decltype(mdspan)::rank() == 3) {
      return subspan(mdspan, std::pair{gcw, mdspan.extent(0) - gcw}, all, all);
    } else {
      return subspan(mdspan, std::pair{gcw, mdspan.extent(0) - gcw}, all, all, all);
    }
  }, state);
}

template <typename To, typename From>
void CopyToInnerRegion(const To& dest, const From& source, int gcw) {
  auto inner_dest = ViewInnerRegion(dest, gcw);
  ForEachRow(Extents(inner_dest),
      [](auto dest, auto source) {
        ForEachMember(
            [](auto dest, auto source) {
              FUB_ASSERT(source.extent(0) == dest.extent(0));
              FUB_ASSERT(source.stride(0) == 1);
              FUB_ASSERT(dest.stride(0) == 1);
              std::copy_n(source.data(), source.extent(0), dest.data());
            },
            dest, source);
      },
      inner_dest, source);
}

} // namespace fub

#endif