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

#ifndef FUB_SIMPLE_GRID_BOUNDARY_CONDITION_PERIODIC_HPP
#define FUB_SIMPLE_GRID_BOUNDARY_CONDITION_PERIODIC_HPP

#include "fub/Equation.hpp"

namespace fub {
namespace simple {

template <typename Equation> class PeriodicBoundary {
public:
  template <typename T> using StateSpan = StateView<T, Equation>;

  explicit PeriodicBoundary(const Equation& equation) : equation_{equation}
  {}

  void SetPhysicalBoundaryConditions(StateSpan<double> states,
                                     double /* time_point */,
                                     std::ptrdiff_t fill_width,
                                     Direction dir = Direction::X) {
    auto extents = Extents(states);
    const std::ptrdiff_t size = extents.extent(int(dir));
    ForEachRow(extents, [&](auto row) {
      ForEachMember(
          [&](auto variable) {
            constexpr int rank = decltype(variable)::rank();
            if constexpr (rank == 1) {
            for (std::ptrdiff_t dest = 0; dest < fill_width; ++dest) {
              const std::ptrdiff_t source = size - 2 * fill_width + dest;
              variable(dest) = variable(source);
              variable(size - dest - 1) = variable(fill_width + dest);
            }
            } else if (rank == 2) {
             for (std::ptrdiff_t dest = 0; dest < fill_width; ++dest) {
              const std::ptrdiff_t source = size - 2 * fill_width + dest;
              for (int i = 0; i < variable.extent(1); ++i) {
              variable(dest, i) = variable(source, i);
              variable(size - dest - 1, i) = variable(fill_width + dest, i);
              }
            } 
            }
          },
          row);
    }, states);
  }

private:
  Equation equation_;
};

} // namespace simple
} // namespace fub

#endif