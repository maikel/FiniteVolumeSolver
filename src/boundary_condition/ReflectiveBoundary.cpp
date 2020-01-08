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

#include "fub/boundary_condition/ReflectiveBoundary.hpp"

namespace fub {
namespace {
template <int Rank>
Index<Rank> ReflectIndex_(Index<Rank> i, const IndexBox<Rank>& fillbox,
                          Direction dir, int side) {
  const std::size_t d = static_cast<std::size_t>(dir);
  if (side == 0) {
    std::ptrdiff_t distance = fillbox.upper[d] - i[d];
    i[d] += 2 * distance - 1;
  } else {
    std::ptrdiff_t distance = i[d] - (fillbox.lower[d] - 1);
    i[d] = i[d] - 2 * distance + 1;
  }
  return i;
}
} // namespace

std::array<std::ptrdiff_t, 1> ReflectIndex(std::array<std::ptrdiff_t, 1> i,
                                           const IndexBox<1>& domain,
                                           Direction dir, int side) {
  return ReflectIndex_(i, domain, dir, side);
}
std::array<std::ptrdiff_t, 2> ReflectIndex(std::array<std::ptrdiff_t, 2> i,
                                           const IndexBox<2>& domain,
                                           Direction dir, int side) {
  return ReflectIndex_(i, domain, dir, side);
}
std::array<std::ptrdiff_t, 3> ReflectIndex(std::array<std::ptrdiff_t, 3> i,
                                           const IndexBox<3>& domain,
                                           Direction dir, int side) {
  return ReflectIndex_(i, domain, dir, side);
}

} // namespace fub
