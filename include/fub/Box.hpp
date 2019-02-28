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

#ifndef FUB_BOX_HPP
#define FUB_BOX_HPP

#include "fub/Direction.hpp"

#include <array>

namespace fub {

template <int Rank> using Index = std::array<std::ptrdiff_t, Rank>;

template <int Rank> struct Box {
  Index<Rank> lower;
  Index<Rank> upper;
};

template <typename Extents>
Box<Extents::rank()> GetBoundaryBox(const Extents& extents, Location loc, int gcw) {
  Box<Extents::rank()> box{{}, AsArray(extents)};
  const int d = static_cast<int>(loc.direction);
  if (loc.side == 0) {
    box.upper[d] = gcw;
  } else {
    box.lower[d] = box.upper[d] - gcw;
  }
  return box;
}


} // namespace fub

#endif