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

#ifndef FUB_GEOMETRY_POLYGON_HPP
#define FUB_GEOMETRY_POLYGON_HPP

#include "fub/core/span.hpp"
#include <boost/container/pmr/vector.hpp>

namespace fub {

class Polygon {
public:
  using Vector = boost::container::pmr::vector<double>;

  Polygon(Vector xs, Vector ys) : xs_{std::move(xs)}, ys_{std::move(ys)} {}

  [[nodiscard]] span<const double> GetXs() const noexcept { return xs_; }
  [[nodiscard]] span<const double> GetYs() const noexcept { return ys_; }

  [[nodiscard]] double ComputeDistanceTo(double x, double y) const noexcept;
  [[nodiscard]] double ComputeDistanceTo(std::array<double, 2> xs) const
      noexcept {
    return ComputeDistanceTo(xs[0], xs[1]);
  }

private:
  Vector xs_;
  Vector ys_;
};

} // namespace fub

#endif
