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

#include "fub/geometry/Polygon.hpp"
#include "fub/ext/Eigen.hpp"

namespace fub {
namespace {

struct Point {
  Array1d x;
  Array1d y;
};

Array1d Dot_(const Point& p, const Point& q) noexcept {
  return p.x * q.x + p.y * q.y;
}

Array1d Distance_(const Point& p, const Point& q) noexcept {
  Point diff{q.x - p.x, q.y - p.y};
  return Dot_(diff, diff).sqrt();
}

Array1d Distance_(const Point& v, const Point& w, const Point& p) {
  const Point w_minus_v{w.x - v.x, w.y - v.y};
  const Array1d length_squared = Dot_(w_minus_v, w_minus_v);
  const Point p_minus_v{p.x - v.x, p.y - v.y};

  const Array1d t = (Dot_(p_minus_v, w_minus_v) / length_squared)
                        .min(Array1d::Constant(1.0))
                        .max(Array1d::Zero());

  const Point proj{(Array1d::Constant(1.0) - t) * v.x + t * w.x,
                   (Array1d::Constant(1.0) - t) * v.y + t * w.y};

  const Array1d distance =
      (length_squared > 0.0).select(Distance_(p, proj), Distance_(p, v));

  return distance;
}

Array1d DistanceN_(Point v, Point w, Point p, int n) {
  Array1d distance = Distance_(v, w, p);
  constexpr double inf = std::numeric_limits<double>::infinity();
  for (int i = n; i < kDefaultChunkSize; ++i) {
    distance[i] = inf;
  }
  return distance;
}

Array<int, 1> CountCrossedLines_(const Point& v, const Point& w,
                                 const Point& p) {
  const Array1d t = ((p.y - v.y) / (w.y - v.y));
  const Array1d xIntercept = (Array1d::Constant(1.0) - t) * v.x + t * w.x;
  Array<bool, 1> in_range =
      ((v.y < p.y && p.y <= w.y) || (w.y < p.y && p.y <= v.y)) &&
      (xIntercept < p.x);
  const Array<int, 1> count =
      in_range.select(Array<int, 1>::Constant(1), Array<int, 1>::Constant(0));
  return count;
}

Array<int, 1> CountCrossedLinesN_(const Point& v, const Point& w,
                                  const Point& p, int n) {
  Array<int, 1> count = CountCrossedLines_(v, w, p);
  for (int i = n; i < kDefaultChunkSize; ++i) {
    count[i] = 0;
  }
  return count;
}

} // namespace

double Polygon::ComputeDistanceTo(double x0, double y0) const noexcept {
  const double* first_x = xs_.data();
  const double* first_y = ys_.data();
  const double* last_x = xs_.data() + xs_.size();

  Array1d min_distance =
      Array1d::Constant(std::numeric_limits<double>::infinity());

  Array<int, 1> num_crossed_lines = Array<int, 1>::Zero();

  const Point p{Array1d::Constant(x0), Array1d::Constant(y0)};

  std::ptrdiff_t count = last_x - first_x;
  while (count > static_cast<std::ptrdiff_t>(kDefaultChunkSize + 1)) {
    const Point v{Array1d::Map(first_x), Array1d::Map(first_y)};
    const Point w{Array1d::Map(first_x + 1), Array1d::Map(first_y + 1)};

    num_crossed_lines += CountCrossedLines_(v, w, p);
    min_distance = min_distance.min(Distance_(v, w, p));

    first_x += kDefaultChunkSize;
    first_y += kDefaultChunkSize;
    count = last_x - first_x;
  }

  Point v{Array1d::Zero(), Array1d::Zero()};
  Point w{Array1d::Zero(), Array1d::Zero()};
  const int n = static_cast<int>(count - 1);
  for (int i = 0; i < n; ++i) {
    v.x[i] = first_x[i];
    v.y[i] = first_y[i];
    w.x[i] = first_x[i + 1];
    w.y[i] = first_y[i + 1];
  }
  num_crossed_lines += CountCrossedLinesN_(v, w, p, n);
  min_distance = min_distance.min(DistanceN_(v, w, p, n));

  const double distance = min_distance.minCoeff();
  const int crossed_lines = num_crossed_lines.sum();
  const int sign = crossed_lines % 2 ? 1 : -1;
  return sign * distance;
}

} // namespace fub
