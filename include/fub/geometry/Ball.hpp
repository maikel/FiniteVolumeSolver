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

#ifndef FUB_GEOMETRY_BALL_HPP
#define FUB_GEOMETRY_BALL_HPP

#include <algorithm>
#include <array>

namespace fub {

template <int Rank> class Ball {
public:
  Ball(const std::array<double, Rank>& offset, double radius);

  double ComputeDistanceTo(const std::array<double, Rank>& x) const {
    std::array<double, Rank> dist;
    std::transform(x.begin(), x.end(), offset_.begin(), dist.begin(),
                   [](double x, double y) { return x - y; });
    return std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0) <
           radius_;
  }

private:
  std::array<double, Rank> offset_;
  double radius_;
};

template <std::size_t Rank>
Ball(const std::array<double, Rank>&, double)->Ball<Rank>;

} // namespace fub

#endif