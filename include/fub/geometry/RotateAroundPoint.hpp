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

#ifndef FUB_GEOMETRY_ROTATE_AROUND_POINT_HPP
#define FUB_GEOMETRY_ROTATE_AROUND_POINT_HPP

#include "fub/ext/Eigen.hpp"

namespace fub {

template <typename Base> class RotateAroundPoint : private Base {
public:
  Eigen::Matrix<double, 2, 2> R_;
  Eigen::Vector2d x0_;

  RotateAroundPoint(const Base& base, double angle,
                    std::array<double, 2> origin)
      : Base(base), x0_{origin[0], origin[1]} {
    const double cos_a = std::cos(-angle);
    const double sin_a = std::sin(-angle);
    R_ << +cos_a, -sin_a, +sin_a, cos_a;
  }

  [[nodiscard]] double ComputeDistanceTo(std::array<double, 2> x) const
      noexcept {
    Eigen::Vector2d y = R_ * (Eigen::Vector2d::Map(x.data()) - x0_) + x0_;
    return Base::ComputeDistanceTo({y[0], y[1]});
  }
};

template <typename Base, typename Args> RotateAroundPoint(const Base&, double, Args&&) -> RotateAroundPoint<Base>;

} // namespace fub

#endif