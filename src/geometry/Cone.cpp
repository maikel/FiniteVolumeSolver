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

#include "fub/geometry/Cone.hpp"

#include <Eigen/Dense>

#include <cmath>

namespace fub {

namespace {
Eigen::Vector3d AsVector(const std::array<double, 3>& x) {
  return {x[0], x[1], x[2]};
}
} // namespace

ConeIF::ConeIF(const std::array<double, 3>& point, double base_radius,
               double height, bool in)
    : base_point_(point), base_radius_(base_radius), height_(height),
      inside_(in) {}

double ConeIF::ComputeDistanceTo(const std::array<double, 3>& x) const
    noexcept {
  const Eigen::Vector3d x0 = AsVector(base_point_);
  const Eigen::Vector3d p = AsVector(x);
  const Eigen::Vector3d normal = AsVector({1.0, 0.0, 0.0});
  const int sign = inside_ - !inside_;

  // Distance to point p along the cone axis
  const double l = normal.dot(p - x0);
  // Distance to point p at l
  const Eigen::Vector3d dp_radial = p - (x0 + l * normal);
  const double r = dp_radial.norm();
  // Radius at l
  const double radius_at_l = (1.0 - l / height_) * base_radius_;

  if (l < 0.0) {
    if (r < base_radius_) {
      return sign * -l;
    }
    return sign * std::sqrt(l * l + (x0 - (p - l * normal)).squaredNorm());
  }

  Eigen::Vector3d apex = x0 + height_ * normal;
  if (height_ < l) {
    return sign * (p - apex).norm();
  }

  const Eigen::Vector3d rim_point = [&]() -> Eigen::Vector3d {
    if (dp_radial.norm() > std::sqrt(std::numeric_limits<double>::epsilon())) {
      return x0 + dp_radial.normalized() * base_radius_;
    }
    return x0 + base_radius_ * Eigen::Vector3d{0.0, 1.0, 0.0};
  }();

  Eigen::Vector3d face_tangential = (apex - rim_point).normalized();
  const double distance_on_face = face_tangential.dot(p - rim_point);

  if (distance_on_face < 0.0) {
    return sign * std::sqrt(l * l + (p - rim_point).squaredNorm());
  }

  Eigen::Vector3d p_on_face = rim_point + distance_on_face * face_tangential;
  if (radius_at_l < r) {
    return sign * (p_on_face - p).norm();
  }

  return sign * -std::min((p_on_face - p).norm(), l);
}

} // namespace fub