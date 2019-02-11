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

#ifndef FUB_CARTESIAN_COORDINATES_HPP
#define FUB_CARTESIAN_COORDINATES_HPP

#include "fub/Eigen.hpp"
#include "fub/core/mdspan.hpp"

namespace fub {

class CartesianCoordinates {
public:
  CartesianCoordinates(const Eigen::Vector3d& lower,
                       const Eigen::Vector3d& upper,
                       DynamicExtents<3> extents) noexcept;

  Eigen::Vector3d operator()(std::ptrdiff_t i, std::ptrdiff_t j = 0,
                             std::ptrdiff_t k = 0) const;

  Eigen::Vector3d CellCoordinates(std::ptrdiff_t i, std::ptrdiff_t j,
                                  std::ptrdiff_t k) const;

  Eigen::Vector3d NodeCoordinates(std::ptrdiff_t i, std::ptrdiff_t j,
                                  std::ptrdiff_t k) const;

  const Eigen::Vector3d& dx() const noexcept { return dx_; }

private:
  Eigen::Vector3d lower_;
  Eigen::Vector3d upper_;
  DynamicExtents<3> extents_;
  Eigen::Vector3d dx_;
};

} // namespace fub

#endif