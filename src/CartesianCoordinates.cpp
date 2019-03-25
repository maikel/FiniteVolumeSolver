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

#include "fub/CartesianCoordinates.hpp"
#include <algorithm>

namespace fub {

CartesianCoordinates::CartesianCoordinates(const Eigen::Vector3d& lower,
                                           const Eigen::Vector3d& upper,
                                           const Eigen::Vector3d& dx,
                                           dynamic_extents<3> extents) noexcept
    : lower_{lower}, upper_{upper}, dx_{dx}, extents_{extents} {}

Eigen::Vector3d CartesianCoordinates::
operator()(std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k) const {
  Eigen::Matrix<std::ptrdiff_t, 3, 1> i_loc{i, j, k};
  Eigen::Matrix<std::ptrdiff_t, 3, 1> di{extents_.extent(0) - 1,
                                      std::max(1L, extents_.extent(1) - 1),
                                      std::max(1L, extents_.extent(2) - 1)};
  Eigen::Array3d lambda =
      i_loc.cast<double>().array() / di.cast<double>().array();
  Eigen::Array3d dx_2 = dx_.array() / 2;
  Eigen::Vector3d result =
      (lower_.array() + dx_2) * (Eigen::Array3d::Constant(1.0) - lambda) +
      (upper_.array() - dx_2) * lambda;
  return result;
}

// Eigen::Vector3d CartesianCoordinates::CellCoordinates(std::ptrdiff_t i,
// std::ptrdiff_t j,
//                                 std::ptrdiff_t k) const;

// Eigen::Vector3d CartesianCoordinates::NodeCoordinates(std::ptrdiff_t i,
// std::ptrdiff_t j,
//                                 std::ptrdiff_t k) const;

} // namespace fub