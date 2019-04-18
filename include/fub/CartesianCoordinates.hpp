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

#include "fub/core/mdspan.hpp"
#include "fub/ext/Eigen.hpp"

namespace fub {

/// This class handles uniform cartesian cell coordinates.
///
/// It stores lower and upper node coordinate bounds, cell widths and cell
/// extents. You can also compute cell and nodes coordinates with this class.
class CartesianCoordinates {
public:
  /// This constructs a coordinate object given domain bounds, cell width and
  /// cell extents.
  ///
  /// \param[in] lower Coordinates of lowest node.
  /// \param[in] upper Coordinates of the highest node.
  /// \param[in] dx Widths of cells
  /// \param[in] extents Number of cells for all coordinate direction.
  CartesianCoordinates(const Eigen::Vector3d& lower,
                       const Eigen::Vector3d& upper, const Eigen::Vector3d& dx,
                       dynamic_extents<3> extents) noexcept;

  /// Returns coordinates of the mid point of a cell in given cell indices.
  ///
  /// \param[in] i The cell index in X direction
  /// \param[in] j The cell index in Y direction
  /// \param[in] k The cell index in Z direction
  ///
  /// \see Same as CartesianCoordinates::CellCoordinates(i, j, k)
  Eigen::Vector3d operator()(std::ptrdiff_t i, std::ptrdiff_t j = 0,
                             std::ptrdiff_t k = 0) const;

  /// Returns coordinates of the mid point of a cell in given cell indices.
  ///
  /// \param[in] i The cell index in X direction
  /// \param[in] j The cell index in Y direction
  /// \param[in] k The cell index in Z direction
  ///
  /// \return A three-dimensional Eigen::Vector3d.
  ///
  /// \see Same as CartesianCoordinates::operator()(i, j, k)
  Eigen::Vector3d CellCoordinates(std::ptrdiff_t i, std::ptrdiff_t j,
                                  std::ptrdiff_t k) const;

  /// Returns coordinates of a nodes in given nodes indices.
  ///
  /// \param[in] i The cell index in X direction
  /// \param[in] j The cell index in Y direction
  /// \param[in] k The cell index in Z direction
  Eigen::Vector3d NodeCoordinates(std::ptrdiff_t i, std::ptrdiff_t j,
                                  std::ptrdiff_t k) const;

  /// Returns the cell width size.
  const Eigen::Vector3d& dx() const noexcept { return dx_; }

  /// Returns the most upper node coordinates (inclusive).
  const Eigen::Vector3d& GetUpper() const noexcept { return upper_; }

  /// Returns the lowest node coordinates (inclusive).
  const Eigen::Vector3d& GetLower() const noexcept { return lower_; }

private:
  Eigen::Vector3d lower_;
  Eigen::Vector3d upper_;
  Eigen::Vector3d dx_;
  dynamic_extents<3> extents_;
};

} // namespace fub

#endif
