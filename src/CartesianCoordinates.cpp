#include "fub/CartesianCoordinates.hpp"
#include <algorithm>

namespace fub {

CartesianCoordinates::CartesianCoordinates(const Eigen::Vector3d& lower,
                                           const Eigen::Vector3d& upper,
                                           const Eigen::Vector3d& dx,
                                           DynamicExtents<3> extents) noexcept
    : lower_{lower}, upper_{upper}, dx_{dx}, extents_{extents} {}

Eigen::Vector3d CartesianCoordinates::
operator()(std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k) const {
  Eigen::Vector3i i_loc{i, j, k};
  Eigen::Vector3i di{extents_.extent(0) - 1,
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