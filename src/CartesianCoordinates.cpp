#include "fub/CartesianCoordinates.hpp"

namespace fub {

CartesianCoordinates::CartesianCoordinates(const Eigen::Vector3d& lower,
                                           const Eigen::Vector3d& upper,
                                           DynamicExtents<3> extents) noexcept
    : lower_{lower}, upper_{upper}, extents_{extents} {
  dx_ = upper_ - lower_;
  dx_[0] /= extents_.extent(0);
  dx_[1] /= extents_.extent(1);
  dx_[2] /= extents_.extent(2);
}

Eigen::Vector3d CartesianCoordinates::
operator()(std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k) const {
  Eigen::Vector3i i_loc{i, j, k};
  Eigen::Vector3i di{extents_.extent(0), extents_.extent(1),
                     extents_.extent(2)};
  Eigen::Vector3d lambda =
      i_loc.cast<double>().array() / di.cast<double>().array();
  Eigen::Vector3d dx_2 = dx_ / 2;
  Eigen::Vector3d result =
      (lower_ + dx_2).array() *
          (Eigen::Vector3d::Constant(1.0) - lambda).array() +
      (upper_ - dx_2).array() * lambda.array();
  return result;
}

// Eigen::Vector3d CartesianCoordinates::CellCoordinates(std::ptrdiff_t i,
// std::ptrdiff_t j,
//                                 std::ptrdiff_t k) const;

// Eigen::Vector3d CartesianCoordinates::NodeCoordinates(std::ptrdiff_t i,
// std::ptrdiff_t j,
//                                 std::ptrdiff_t k) const;

} // namespace fub