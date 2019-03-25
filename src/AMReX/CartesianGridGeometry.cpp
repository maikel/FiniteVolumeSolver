#include "fub/grid/AMReX/CartesianGridGeometry.hpp"

namespace fub {
namespace amrex {

CartesianCoordinates GetCartesianCoordinates(const ::amrex::Geometry& geom,
                                             const ::amrex::Box&) {
  Eigen::Vector3d lower{0, 0, 0};
  Eigen::Vector3d upper{0, 0, 0};
  Eigen::Vector3d dx{0, 0, 0};
  std::array<std::ptrdiff_t, 3> extents{1, 1, 1};
  const ::amrex::Box& box = geom.Domain();
  geom.LoNode(box.smallEnd(), &lower[0]);
  geom.HiNode(box.bigEnd(), &upper[0]);
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[static_cast<std::size_t>(i)] = box.length(i);
    dx[i] = geom.CellSize(i);
  }
  return CartesianCoordinates(lower, upper, dx, dynamic_extents<3>(extents));
}

} // namespace amrex
} // namespace fub
