#include "fub/grid/AMReX/CartesianGridGeometry.hpp"

namespace fub {
namespace amrex {

CartesianCoordinates
GetCartesianCoordinates(const ::amrex::Geometry& geom,
                        const ::amrex::Box& box) {
  Eigen::Vector3d lower{0, 0, 0,};
  geom.LoNode(box.smallEnd(), lower.data());
  Eigen::Vector3d upper{0, 0, 0};
  geom.HiNode(box.bigEnd(), upper.data());
  Eigen::Vector3d dx{0, 0, 0};
  std::array<std::ptrdiff_t, 3> extents{1, 1, 1};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[i] = box.length(i);
    dx[i] = geom.CellSize(i);
  }
  return CartesianCoordinates(lower, upper, dx, dynamic_extents<3>(extents));
}

}
}