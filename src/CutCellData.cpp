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

#include "fub/CutCellData.hpp"
#include "fub/PatchDataView.hpp"

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>

namespace fub {

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_2 = K::Point_2;
using Polygon_2 = CGAL::Polygon_2<K>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;


HGridIntegrationPoints<2> GetHGridIntegrationPoints(const CutCellData<2>& geom, Index<2> index, const Coordinates<2>& dx, Direction dir)
{
  HGridIntegrationPoints<2> integration_points{};
  integration_points.iB = index;
  integration_points.xB = GetAbsoluteBoundaryCentroid(geom, index, dx);
  const auto d = static_cast<std::size_t>(dir);
  for (int i = 0; i < HGridIntegrationPoints<2>::kMaxSources; ++i) {
    integration_points.volume[i] = geom.hgrid_integration_points[d](index, i*5 + 0);
    integration_points.xM[i][0] = geom.hgrid_integration_points[d](index, i*5 + 1);
    integration_points.xM[i][1] = geom.hgrid_integration_points[d](index, i*5 + 2);

    const double* ptr_x = &geom.hgrid_integration_points[d](index, i*5 + 3);
    std::memcpy(&integration_points.index[i][0], ptr_x, sizeof(double));

    const double* ptr_y = &geom.hgrid_integration_points[d](index, i*5 + 4);
    std::memcpy(&integration_points.index[i][1], ptr_y, sizeof(double));
  }
  return integration_points;
}

Coordinates<2> Centroid(const HGridIntegrationPoints<2>& integration)
{
  K::FT total_area = 0;
  K::FT total_center_x = 0;
  K::FT total_center_y = 0;
  constexpr auto size = HGridIntegrationPoints<2>::kMaxSources;
  for (int i = 0; i < size && integration.volume[i] > 0.0; ++i)
  {
    total_center_x += integration.volume[i] * integration.xM[i][0];
    total_center_y += integration.volume[i] * integration.xM[i][1];
    total_area += integration.volume[i];
  }
  total_center_x /= total_area;
  total_center_y /= total_area;
  const double x = total_center_x.exact().convert_to<double>();
  const double y = total_center_y.exact().convert_to<double>();
  return Coordinates<2>{x, y};
}

Eigen::Vector2d GetBoundaryNormal(const CutCellData<2>& ccdata,
                                  const std::array<std::ptrdiff_t, 2>& index) {
  Eigen::Vector2d normal;
  normal[0] = -ccdata.boundary_normals(index[0], index[1], 0);
  normal[1] = -ccdata.boundary_normals(index[0], index[1], 1);
  return normal;
}

Eigen::Vector3d GetBoundaryNormal(const CutCellData<3>& ccdata,
                                  const std::array<std::ptrdiff_t, 3>& index) {
  Eigen::Vector3d normal;
  normal[0] = -ccdata.boundary_normals(index[0], index[1], index[2], 0);
  normal[1] = -ccdata.boundary_normals(index[0], index[1], index[2], 1);
  normal[2] = -ccdata.boundary_normals(index[0], index[1], index[2], 2);
  return normal;
}

Eigen::Vector2d GetBoundaryCentroid(const CutCellData<2>& ccdata,
                                    const std::array<std::ptrdiff_t, 2>& index) {
  Eigen::Vector2d centroid;
  centroid[0] = ccdata.boundary_centeroids(index[0], index[1], 0);
  centroid[1] = ccdata.boundary_centeroids(index[0], index[1], 1);
  return centroid;
}

Eigen::Vector3d GetBoundaryCentroid(const CutCellData<3>& ccdata,
                                    const std::array<std::ptrdiff_t, 3>& index) {
  Eigen::Vector3d centroid;
  centroid[0] = ccdata.boundary_centeroids(index[0], index[1], index[2], 0);
  centroid[1] = ccdata.boundary_centeroids(index[0], index[1], index[2], 1);
  centroid[2] = ccdata.boundary_centeroids(index[0], index[1], index[2], 2);
  return centroid;
}

Eigen::Vector2d GetVolumeCentroid(const CutCellData<2>& ccdata,
                                    const std::array<std::ptrdiff_t, 2>& index) {
  Eigen::Vector2d centroid;
  centroid[0] = ccdata.volume_centeroid(index[0], index[1], 0);
  centroid[1] = ccdata.volume_centeroid(index[0], index[1], 1);
  return centroid;
}

Eigen::Vector3d GetVolumeCentroid(const CutCellData<3>& ccdata,
                                    const std::array<std::ptrdiff_t, 3>& index) {
  Eigen::Vector3d centroid;
  centroid[0] = ccdata.volume_centeroid(index[0], index[1], index[2], 0);
  centroid[1] = ccdata.volume_centeroid(index[0], index[1], index[2], 1);
  centroid[2] = ccdata.volume_centeroid(index[0], index[1], index[2], 2);
  return centroid;
}

Eigen::Vector2d GetAbsoluteUnshieldedCentroid(const CutCellData<2>& geom,
                                      const Index<2>& face,
                                      const Eigen::Vector2d& dx,
                                      Direction dir)
{
  const int dir_x = static_cast<int>(dir);
  const int dir_y = dir_x == 0 ? 1 : 0;
  const Eigen::Vector2d node_offset = GetOffset(face).array() * dx.array();
  Eigen::Vector2d face_offset = 0.5 * dx;
  face_offset[dir_x] = 0.0;
  const Eigen::Vector2d regular_centroid = node_offset + face_offset;
  const double betaSL = geom.shielded_left_fractions[dir_x](face[0], face[1]);
  const double betaSR = geom.shielded_right_fractions[dir_x](face[0], face[1]);
  
  // find unshielded from left interval
  const Index<2> iL = LeftTo(face, dir);
  std::array<double, 2> betaUL_interval{0.0, 1.0};
  if (betaSL > 0.0) {
    const double ny = geom.boundary_normals(iL[0], iL[1], dir_y);
    if (ny < 0.0) {
      betaUL_interval[0] = betaSL;
    } else {
      betaUL_interval[1] = 1.0 - betaSL;
    }
  }
  // find unshielded from right interval
  const Index<2> iR = RightTo(face, dir);
  std::array<double, 2> betaUR_interval{0.0, 1.0};
  if (betaSR > 0.0) {
    const double ny = geom.boundary_normals(iR[0], iR[1], dir_y);
    if (ny < 0.0) {
      betaUR_interval[0] = betaSR;
    } else {
      betaUR_interval[1] = 1.0 - betaSR;
    }
  }

  std::array<double, 2> betaU_interval = Intersect(betaUL_interval, betaUR_interval);
  const double centroid_y = 0.5 * (betaU_interval[0] + betaU_interval[1]) - 0.5;
  Eigen::Vector2d centroid = regular_centroid;
  centroid[dir_y] = regular_centroid[dir_y] + centroid_y * dx[dir_y];
  return centroid;
}

Eigen::Vector2d
GetAbsoluteSinglyShieldedFromRightVolumeCentroid(const CutCellData<2>& geom,
                                    const Index<2>& face, Side side,
                                    Direction dir, const Eigen::Vector2d& dx)
{
  const Index<2> iR = RightTo(face, dir);
  FUB_ASSERT(IsCutCell(geom, iR));
  Eigen::Vector2d xB = GetAbsoluteBoundaryCentroid(geom, iR, dx);
  const auto d = static_cast<std::size_t>(dir);
  FUB_ASSERT(geom.shielded_right_fractions[d](face) > 0.0);
  const double alpha = 0.5 + geom.boundary_centeroids(iR, d);
  FUB_ASSERT(side == Side::Lower || side == Side::Upper);
  if (side == Side::Upper) {
    xB[d] -= 0.5 * alpha * dx[d];
  } else if (side == Side::Lower) {
    xB[d] -= (0.5 + alpha) * dx[d];
  }
  return xB;
}

Eigen::Vector2d
GetAbsoluteSinglyShieldedFromLeftVolumeCentroid(const CutCellData<2>& geom,
                                    const Index<2>& face, Side side,
                                    Direction dir, const Eigen::Vector2d& dx)
{
  const Index<2> iL = LeftTo(face, dir, 1);
  FUB_ASSERT(IsCutCell(geom, iL));
  Eigen::Vector2d xB = GetAbsoluteBoundaryCentroid(geom, iL, dx);
  const auto d = static_cast<std::size_t>(dir);
  FUB_ASSERT(geom.shielded_left_fractions[d](face) > 0.0);
  const double alpha = 0.5 - geom.boundary_centeroids(iL, d);
  FUB_ASSERT(side == Side::Lower || side == Side::Upper);
  if (side == Side::Lower) {
    xB[d] += 0.5 * alpha * dx[d];
  } else if (side == Side::Upper) {
    xB[d] += (0.5 + alpha) * dx[d];
  }
  return xB;
}

namespace {
template <int Rank>
bool IsCutCell_(const CutCellData<Rank>& geom,
                const std::array<std::ptrdiff_t, static_cast<std::size_t>(Rank)>& index) {
  if (geom.volume_fractions(index) <= 0.0) {
    return false;
  }
  if (geom.volume_fractions(index) < 1.0) {
    return true;
  }
  FUB_ASSERT(geom.volume_fractions(index) == 1.0);
  auto left = index;
  for (std::size_t i = 0; i < static_cast<std::size_t>(Rank); ++i) {
    auto right = Shift(left, Direction(i), 1);
    if (geom.face_fractions[i](left) < 1.0 ||
        geom.face_fractions[i](right) < 1.0) {
      return true;
    }
  }
  return false;
}
} // namespace

bool IsCutCell(const CutCellData<2>& geom,
               const std::array<std::ptrdiff_t, 2>& index) {
  return IsCutCell_<2>(geom, index);
}

bool IsCutCell(const CutCellData<3>& geom,
               const std::array<std::ptrdiff_t, 3>& index) {
  return IsCutCell_<3>(geom, index);
}

} // namespace fub
