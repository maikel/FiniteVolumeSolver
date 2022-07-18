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

namespace fub {

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
