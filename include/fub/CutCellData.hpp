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

#ifndef FUB_CUTCELL_DATA_HPP
#define FUB_CUTCELL_DATA_HPP

#include "fub/PatchDataView.hpp"
#include "fub/ext/Eigen.hpp"

namespace fub {

inline std::array<double, 2> Intersect(const std::array<double, 2>& i1,
                                       const std::array<double, 2>& i2) {
  return {std::max(i1[0], i2[0]), std::min(i1[1], i2[1])};
}

template <int Rank>
struct HGridIntegrationPoints {
    static constexpr int kMaxSources = 9;

    std::array<Index<Rank>, kMaxSources> index{};
    std::array<double, kMaxSources> volume{};
    std::array<Coordinates<Rank>, kMaxSources> xM{};

    Index<Rank> iB{};
    Coordinates<Rank> xB{};
};

Coordinates<2> Centroid(const HGridIntegrationPoints<2>& polygon);

template <int Rank> struct CutCellData {
  static constexpr auto sRank = static_cast<std::size_t>(Rank);

  // The next member variables are given by AMReX
  PatchDataView<const double, Rank> volume_fractions;
  PatchDataView<const double, Rank + 1> volume_centeroid;
  std::array<PatchDataView<const double, Rank>, sRank> face_fractions;
  PatchDataView<const double, Rank + 1> boundary_normals;
  PatchDataView<const double, Rank + 1> boundary_centeroids;
  // The following members need to be computed from the AMReX EB database
  PatchDataView<const double, Rank + 1> distance_to_boundary_left;
  PatchDataView<const double, Rank + 1> distance_to_boundary_right;
  std::array<PatchDataView<const double, Rank>, sRank> unshielded_fractions;
  std::array<PatchDataView<const double, Rank>, sRank> shielded_left_fractions;
  std::array<PatchDataView<const double, Rank>, sRank> shielded_right_fractions;
  std::array<PatchDataView<const double, Rank>, sRank>
      doubly_shielded_fractions;
  std::array<PatchDataView<const double, Rank>, sRank> unshielded_fractions_rel;
  std::array<PatchDataView<const double, Rank>, sRank>
      shielded_left_fractions_rel;
  std::array<PatchDataView<const double, Rank>, sRank>
      shielded_right_fractions_rel;
  std::array<PatchDataView<const double, Rank>, sRank>
      doubly_shielded_fractions_rel;

  std::array<PatchDataView<const double, Rank + 1>, sRank> hgrid_integration_points;
  PatchDataView<const double, Rank + 1> hgrid_total_inner_integration_points;
};

HGridIntegrationPoints<2> GetHGridIntegrationPoints(const CutCellData<2>& geom, Index<2> index, const Coordinates<2>& dx, Direction dir);

HGridIntegrationPoints<2> GetHGridTotalInnerIntegrationPoints(const CutCellData<2>& geom, Index<2> index, const Coordinates<2>& dx);

[[nodiscard]] bool IsCutCell(const CutCellData<2>& geom,
                             const std::array<std::ptrdiff_t, 2>& index);

[[nodiscard]] bool IsCutCell(const CutCellData<3>& geom,
                             const std::array<std::ptrdiff_t, 3>& index);

[[nodiscard]] Eigen::Vector2d
GetBoundaryNormal(const CutCellData<2>& ccdata,
                  const std::array<std::ptrdiff_t, 2>& index);

[[nodiscard]] Eigen::Vector3d
GetBoundaryNormal(const CutCellData<3>& ccdata,
                  const std::array<std::ptrdiff_t, 3>& index);

[[nodiscard]] Eigen::Vector2d
GetBoundaryCentroid(const CutCellData<2>& ccdata,
                    const std::array<std::ptrdiff_t, 2>& index);

[[nodiscard]] Eigen::Vector3d
GetBoundaryCentroid(const CutCellData<3>& ccdata,
                    const std::array<std::ptrdiff_t, 3>& index);

[[nodiscard]] Eigen::Vector2d
GetVolumeCentroid(const CutCellData<2>& ccdata,
                  const std::array<std::ptrdiff_t, 2>& index);

[[nodiscard]] Eigen::Vector3d
GetVolumeCentroid(const CutCellData<3>& ccdata,
                  const std::array<std::ptrdiff_t, 3>& index);

[[nodiscard]] Eigen::Vector2d GetUnshieldedCentroid(const CutCellData<2>& geom,
                                                    const Index<2>& face,
                                                    Direction dir);

[[nodiscard]] Eigen::Vector2d
GetAbsoluteUnshieldedCentroid(const CutCellData<2>& geom, const Index<2>& face,
                              const Eigen::Vector2d& dx, Direction dir);

[[nodiscard]] Eigen::Vector2d
GetUnshieldedVolumeCentroid(const CutCellData<2>& geom, const Index<2>& face,
                            Side side, Direction dir);

[[nodiscard]] Eigen::Vector2d
GetAbsoluteUnshieldedVolumeCentroid(const CutCellData<2>& geom,
                                    const Index<2>& face, Side side,
                                    Direction dir, const Eigen::Vector2d& dx);

[[nodiscard]] Eigen::Vector2d
GetAbsoluteSinglyShieldedFromRightVolumeCentroid(const CutCellData<2>& geom,
                                    const Index<2>& face, Side side,
                                    Direction dir, const Eigen::Vector2d& dx);

[[nodiscard]] Eigen::Vector2d
GetAbsoluteSinglyShieldedFromLeftVolumeCentroid(const CutCellData<2>& geom,
                                    const Index<2>& face, Side side,
                                    Direction dir, const Eigen::Vector2d& dx);

template <std::size_t Rank>
Eigen::Matrix<double, Rank, 1>
GetOffset(const std::array<std::ptrdiff_t, Rank>& index) {
  static constexpr int iRank = static_cast<int>(Rank);
  Eigen::Matrix<double, iRank, 1> offset;
  for (int i = 0; i < iRank; ++i) {
    offset[i] = static_cast<double>(index[i]);
  }
  return offset;
}

template <int Rank>
Eigen::Matrix<double, Rank, 1>
GetAbsoluteVolumeCentroid(const CutCellData<Rank>& geom,
                          const Index<Rank>& index,
                          const Eigen::Matrix<double, Rank, 1>& dx) {
  const Eigen::Matrix<double, Rank, 1> relative_xM =
      GetVolumeCentroid(geom, index);
  const Eigen::Matrix<double, Rank, 1> offset = GetOffset<Rank>(index);
  const Eigen::Matrix<double, Rank, 1> xM =
      (offset + relative_xM).array() * dx.array() + 0.5 * dx.array();
  return xM;
}

template <int Rank>
Eigen::Matrix<double, Rank, 1>
GetAbsoluteBoundaryCentroid(const CutCellData<Rank>& geom,
                            const Index<Rank>& index,
                            const Eigen::Matrix<double, Rank, 1>& dx) {
  const Eigen::Matrix<double, Rank, 1> relative_xB =
      GetBoundaryCentroid(geom, index);
  const Eigen::Matrix<double, Rank, 1> offset = GetOffset<Rank>(index);
  const Eigen::Matrix<double, Rank, 1> xB =
      (offset + relative_xB).array() * dx.array() + 0.5 * dx.array();
  return xB;
}

} // namespace fub

#endif
