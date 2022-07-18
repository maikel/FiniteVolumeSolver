#include "fub/cutcell_method/HGridReconstruction.hpp"

namespace fub {
namespace {
template <typename Vector>
Vector ComputeReflectedCoordinates(const Vector& offset,
                                   const Vector& boundary_normal) {
  Vector reflected =
      offset - 2.0 * offset.dot(boundary_normal) * boundary_normal;
  return reflected;
}

template <typename T, std::size_t N>
T& at(const std::array<T, N>& array, Direction dir) {
  FUB_ASSERT(std::size_t(dir) < N);
  return array[std::size_t(dir)];
}

template <typename T, std::ptrdiff_t N>
T& at(const span<T, N>& array, Direction dir) {
  FUB_ASSERT(std::size_t(dir) < N);
  return array[std::size_t(dir)];
}

} // namespace

void HGridReconstruction::ComputeRegularHGridValues(
    span<double> hgrid_left, span<double> hgrid_right, const PatchData& patch,
    const Index<2>& face) const {
  const int d = static_cast<int>(patch.dir);
  FUB_ASSERT(0 <= d && d < 2);
  const double unshielded_fraction = patch.geometry.face_fractions[d](face);
  FUB_ASSERT(0.0 < unshielded_fraction);
  Index<2> iL(LeftTo(face, patch.dir));
  Index<2> iR(RightTo(face, patch.dir));

  // If no unshielded fraction exist, we return 0.0.
  if (unshielded_fraction <= 0.0) {
    std::fill(hgrid_left.begin(), hgrid_left.end(), 0.0);
    std::fill(hgrid_right.begin(), hgrid_right.end(), 0.0);
    return;
  }

  const int n_components = patch.data.Extent(2);
  FUB_ASSERT(hgrid_left.size() >= n_components);
  FUB_ASSERT(hgrid_right.size() >= n_components);
  // If the face is completely unshielded we do not need to reconstruct any
  // states.
  if (unshielded_fraction >= 1.0) {
    for (int component = 0; component < n_components; ++component) {
      const double qL = patch.data(iL, component);
      const double qR = patch.data(iR, component);
      hgrid_left[component] = qL;
      hgrid_right[component] = qR;
    }
    return;
  }

  // Now we compute the relevant centroids and do the spatial reconstruction at
  // those coordinates.
  const Eigen::Vector2d mL = GetAbsoluteUnshieldedVolumeCentroid(
      patch.geometry, face, Side::Lower, patch.dir, patch.dx);
  const Eigen::Vector2d xL =
      GetAbsoluteVolumeCentroid(patch.geometry, iL, patch.dx);
  const Eigen::Vector2d mR = GetAbsoluteUnshieldedVolumeCentroid(
      patch.geometry, face, Side::Upper, patch.dir, patch.dx);
  const Eigen::Vector2d xR =
      GetAbsoluteVolumeCentroid(patch.geometry, iR, patch.dx);

  const Eigen::Vector2d dxL = mL - xL;
  const Eigen::Vector2d dxR = mR - xR;
  for (int component = 0; component < patch.data.Extent(2); ++component) {
    const double qL = patch.data(iL, component);
    const double qR = patch.data(iR, component);
    Eigen::Vector2d dqL_dx{patch.gradients_x(iL, component),
                           patch.gradients_y(iL, component)};
    Eigen::Vector2d dqR_dx{patch.gradients_x(iR, component),
                           patch.gradients_y(iR, component)};
    const double wL = qL + dqL_dx.dot(dxL);
    const double wR = qR + dqR_dx.dot(dxR);

    hgrid_left[component] = wL;
    hgrid_right[component] = wR;
  }
}

void HGridReconstruction::ComputeInnerBoundaryHGridValues(span<double> hgrid,
                                       const PatchData& data,
                                       const Index<2>& face) const
{
  
}

// double TotalLength(const AuxiliaryReconstructionData& aux) noexcept {
//   double total = 0.0;
//   for (int i = 0; i < AuxiliaryReconstructionData::kMaxSources; ++i) {
//     total += aux.end[i] - aux.start[i];
//   }
//   return total;
// }

// std::array<double, 2> FindAdmissibleInterval(const Coordinates<Rank>& xM,
//                                              const Coordinates<Rank>& xB,
//                                              const Coordinates<Rank>& slope,
//                                              const Coordinates<Rank>& dx) {
//   double lower = std::numeric_limits<double>::lowest();
//   double upper = std::numeric_limits<double>::max();
//   for (int d = 0; d < Rank; ++d) {
//     if (slope[d] != 0) {
//       // xM - (xB + slope * lambda) <= 0.5 * dx
//       const double lambda_lower = (xM[d] - xB[d] - 0.5 * dx[d]) / slope[d];
//       // (xB + slope * lambda) - xM <= 0.5 * dx
//       const double lambda_upper = (0.5 * dx[d] + xM[d] - xB[d]) / slope[d];
//       if (slope[d] > 0) {
//         upper = std::min(upper, lambda_upper);
//         lower = std::max(lower, lambda_lower);
//       } else {
//         upper = std::min(upper, lambda_lower);
//         lower = std::max(lower, lambda_upper);
//       }
//     } else if (std::abs(xM[d] - xB[d]) > 0.5 * dx[d]) {
//       lower = std::numeric_limits<double>::max();
//       upper = std::numeric_limits<double>::lowest();
//     }
//   }
//   return {lower, upper};
// }

// double Length(const std::array<double, 2>& interval) {
//   return interval[1] - interval[0];
// }

// IndexBox<Rank> MakeIndexBox(const Index<Rank>& signs) noexcept {
//   IndexBox<Rank> box{};
//   for (int d = 0; d < Rank; ++d) {
//     if (signs[d] >= 0) {
//       box.upper[d] = 2 * signs[d] + 1;
//     } else {
//       box.upper[d] = 1;
//       box.lower[d] = 2 * signs[d];
//     }
//   }
//   return box;
// }

// AuxiliaryReconstructionData Sort(const AuxiliaryReconstructionData& aux_data)
// {
//   std::array<int, AuxiliaryReconstructionData::kMaxSources> indices;
//   const int count = aux_data.n_sources;
//   std::iota(indices.begin(), indices.begin() + count, 0);
//   std::sort(indices.begin(), indices.begin() + count, [&](int i, int j) {
//     return aux_data.start[i] < aux_data.start[j];
//   });
//   AuxiliaryReconstructionData sorted_aux{};
//   sorted_aux.slope = aux_data.slope;
//   sorted_aux.xB = aux_data.xB;
//   for (int i = 0; i < count; ++i) {
//     sorted_aux.sources[i] = aux_data.sources[indices[i]];
//     sorted_aux.start[i] = aux_data.start[indices[i]];
//     sorted_aux.end[i] = aux_data.end[indices[i]];
//   }
//   for (int i = count; i < AuxiliaryReconstructionData::kMaxSources; ++i) {
//     sorted_aux.sources[i] = sorted_aux.sources[count - 1];
//     sorted_aux.start[i] = sorted_aux.end[count - 1];
//     sorted_aux.end[i] = sorted_aux.end[count - 1];
//   }
//   sorted_aux.n_sources = count;
//   return sorted_aux;
// }

// std::array<AuxiliaryReconstructionData, 2>
// SplitAt(const AuxiliaryReconstructionData& aux, double dx) {
//   std::array<AuxiliaryReconstructionData, 2> datas{};
//   datas[0].slope = aux.slope;
//   datas[0].xB = aux.xB;
//   datas[1].slope = aux.slope;
//   datas[1].xB = aux.xB;
//   int i = 0;
//   while (i < AuxiliaryReconstructionData::kMaxSources) {
//     datas[0].sources[i] = aux.sources[i];
//     datas[0].start[i] = aux.start[i];
//     datas[0].n_sources = i + 1;
//     if (dx < aux.end[i]) {
//       datas[0].end[i] = dx;
//       break;
//     } else {
//       datas[0].end[i] = aux.end[i];
//     }
//     ++i;
//   }
//   datas[1].sources[0] = aux.sources[i];
//   datas[1].start[0] = dx;
//   datas[1].end[0] = aux.end[i];
//   for (int j = i + 1; j < AuxiliaryReconstructionData::kMaxSources; ++j) {
//     datas[1].sources[j - i] = aux.sources[j];
//     datas[1].start[j - i] = aux.start[j];
//     datas[1].end[j - i] = aux.end[j];
//     datas[1].n_sources = aux.n_sources - datas[0].n_sources;
//   }
//   return datas;
// }

// void SetZero(Conservative& cons) {
//   ForEachComponent([](auto&& u) { u = 0.0; }, cons);
// }

// AuxiliaryReconstructionData GetAuxiliaryData(const Index<Rank>& index,
//                                              const CutCellData<Rank>& geom,
//                                              const Coordinates<Rank>& dx,
//                                              const Coordinates<Rank>& slope,
//                                              double required_length) {
//   Index<Rank> signs{};
//   for (int d = 0; d < Rank; ++d) {
//     signs[d] = (slope[d] > 0.0) - (slope[d] < 0.0);
//   }
//   // Fill aux data for the cut-cell itself
//   const Coordinates<Rank> xB =
//       GetBoundaryCentroid(geom, index).array() * dx.array();
//   AuxiliaryReconstructionData aux_data{};
//   aux_data.slope = slope;
//   aux_data.xB = GetAbsoluteBoundaryCentroid(geom, index, dx);
//   int count = 0;
//   ForEachIndex(MakeIndexBox(signs), [&](auto... is) {
//     Coordinates<Rank> x_i{double(is)...};
//     const Coordinates<Rank> xM = x_i.array() * dx.array();
//     const auto interval = FindAdmissibleInterval(xM, xB, slope, dx);
//     const auto [a, b] = Intersect(interval, {0.0, required_length});
//     if (b - a > 0.0) {
//       FUB_ASSERT(count < AuxiliaryReconstructionData::kMaxSources);
//       Index<Rank> i{is...};
//       Index<Rank> global{};
//       std::transform(i.begin(), i.end(), index.begin(), global.begin(),
//                      std::plus<>());
//       aux_data.sources[count] = global;
//       aux_data.start[count] = a;
//       aux_data.end[count] = b;
//       ++count;
//     }
//   });
//   FUB_ASSERT(count > 0 && count <= AuxiliaryReconstructionData::kMaxSources);
//   aux_data.n_sources = count;
//   AuxiliaryReconstructionData sorted_aux = Sort(aux_data);
//   return sorted_aux;
// }

// AuxiliaryReconstructionData GetAuxiliaryData(const Index<Rank>& index,
//                                              const CutCellData<Rank>& geom,
//                                              const Coordinates<Rank>& dx,
//                                              Direction dir) {
//   const Coordinates<Rank> normal = GetBoundaryNormal(geom, index);
//   const int dir_v = static_cast<int>(dir);
//   const Coordinates<Rank> unit_vector =
//       ((normal[dir_v] < 0) - (normal[dir_v] >= 0)) *
//       Eigen::Matrix<double, Rank, 1>::Unit(dir_v);
//   const Coordinates<Rank> slope =
//       unit_vector - 2.0 * unit_vector.dot(normal) * normal;
//   return GetAuxiliaryData(index, geom, dx, slope, dx[int(dir)]);
// }

} // namespace fub