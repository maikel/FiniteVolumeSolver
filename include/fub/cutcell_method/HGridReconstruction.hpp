#ifndef FUB_CUTCELL_H_GRID_RECONSTRUCTION_HPP
#define FUB_CUTCELL_H_GRID_RECONSTRUCTION_HPP

#include "fub/CutCellData.hpp"
#include "fub/Direction.hpp"
#include "fub/PatchDataView.hpp"
#include "fub/ext/Eigen.hpp"
#include <array>

namespace fub {
// This struct computes indices that are used to integrator mirror data on an
// embedded boundary. This mirror data is being used in the first order
// reconstruction of the half time step approximation.
template <int Rank> class AuxiliaryReconstructionData {
public:
  AuxiliaryReconstructionData(const Index<Rank>& index,
                              const CutCellData<Rank>& geom,
                              const Coordinates<Rank>& dx, Direction dir);

  AuxiliaryReconstructionData(const Index<Rank>& index,
                              const CutCellData<Rank>& geom,
                              const Coordinates<Rank>& dx, Direction dir,
                              double total_length);

  std::array<AuxiliaryReconstructionData, 2> SplitAt(double length);

  span<const Index<Rank>> GetSources() const noexcept;
  span<const double> GetStarts() const noexcept;
  span<const double> GetEnds() const noexcept;

  const Coordinates<Rank>& GetSlope() const noexcept;
  const Coordinates<Rank>& GetVolumeCentroid() const noexcept;

private:
  static constexpr int kMaxSources = 6;
  AuxiliaryReconstructionData() = default;

  std::array<Index<Rank>, kMaxSources> sources{};
  std::array<double, kMaxSources> start{};
  std::array<double, kMaxSources> end{};

  int n_sources{};

  Coordinates<Rank> slope{};
  Coordinates<Rank> xB{};
};

struct HGridReconstruction {
  struct PatchData {
    PatchDataView<const double, 3> data;
    PatchDataView<const double, 3> gradients_x;
    PatchDataView<const double, 3> gradients_y;
    CutCellData<2> geometry;
    Eigen::Vector2d dx;
    Direction dir;
  };

  /// \brief Computes h-grid cell states on a regular face fraction.
  ///
  /// \param[in] i  The cell index of a cut-cell
  /// \param[in] data  The patch data view of a component, such as density
  /// \param[in] geom  Geometry informations for the embedded boundary
  /// \param[in] dir  The relevant split direction, either X, Y or Z.
  ///
  /// \return  Returns states that are left and right to the face that is
  /// opposed to the embedded boundary but only for the unshielded face
  /// fraction.
  ///
  /// \see ComputeSinglyShieldedHGridValues to compute the h-grid states left
  /// and right to the singly shielded face fraction.
  ///
  /// \see ComputeBoundaryHGridValues to compute the h-grid states left and
  /// right to the embedded boundary itself.
  void ComputeRegularHGridValues(span<double> hgrid_left,
                                 span<double> hgrid_right,
                                 const PatchData& data,
                                 const Index<2>& face) const;

  /// \brief Computes h-grid cell states on a singly shielded face fraction.
  ///
  /// \param[in] i  The cell index of a cut-cell
  /// \param[in] data  The patch data view of a component, such as density
  /// \param[in] geom  Geometry informations for the embedded boundary
  /// \param[in] dir  The relevant split direction, either X, Y or Z.
  ///
  /// \return  Returns states that are left and right to the face that is
  /// opposed to the embedded boundary.
  ///
  /// \see ComputeBoundaryHGridValues to compute the h-grid states left and
  /// right to the embedded boundary itself.
  void ComputeSinglyShieldedHGridValues(const PatchData& data,
                                        const Index<2>& face) const;

  /// \brief Computes h-grid cell states on an embedded boundary.
  ///
  /// \param[in] i  The cell index of a cut-cell
  /// \param[in] data  The patch data view of a component, such as density
  /// \param[in] geom  Geometry informations for the embedded boundary
  /// \param[in] dir  The relevant split direction, either X, Y or Z.
  ///
  /// \return  Returns states that are left and right to the face that is
  /// opposed to the embedded boundary.
  ///
  /// \see ComputeBoundaryHGridValues to compute the h-grid states left and
  /// right to the embedded boundary itself.
  void ComputeInnerBoundaryHGridValues(span<double> hgrid,
                                       const PatchData& data,
                                       const Index<2>& face) const;
};
} // namespace fub

#endif