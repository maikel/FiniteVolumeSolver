#ifndef FUB_CUTCELL_H_GRID_RECONSTRUCTION_HPP
#define FUB_CUTCELL_H_GRID_RECONSTRUCTION_HPP

#include "fub/CutCellData.hpp"
#include "fub/Direction.hpp"
#include "fub/PatchDataView.hpp"
#include "fub/ext/Eigen.hpp"
#include <array>

namespace fub {

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