// Copyright (c) 2018 Maikel Nadolski
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

#ifndef FUB_SOLVER_BOUNDARY_CONDITION_HPP
#define FUB_SOLVER_BOUNDARY_CONDITION_HPP

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "fub/solver/Direction.hpp"

namespace fub {
/// \ingroup Abstract
/// \brief This is an abstract interface class for general purpose boundary
/// conditions.
struct BoundaryCondition {
  virtual ~BoundaryCondition() = default;

  /// \brief Fill ghost cell values of a given patch.
  ///
  /// This function will be called for each patch which touches the
  /// computational domain. The vector ghost_width_to_fill indicates which
  /// border to fill.
  ///
  /// \param[in,out] patch  The patch which ghost cells have to be filled.
  /// \param[in] fill_time  The time point of the current simulation.
  /// \param[in] ghost_width_to_fill  A vector which indicates the ghost layer
  ///                                 width which has to be filled by this
  ///                                 routine.
  virtual void setPhysicalBoundaryConditions(
      const SAMRAI::hier::Patch& patch, double fill_time,
      const SAMRAI::hier::IntVector& ghost_width_to_fill) const = 0;

  /// \brief Returns the required stencil in each direction.
  ///
  /// Settings this will ensure that ghost cells are filled with the required
  /// amount before the member function `setPhysicalBoundaryConditions` will be
  /// called.
  ///
  /// \param[in] dim The dimension of the patch hierarchy.
  virtual SAMRAI::hier::IntVector
  getStencilWidth(const SAMRAI::tbox::Dimension& dim) const = 0;

  virtual void
  PreFillGhostLayer(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&,
                    double, Direction) {}
};

} // namespace fub

#endif