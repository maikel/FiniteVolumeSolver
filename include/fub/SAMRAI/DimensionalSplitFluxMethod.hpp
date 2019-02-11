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

#ifndef FUB_SAMRAI_DIMENSIONAL_SPLIT_FLUX_METHOD_HPP
#define FUB_SAMRAI_DIMENSIONAL_SPLIT_FLUX_METHOD_HPP

#include "fub/Direction.hpp"

#include "SAMRAI/hier/Patch.h"

#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/pdat/CellData.h>

#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/Dimension.h"

namespace fub {
/// \ingroup Abstract
/// \brief This is an abstract interface for flux methods.
///
/// A flux method has the task to compute flux data between cells. Flux data is
/// normally located at faces.
struct DimensionalSplitFluxMethod {
  virtual ~DimensionalSplitFluxMethod() = default;

  /// \brief Returns a time step value which estimates a maximal stable time
  /// step size.
  ///
  /// The actual time step size may vary from but not exceed the return value.
  /// To make a stable estimation this function may assume a standard forward
  /// euler update for the conservative variables.
  ///
  /// \param[in] states A reference to complete state data.
  /// \param[in] patch  The patch in question.
  virtual double
  ComputeStableDtOnPatch(span<const SAMRAI::pdat::CellData<double>*> states,
                         const SAMRAI::hier::Patch& patch, int dir) const = 0;

  /// \brief Fill the flux data on a patch based on the given states.
  ///
  /// \param[in] fluxes  A reference to the flux data.
  /// \param[in] states  A reference to the complete state data.
  /// \param[in] patch   The patch which may contain additional information
  ///                    related to the area in question.
  /// \param[in] dt      The time step size which will be taken.
  /// \param[in] dir     The split direction.
  virtual void
  ComputeFluxesOnPatch(span<SAMRAI::pdat::SideData<double>*> fluxes,
                       span<const SAMRAI::pdat::CellData<double>*> states,
                       const SAMRAI::hier::Patch& patch, double dt,
                       int dir) const = 0;

  /// \brief Returns the required stencil in each direction.
  ///
  /// Settings this will ensure that ghost cells are filled with the required
  /// amount before the member function `SetPhysicalBoundaryConditions` will be
  /// called.
  ///
  /// \param[in] dim The dimension of the patch hierarchy.
  virtual SAMRAI::hier::IntVector
  GetStencilWidth(const SAMRAI::tbox::Dimension& dim) const = 0;
};

} // namespace fub

#endif