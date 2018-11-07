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

#include "fub/solver/BoundaryCondition.hpp"

#include "SAMRAI/hier/PatchHierarchy.h"

#include <memory>

namespace fub {

struct SystemSolver {
  /// \brief Allocate patch data on a level which is required for this
  /// intergrator.
  ///
  /// This method will be typically called by some
  /// SAMRAI::mesh::GriddingAlgorithm in a intialisation or regridding step.
  ///
  /// \param[in] level  The patch level which patches shall allocate data.
  virtual void
  allocatePatchDataOnPatchLevel(SAMRAI::hier::PatchLevel& level) const = 0;

  /// \brief Estimates the next stable time size.
  ///
  /// \param[in] hierarchy  The patch hierarchy where to operate on.
  /// \param[in] boundary_condition  The boundary condition needed to fill the
  ///                                ghost cell layer of patches.
  /// \param[in] time_point  The current time point of the simulation.
  ///
  /// \return A time step size value.
  virtual double computeStableDt(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      const BoundaryCondition& boundary_condition, double time_point) const = 0;

  /// \brief Advanves a patch hierarchy in time.
  ///
  /// This method applies the splitting method for each spatial diretion of the
  /// hierarchy.
  ///
  /// \param[in,out] hierarchy  The patch hierarchy where to operate on.
  /// \param[in] boundary_condition  The boundary condition needed to fill the
  ///                                ghost cell layer of patches.
  /// \param[in] time_point  The current time point of the simulation.
  /// \param[in] time_step_size  The time step size which will be taken.
  virtual void
  advanceTime(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
              const BoundaryCondition& boundary_condition, double time_point,
              double time_step_size) const = 0;
};

} // namespace fub