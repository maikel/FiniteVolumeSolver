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

#ifndef FUB_SAMRAI_DIMENSIONAL_SPLIT_INTEGRATOR_HPP
#define FUB_SAMRAI_DIMENSIONAL_SPLIT_INTEGRATOR_HPP

#include "fub/SAMRAI/utility.hpp"
#include "fub/SAMRAI/BoundaryCondition.hpp"
#include "fub/Direction.hpp"
#include "fub/SAMRAI/InitialCondition.hpp"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

namespace fub {
/// \ingroup Abstract
/// \brief This class is an abstract interface for time integrators which will
/// be used in a dimensional split setting.
///
/// The dimensional split intergrators are different from standard ones because
/// they take an additional parameter of type Direction in their methods.
/// That parameter indicates in which split direction a specific method shall
/// take action in.
class DimensionalSplitTimeIntegrator {
public:
  virtual ~DimensionalSplitTimeIntegrator() = default;

  /// \brief Allocate patch data on a level which is required for this
  /// intergrator.
  ///
  /// This method will be typically called by some
  /// SAMRAI::mesh::GriddingAlgorithm in a intialisation or regridding step.
  ///
  /// \param[in] level  The patch level which patches shall allocate data.
  void AllocatePatchDataOnPatchLevel(SAMRAI::hier::PatchLevel& level) const;

  /// \brief Returns a time step value which estimates a maximal stable time
  /// step size for a given hierarchy.
  ///
  /// The actual time step size may vary from but not exceed the return value.
  ///
  /// \param[in, out] hierarchy The patch in question.
  /// \param[in] time_point The current time point of the simulation.
  ///
  /// \return a time step size value.
  double ComputeStableDt(const SAMRAI::hier::PatchHierarchy& hierarchy,
                         double time_point, Direction dir) const;

  /// \brief Fill ghost cells of all patches in the patch hierarchy.
  ///
  /// This function will fill ghost cells on patches across all patch levels in
  /// the hierarchy. Ghost cells which lie outside of the domain will be filled
  /// using the specified boundary condition. Interpolation of ghost cells
  /// between the levels has to be handled via the RefineAlgorithm and
  /// CoarsenAlgorithm returned by \ref GetFillGhostLayerRefineAlgorithm and
  /// \ref GetFillGhostLayerCoarsenAlgorithm.
  ///
  /// \param[in, out] hierarchy      The patch hierarchy in question.
  /// \param[in] boundary_condition  The boundary condition which handles ghost
  ///                                cells outside of the domain.
  /// \param[in] time_point          The current time point of the simulation.
  /// \param[in] dir                 The direction in which the ghost cells
  ///                                shall be filled.
  void
  FillGhostLayer(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                 BoundaryCondition& boundary_condition, double time_point,
                 Direction dir) const;

  /// Performs a split step for all patches in the hierarchy.
  ///
  /// This method advances the given PatchHierarchy `hierarchy` in time from
  /// `time_point` to `time_point + time_step_size` in a dimensional split
  /// fashion.
  ///
  /// \param[in, out] hierarchy  The PatchHierarchy which holds the data.
  /// \param[in] time_point  The initial time to start the time step from.
  /// \param[in] time_step_size  The time step size of the coarsest level in the
  ///                            hierarchy.
  /// \param[in] dir  The direction of the dimensional split step.
  void
  AdvanceTime(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
              double time_point, double time_step_size, Direction dir) const;

private:
  /// \brief Inidicates which patch data ids shall be allocated.
  ///
  /// This function will be used in \ref AllocatePatchDataOnPatchLevel to
  /// allocate the flagged data ids automatically for each new patch.
  ///
  /// \param[out] which_to_allocate The ComponentSelector which shall be
  /// modified.
  virtual void FlagPatchDataIdsToAllocate(
      SAMRAI::hier::ComponentSelector& which_to_allocate) const = 0;

  /// \brief Returns a time step value which estimates a maximal stable time
  /// step size for a given patch.
  ///
  /// The actual time step size may vary from but not exceed the return value.
  ///
  /// \param[in,out] patch  The patch in question.
  /// \param[in] time_point  The current time point of the simulation.
  ///
  /// \return a time step size value.
  virtual double ComputeStableDtOnPatch(const SAMRAI::hier::Patch& patch,
                                        double time_point,
                                        Direction dir) const = 0;

  /// \brief Advances a specified patch in time.
  ///
  /// An concrete implementation shall do its time integration of patch in this
  /// function. You may assume that the neccessary ghost cells are filled and
  /// interpolated between patch levels.
  ///
  /// \param[in,out] patch       The patch where data and geomerty lives on.
  /// \param[in] time_point      The current time point of the simulation.
  /// \param[in] time_step_size  The total time step size which shall be taken.
  /// \param[in] dir             The split direction in which shall be
  ///                            integrated.
  virtual void AdvanceTimeOnPatch(const SAMRAI::hier::Patch& patch,
                                  double time_point, double time_step_size,
                                  Direction dir) const = 0;

  /// \brief Returns a RefineAlgorithm which describes which patch data and how
  /// that data is getting synchronised on ghost cells.
  ///
  /// This method is called from \ref FillGhostLayer and it will be used for
  /// communication pass from coarse to fine patch levels.
  ///
  /// \param[in] dir  The split direction in which to fill ghost cell values.
  ///
  /// \return a RefineAlgorithm
  virtual std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>
  GetFillGhostLayerRefineAlgorithm(Direction dir) const = 0;

  /// \brief Returns a CoarsenAlgorithm which describes which patch data and how
  /// that data is getting synchronised on ghost cells.
  ///
  /// This method is called from \ref FillGhostLayer and it will be used for
  /// communication pass from fine to coarse patch levels.
  ///
  /// \param[in] dir  The split direction in which to fill ghost cell values.
  ///
  /// \return a CoarsenAlgorithm
  virtual std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
  GetFillGhostLayerCoarsenAlgorithm(Direction dir) const = 0;

  /// \brief Apply any post processing for a patch level.
  ///
  /// This method will be called after all patches of the given patch level have
  /// been advanced via prior calls to \ref AdvanceTimeOnPatch.
  ///
  /// This method may be used for any additional correction passes. For example
  /// if the concrete implementation needs a conservation fixup between two
  /// patch levels.
  ///
  /// \param[in,out] level
  ///
  /// \note A default no-op implementation is given for convenience such that
  /// the user does not have to provide this method.
  virtual void PostprocessAdvanceLevelState(
      const std::shared_ptr<SAMRAI::hier::PatchLevel>& /* level */) const {}
};

/// \ingroup Solver
/// \brief Initialize a hierarchy and allocate data required by the integrator.
///
/// This method uses the SAMRAI::mesh::GriddingAlgorithm class with some default
/// strategies to generate a coarsest patch level on `hierarchy`. It will
/// construct the level using the following strategies:
///
/// - SAMRAI::mesh::TileClustering for the BoxGenerator
/// - SAMRAI::mesh::CascadePartitioner for the LoadBalancer
///
/// \param[out] hierarchy The hierarchy where the coarsest level will be
///             generated.
/// \param[in] integrator The integrator which tells which variabels to allocate
/// \param[in] initial_condition The initial condition sets values of newly
///            generated patches.
void InitializePatchHierarchy(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    const DimensionalSplitTimeIntegrator& integrator,
    const InitialCondition& initial_condition);

} // namespace fub

#endif