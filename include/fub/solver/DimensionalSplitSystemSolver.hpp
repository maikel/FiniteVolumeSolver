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

/// \defgroup Solver Solver
/// Classes and methods related to numerical finite volume solvers.

#ifndef FUB_SOLVER_HYPERBOLIC_SYSTEM_SOLVER_HPP
#define FUB_SOLVER_HYPERBOLIC_SYSTEM_SOLVER_HPP

#include "fub/solver/DimensionalSplitTimeIntegrator.hpp"
#include "fub/solver/SplittingMethod.hpp"
#include "fub/solver/SystemSolver.hpp"

namespace fub {
/// \ingroup Solver
/// \brief This is a solver for hyperbolic systems which uses dimensional
/// splitting.
class DimensionalSplitSystemSolver : public SystemSolver {
public:
  /// \brief Constructs a solver object from an dimensional split integrator and
  /// a splitting method.
  DimensionalSplitSystemSolver(const DimensionalSplitTimeIntegrator& integrator,
                         const SplittingMethod& splitting)
      : integrator_{&integrator}, splitting_method_{&splitting} {}

  void allocatePatchDataOnPatchLevel(SAMRAI::hier::PatchLevel& level) const override;

  /// \brief Estimates the next stable time size.
  ///
  /// \param[in] hierarchy  The patch hierarchy where to operate on.
  /// \param[in] boundary_condition  The boundary condition needed to fill the
  ///                                ghost cell layer of patches.
  /// \param[in] time_point  The current time point of the simulation.
  ///
  /// \return A time step size value.
  double computeStableDt(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      const BoundaryCondition& boundary_condition, double time_point) const override;

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
  void
  advanceTime(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
              const BoundaryCondition& boundary_condition, double time_point,
              double time_step_size) const override;

  /// \brief Returns the underlying time integrator.
  const DimensionalSplitTimeIntegrator& getTimeIntegrator() const noexcept {
    return *integrator_;
  }

  /// \brief Returns the underlying splitting method.
  const SplittingMethod& getSplittingMethod() const noexcept {
    return *splitting_method_;
  }

private:
  /// \brief This is an observing pointer to the underlying split time
  /// integrator.
  ///
  /// \note The developer has to make sure that the integrator has to outlive
  /// this object.
  const DimensionalSplitTimeIntegrator* integrator_;

  /// \brief This is an observing pointer to the underlying split time
  /// integrator.
  ///
  /// \note The developer has to make sure that the integrator has to outlive
  /// this object.
  const SplittingMethod* splitting_method_;
};

} // namespace fub

#endif