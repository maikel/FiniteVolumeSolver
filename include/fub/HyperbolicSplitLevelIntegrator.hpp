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

#ifndef FUB_HYPERBOLIC_SPLIT_LEVEL_INTEGRATOR_HPP
#define FUB_HYPERBOLIC_SPLIT_LEVEL_INTEGRATOR_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/Execution.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/ext/outcome.hpp"

#include <stdexcept>

#include <mpi.h>

namespace fub {
/// This Level Integrator applies a very general anisotropic AMR integration
/// scheme in context of dimensional splitting.
///
/// The time integration is split into multiple intermediate steps where each is
/// supposed to do a certain task. The detailed implementation of these tasks
/// happens in the integrator context object.
template <typename Equation, typename Context>
class HyperbolicSplitLevelIntegrator : private Context {
public:
  HyperbolicSplitLevelIntegrator(const Equation& equation,
                                 const Context& context)
      : Context(context), equation_{equation} {}

  HyperbolicSplitLevelIntegrator(const Equation& equation,
                                 Context&& context) // NOLINT
      : Context(context), equation_{equation} {}

  const Context& GetContext() const noexcept { return *this; }
  Context& GetContext() noexcept { return *this; }

  const Equation& GetEquation() const noexcept { return equation_; }

  using Context::ResetHierarchyConfiguration;

  /// Returns the total refinement ratio between specified coarse to fine level
  /// number.
  int GetTotalRefineRatio(int fine_level, int coarse_level = 0) const;

  /// Returns a stable dt across all levels and in one spatial direction.
  ///
  /// This function takes the refinement level into account.
  /// For stability it is advised to use some additional CFL condition.
  ///
  /// \param[in] dir  The direction into which the time step size is calculated.
  Duration ComputeStableDt(Direction dir);

  /// Advance a specified patch level and all finer levels in the specified
  /// split direction by time dt.
  ///
  /// This method subcycles finer levels.
  ///
  /// \param[in] level_num  An integer denoting the patch level where 0 is the
  /// coarsest level.
  ///
  /// \param[in] direction The dimensional split direction which will be used to
  /// advance.
  ///
  /// \param[in] dt A stable time step size for the level_num-th patch level.
  Result<void, TimeStepTooLarge>
  AdvanceLevel(int this_level, Direction direction, Duration dt, int subcycle);

private:
  Equation equation_;
};

// Implementation

template <typename Eq, typename Context>
int HyperbolicSplitLevelIntegrator<Eq, Context>::GetTotalRefineRatio(
    int fine_level, int coarse_level) const {
  int refine_ratio = 1;
  for (int level = fine_level; level > coarse_level; --level) {
    refine_ratio *= Context::GetRatioToCoarserLevel(level).max();
  }
  return refine_ratio;
}

template <typename Eq, typename Context>
Duration
HyperbolicSplitLevelIntegrator<Eq, Context>::ComputeStableDt(Direction dir) {
  int refine_ratio = 1;
  double coarse_dt = std::numeric_limits<double>::infinity();
  for (int level_num = 0; Context::LevelExists(level_num); ++level_num) {
    if (level_num > 0) {
      Context::FillGhostLayerTwoLevels(level_num, level_num - 1, dir);
    } else {
      Context::FillGhostLayerSingleLevel(level_num, dir);
    }
    const double level_dt = Context::ComputeStableDt(level_num, dir).count();
    refine_ratio *= Context::GetRatioToCoarserLevel(level_num).max();
    coarse_dt = std::min(refine_ratio * level_dt, coarse_dt);
  }
  const MPI_Comm comm = Context::GetMpiCommunicator();
  double global_min_dt{0};
  MPI_Allreduce(&coarse_dt, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, comm);
  return Duration(global_min_dt);
}

template <typename Eq, typename Context>
Result<void, TimeStepTooLarge>
HyperbolicSplitLevelIntegrator<Eq, Context>::AdvanceLevel(int this_level,
                                                          Direction direction,
                                                          Duration dt,
                                                          int subcycle) {
  // PreAdvanceLevel might regrid this and all finer levels.
  // The Context must make sure that scratch data is allocated
  Context::PreAdvanceLevel(this_level, direction, dt, subcycle);

  // Fill the ghost layer which is needed for this split direction
  if (subcycle == 0 && this_level > 0) {
    Context::FillGhostLayerTwoLevels(this_level, this_level - 1, direction);
  } else {
    Context::FillGhostLayerSingleLevel(this_level, direction);
  }

  const double level_dt =
      Context::ComputeStableDt(this_level, direction).count();
  if (level_dt < dt.count()) {
    const int refine_ratio = GetTotalRefineRatio(this_level);
    return TimeStepTooLarge{Duration(refine_ratio * level_dt)};
  }

  // Compute fluxes in the specified direction
  Context::ComputeNumericFluxes(this_level, dt, direction);

  // Accumulate Fluxes on the coarse-fine interface to the next coarser level.
  if (this_level > 0) {
    Context::AccumulateCoarseFineFluxes(this_level, dt, direction);
  }

  // If a finer level exists in the hierarchy, we subcycle that finer level
  // multiple times and use the fine fluxes on coarse-fine interfaces
  const int next_level = this_level + 1;
  if (Context::LevelExists(next_level)) {
    Context::ResetCoarseFineFluxes(next_level, this_level, direction);
    const int refine_ratio =
        Context::GetRatioToCoarserLevel(next_level, direction);
    for (int r = 0; r < refine_ratio; ++r) {
      auto result = AdvanceLevel(next_level, direction, dt / refine_ratio, r);
      if (!result) {
        return result.as_failure();
      }
    }
    Context::ApplyFluxCorrection(next_level, this_level, dt, direction);
  }

  // Use the updated fluxes to update cons variables at the "SCRATCH" context.
  Context::UpdateConservatively(this_level, dt, direction);

  // Coarsen inner regions from next finer level to this level.
  if (Context::LevelExists(next_level)) {
    Context::CoarsenConservatively(next_level, this_level, direction);
  }

  // The conservative update and and the coarsening happened on conservative
  // variables. We have to reconstruct the missing variables in the complete
  // state.
  Context::CompleteFromCons(this_level, dt, direction);

  // Apply any further context related work after advancing this level.
  // This function can also indicate if some error occured.
  // For example the context could detect unphysical states and return a
  // TooLargeTimeStep error condition.
  return Context::PostAdvanceLevel(this_level, direction, dt, subcycle);
}

} // namespace fub

#endif
