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
#include "fub/Meta.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/ext/outcome.hpp"
#include "fub/solver/SolverFacade.hpp"

#include <stdexcept>

#include <mpi.h>

namespace fub {

/// This Level Integrator applies a very general AMR integration scheme in
/// context of dimensional splitting.
///
/// The time integration is split into multiple intermediate steps where each is
/// supposed to do a certain task. The detailed implementation of these tasks
/// happens in the integrator context object.
template <typename LevelIntegrator>
class SubcycleFineFirstSolver : public SolverFacade<LevelIntegrator> {
private:
  using Base = SolverFacade<LevelIntegrator>;

public:
  static constexpr int Rank = LevelIntegrator::Rank;

  explicit SubcycleFineFirstSolver(LevelIntegrator&& integrator)
      : Base(std::move(integrator)) {}

  explicit SubcycleFineFirstSolver(const LevelIntegrator& integrator)
      : Base(integrator) {}

  /// Returns the total refinement ratio between specified coarse to fine level
  /// number.
  int GetTotalRefineRatio(int fine_level, int coarse_level = 0) const;

  /// Returns a stable time step size for the coarsest refinement level.
  ///
  /// For stability it is advised to multiply some additional CFL number < 1.0.
  Duration ComputeStableDt();

  /// Advance a specified patch level and all finer levels by time `dt`.
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
  Result<void, TimeStepTooLarge> AdvanceLevel(int level_number, Duration dt,
                                              std::pair<int, int> subcycle);

  Result<void, TimeStepTooLarge> AdvanceHierarchy(Duration dt);
};

// Implementation

template <typename LevelIntegrator>
int SubcycleFineFirstSolver<LevelIntegrator>::GetTotalRefineRatio(
    int fine_level, int coarse_level) const {
  int refine_ratio = 1;
  for (int level = fine_level; level > coarse_level; --level) {
    refine_ratio *= Base::GetRatioToCoarserLevel(level).max();
  }
  return refine_ratio;
}

template <typename LevelIntegrator>
Duration SubcycleFineFirstSolver<LevelIntegrator>::ComputeStableDt() {
  Duration min_dt(std::numeric_limits<double>::max());
  for (int level_num = 0; Base::LevelExists(level_num); ++level_num) {
    int ratio = GetTotalRefineRatio(level_num);
    min_dt = std::min(min_dt, ratio * Base::ComputeStableDt(level_num));
  }
  MPI_Comm comm = Base::GetMpiCommunicator();
  const double local_dt = min_dt.count();
  double global_min_dt{0};
  MPI_Allreduce(&local_dt, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, comm);
  return Duration(global_min_dt);
}

template <typename LevelIntegrator>
Result<void, TimeStepTooLarge>
SubcycleFineFirstSolver<LevelIntegrator>::AdvanceLevel(
    int this_level, Duration dt, std::pair<int, int> subcycle) {
  // PreAdvanceLevel might regrid all finer levels.
  // The Context must make sure that scratch data is allocated
  Base::PreAdvanceLevel(this_level, dt, subcycle);

  auto scale_dt_on_error = [this](Result<void, TimeStepTooLarge> result) {
    TimeStepTooLarge error = result.error();
    int ratio = GetTotalRefineRatio(error.level);
    error.level = 0;
    error.dt *= ratio;
    return error;
  };

  // If a finer level exists in the hierarchy, we subcycle that finer level
  // multiple times and use the fine fluxes on coarse-fine interfaces
  const int next_level = this_level + 1;
  if (Base::LevelExists(next_level)) {
    Base::ResetCoarseFineFluxes(next_level, this_level);
    const int refine_ratio = Base::GetRatioToCoarserLevel(next_level).max();
    for (int r = 0; r < refine_ratio; ++r) {
      auto result =
          AdvanceLevel(next_level, dt / refine_ratio, {r, refine_ratio});
      if (!result) {
        return scale_dt_on_error(result);
      }
    }
  }

  Result<void, TimeStepTooLarge> result =
      Base::AdvanceLevelNonRecursively(this_level, dt, subcycle);
  if (!result) {
    return scale_dt_on_error(result);
  }
  if (Base::LevelExists(next_level)) {
    Base::ApplyFluxCorrection(next_level, this_level, dt);
  }

  // Coarsen inner regions from next finer level to this level.
  if (Base::LevelExists(next_level)) {
    Base::CoarsenConservatively(next_level, this_level);

    // The conservative update and the coarsening happened on conservative
    // variables. We have to reconstruct the missing variables in the complete
    // state.
    Base::CompleteFromCons(this_level, dt);
  }

  Base::CopyScratchToData(this_level);

  // Apply any further context related work after advancing this level.
  // This function can also indicate if some error occured.
  // For example the context could detect unphysical states and return a
  // TooLargeTimeStep error condition.
  result = Base::PostAdvanceLevel(this_level, dt, subcycle);
  if (!result) {
    return scale_dt_on_error(result);
  }
  return result;
}

template <typename LevelIntegrator>
Result<void, TimeStepTooLarge>
SubcycleFineFirstSolver<LevelIntegrator>::AdvanceHierarchy(Duration dt) {
  return AdvanceLevel(0, dt, {0, 1});
}

} // namespace fub

#endif
