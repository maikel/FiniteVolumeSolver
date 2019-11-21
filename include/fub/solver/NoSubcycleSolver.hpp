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

#ifndef FUB_NO_SUBCYCLE_SOLVER_HPP
#define FUB_NO_SUBCYCLE_SOLVER_HPP

#include "fub/Meta.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/ext/outcome.hpp"

namespace fub {

/// This is a solver class which does a AMR integration scheme without
/// subcycling on finer levels. Which means, that this uses one time step size
/// for each refinement level on the patch hierarchy.
template <typename LevelIntegrator>
class NoSubcycleSolver : public SolverFacade<LevelIntegrator> {
private:
  using Base = SolverFacade<LevelIntegrator>;

public:
  explicit NoSubcycleSolver(LevelIntegrator&& integrator)
      : Base(std::move(integrator)) {}

  explicit NoSubcycleSolver(const LevelIntegrator& integrator)
      : Base(integrator) {}

  [[nodiscard]] Duration ComputeStableDt();

  [[nodiscard]] Result<void, TimeStepTooLarge>
  AdvanceLevel(int level, Duration time_step_size);

  [[nodiscard]] Result<void, TimeStepTooLarge>
  AdvanceHierarchy(Duration time_step_size);
};

template <typename LevelIntegrator>
Duration NoSubcycleSolver<LevelIntegrator>::ComputeStableDt() {
  Duration min_dt(std::numeric_limits<double>::max());
  for (int level_num = 0; Base::LevelExists(level_num); ++level_num) {
    min_dt = std::min(min_dt, Base::ComputeStableDt(level_num));
  }
  MPI_Comm comm = Base::GetMpiCommunicator();
  const double local_dt = min_dt.count();
  double global_min_dt{0};
  MPI_Allreduce(&local_dt, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, comm);
  return Duration(global_min_dt);
}

template <typename LevelIntegrator>
Result<void, TimeStepTooLarge>
NoSubcycleSolver<LevelIntegrator>::AdvanceHierarchy(Duration time_step_size) {
  return AdvanceLevel(0, time_step_size);
}

template <typename LevelIntegrator>
Result<void, TimeStepTooLarge>
NoSubcycleSolver<LevelIntegrator>::AdvanceLevel(int level,
                                                Duration time_step_size) {
  Base::PreAdvanceLevel(level, time_step_size, {0,1});
  const int next_level = level + 1;
  if (Base::LevelExists(next_level)) {
    Base::ResetCoarseFineFluxes(next_level, level);
    Result<void, TimeStepTooLarge> result =
        AdvanceLevel(next_level, time_step_size);
    if (!result) {
      return result;
    }
  }
  Result<void, TimeStepTooLarge> result =
      Base::AdvanceLevelNonRecursively(level, time_step_size, {0,1});
  if (!result) {
    return result;
  }
  if (Base::LevelExists(next_level)) {
    Base::ApplyFluxCorrection(next_level, level, time_step_size);
  }
  if (Base::LevelExists(next_level)) {
    Base::CoarsenConservatively(next_level, level);
    Base::CompleteFromCons(level, time_step_size);
  }
  return Base::PostAdvanceLevel(level, time_step_size, {0,1});
}

} // namespace fub

#endif