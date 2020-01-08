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

#ifndef FUB_DIMENSIONAL_SPLIT_LEVEL_INTEGRATOR_HPP
#define FUB_DIMENSIONAL_SPLIT_LEVEL_INTEGRATOR_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/Execution.hpp"
#include "fub/Meta.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/core/algorithm.hpp"
#include "fub/ext/Mpi.hpp"
#include "fub/ext/outcome.hpp"

#include "fub/split_method/GodunovSplitting.hpp"

#include <utility>

namespace fub {

/// This Level Integrator applies a very general AMR integration scheme in
/// context of dimensional splitting.
///
/// The time integration is split into multiple intermediate steps where each is
/// supposed to do a certain task. The detailed implementation of these tasks
/// happens in the integrator context object.
template <int R, typename IntegratorContext,
          typename SplitMethod = GodunovSplitting>
class DimensionalSplitLevelIntegrator : public IntegratorContext,
                                        private SplitMethod {
public:
  static constexpr int Rank = R;

  DimensionalSplitLevelIntegrator(IntegratorContext context,
                                  SplitMethod splitting = SplitMethod())
      : IntegratorContext(std::move(context)),
        SplitMethod(std::move(splitting)) {}

  DimensionalSplitLevelIntegrator(int_constant<R>, IntegratorContext context,
                                  SplitMethod splitting = SplitMethod())
      : DimensionalSplitLevelIntegrator(std::move(context),
                                        std::move(splitting)) {}

  const IntegratorContext& GetContext() const noexcept { return *this; }
  IntegratorContext& GetContext() noexcept { return *this; }

  const SplitMethod& GetSplitMethod() const noexcept { return *this; }

  void PreAdvanceHierarchy();

  void PostAdvanceHierarchy([[maybe_unused]] Duration time_step_size);

  void PreAdvanceLevel([[maybe_unused]] int level,
                       [[maybe_unused]] Duration time_step_size,
                       [[maybe_unused]] std::pair<int, int> subcycle);

  Result<void, TimeStepTooLarge>
  PostAdvanceLevel([[maybe_unused]] int level,
                   [[maybe_unused]] Duration time_step_size,
                   [[maybe_unused]] std::pair<int, int> subcycle);

  using IntegratorContext::GetCycles;
  using IntegratorContext::GetMpiCommunicator;
  using IntegratorContext::GetTimePoint;

  using IntegratorContext::FillGhostLayerSingleLevel;
  using IntegratorContext::FillGhostLayerTwoLevels;

  using IntegratorContext::GetRatioToCoarserLevel;
  using IntegratorContext::LevelExists;

  using IntegratorContext::CoarsenConservatively;
  using IntegratorContext::CompleteFromCons;
  using IntegratorContext::ComputeNumericFluxes;
  using IntegratorContext::CopyDataToScratch;
  using IntegratorContext::CopyScratchToData;
  using IntegratorContext::UpdateConservatively;

  using IntegratorContext::AccumulateCoarseFineFluxes;
  using IntegratorContext::ApplyFluxCorrection;
  using IntegratorContext::ResetCoarseFineFluxes;

  /// Returns a stable dt on a specified level across all spatial directions.
  ///
  /// For stability it is advised to multiply some additional CFL factor < 1.0.
  Duration ComputeStableDt(int level_number);

  /// Advance a specified patch level and all finer levels by time `dt`.
  ///
  /// This method subcycles finer levels.
  ///
  /// \param[in] level_num  An integer denoting the patch level where 0 is the
  /// coarsest level.
  ///
  /// \param[in] dt  A stable time step size for the level_num-th patch level.
  ///
  /// \param[in] subcycle  The ith subcycle which we are currently in, starting
  /// at 0.
  Result<void, TimeStepTooLarge>
  AdvanceLevelNonRecursively(int level_number, Duration dt,
                             std::pair<int, int> subcycle);
};

template <int R, typename IntegratorContext, typename SplitMethod>
void DimensionalSplitLevelIntegrator<R, IntegratorContext,
                                     SplitMethod>::PreAdvanceHierarchy() {
  if constexpr (is_detected<meta::PreAdvanceHierarchy, IntegratorContext&>()) {
    IntegratorContext::PreAdvanceHierarchy();
  }
}
template <int R, typename IntegratorContext, typename SplitMethod>
void DimensionalSplitLevelIntegrator<R, IntegratorContext, SplitMethod>::
    PostAdvanceHierarchy([[maybe_unused]] Duration time_step_size) {
  if constexpr (is_detected<meta::PostAdvanceHierarchy, IntegratorContext&>()) {
    IntegratorContext::PostAdvanceHierarchy();
  } else if constexpr (is_detected<meta::PostAdvanceHierarchy,
                                   IntegratorContext&, Duration>()) {
    IntegratorContext::PostAdvanceHierarchy(time_step_size);
  }
}
template <int R, typename IntegratorContext, typename SplitMethod>
void DimensionalSplitLevelIntegrator<R, IntegratorContext, SplitMethod>::
    PreAdvanceLevel(int level, Duration time_step_size,
                    std::pair<int, int> subcycle) {
  if constexpr (is_detected<meta::PreAdvanceLevel, IntegratorContext&, int,
                            Duration, std::pair<int, int>>()) {
    IntegratorContext::PreAdvanceLevel(level, time_step_size, subcycle);
  }
}
template <int R, typename IntegratorContext, typename SplitMethod>
Result<void, TimeStepTooLarge>
DimensionalSplitLevelIntegrator<R, IntegratorContext, SplitMethod>::
    PostAdvanceLevel(int level, Duration time_step_size,
                     std::pair<int, int> subcycle) {
  if constexpr (is_detected<meta::PostAdvanceLevel, IntegratorContext&, int,
                            Duration, std::pair<int, int>>()) {
    return IntegratorContext::PostAdvanceLevel(level, time_step_size, subcycle);
  }
  return boost::outcome_v2::success();
}

template <int Rank, typename Context, typename SplitMethod>
Duration
DimensionalSplitLevelIntegrator<Rank, Context, SplitMethod>::ComputeStableDt(
    int level) {
  if (level > 0) {
    Context::FillGhostLayerTwoLevels(level, level - 1);
  } else {
    Context::FillGhostLayerSingleLevel(level);
  }
  Duration min_dt(std::numeric_limits<double>::max());
  for (int d = 0; d < Rank; ++d) {
    const Direction dir = static_cast<Direction>(d);
    min_dt = std::min(min_dt, Context::ComputeStableDt(level, dir));
  }
  return min_dt;
}

template <int Rank, typename Context, typename SplitMethod>
Result<void, TimeStepTooLarge>
DimensionalSplitLevelIntegrator<Rank, Context, SplitMethod>::
    AdvanceLevelNonRecursively(int this_level, Duration dt,
                               std::pair<int, int>) {
  auto AdvanceLevel_Split = [&](Direction dir) {
    return [&, this_level, dir, count_split_steps = 0](
               Duration split_dt) mutable -> Result<void, TimeStepTooLarge> {
      if (count_split_steps > 0) {
        // Apply boundary condition for the physical boundary only.
        Context::ApplyBoundaryCondition(this_level, dir);
      }
      count_split_steps += 1;

      // Check stable time step size and if the CFL condition is violated then
      // restart the coarse time step
      const Duration local_dt = Context::ComputeStableDt(this_level, dir);
      MPI_Comm comm = GetMpiCommunicator();
      const Duration level_dt = MinAll(comm, local_dt);
      if (level_dt < split_dt) {
        return TimeStepTooLarge{level_dt, this_level};
      }

      // Compute fluxes in the specified direction
      Context::ComputeNumericFluxes(this_level, split_dt, dir);

      // Use the updated fluxes to update cons variables at the "SCRATCH"
      // context.
      Context::UpdateConservatively(this_level, split_dt, dir);

      // The conservative update and happened on conservative variables.
      // We have to reconstruct the missing variables in the complete state.
      Context::CompleteFromCons(this_level, split_dt);

      // We have to accumulate the fluxes now,
      const double scale = split_dt.count();
      Context::AccumulateCoarseFineFluxes(this_level, scale, dir);

      return boost::outcome_v2::success();
    };
  };

  // Depending on the space dimension of this dimensional split operator we
  // apply the configured splitting method with multiple operators

  Result<void, TimeStepTooLarge> result = std::apply(
      [&](auto... directions) {
        return GetSplitMethod().Advance(dt, AdvanceLevel_Split(directions)...);
      },
      MakeSplitDirections<Rank>());

  if (result) {
    Context::CopyScratchToData(this_level);
  }
  return result;
}

} // namespace fub

#endif