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
#include "fub/split_method/GodunovSplitting.hpp"

#include <stdexcept>

#include <mpi.h>

namespace fub {

template <typename Context, typename... Args>
using PreAdvanceHierarchy = decltype(
    std::declval<Context>().PreAdvanceHierarchy(std::declval<Args>()...));

template <typename Context, typename... Args>
using PostAdvanceHierarchy = decltype(
    std::declval<Context>().PostAdvanceHierarchy(std::declval<Args>()...));

template <typename T>
using GriddingAlgorithm =
    std::decay_t<decltype(*std::declval<T&>().GetGriddingAlgorithm())>;

/// This Level Integrator applies a very general AMR integration scheme in
/// context of dimensional splitting.
///
/// The time integration is split into multiple intermediate steps where each is
/// supposed to do a certain task. The detailed implementation of these tasks
/// happens in the integrator context object.
template <int R, typename Context, typename SplitMethod = GodunovSplitting>
class DimensionalSplitLevelIntegrator : public Context, private SplitMethod {
public:
  static constexpr int Rank = R;

  DimensionalSplitLevelIntegrator(const Context& context,
                                  const SplitMethod& splitting = SplitMethod());

  DimensionalSplitLevelIntegrator(int_constant<R>, const Context& context,
                                  const SplitMethod& splitting = SplitMethod());

  const Context& GetContext() const noexcept;
  Context& GetContext() noexcept;

  const SplitMethod& GetSplitMethod() const noexcept;

  void PreAdvanceHierarchy() {
    if constexpr (is_detected<::fub::PreAdvanceHierarchy, Context&>()) {
      Context::PreAdvanceHierarchy();
    }
  }

  void PostAdvanceHierarchy() {
    if constexpr (is_detected<::fub::PostAdvanceHierarchy, Context&>()) {
      Context::PostAdvanceHierarchy();
    }
  }

  /// Returns the total refinement ratio between specified coarse to fine level
  /// number.
  int GetTotalRefineRatio(int fine_level, int coarse_level = 0) const;

  /// Returns a stable dt across all levels and in one spatial direction.
  ///
  /// This function takes the refinement level into account.
  /// For stability it is advised to use some additional CFL condition.
  Duration ComputeStableDt(int level_number);
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
                                              int subcycle);

  Result<void, TimeStepTooLarge>
  AdvanceLevelNonRecursively(int level_number, Duration dt, int subcycle);

  Result<void, TimeStepTooLarge> AdvanceHierarchy(Duration dt);
};

// Implementation
template <int Rank, typename Context, typename SplitMethod>
DimensionalSplitLevelIntegrator<Rank, Context, SplitMethod>::
    DimensionalSplitLevelIntegrator(const Context& context,
                                    const SplitMethod& splitting)
    : Context(context), SplitMethod(splitting) {}

template <int R, typename Context, typename SplitMethod>
DimensionalSplitLevelIntegrator<R, Context, SplitMethod>::
    DimensionalSplitLevelIntegrator(int_constant<R>, const Context& context,
                                    const SplitMethod& splitting)
    : DimensionalSplitLevelIntegrator(context, splitting) {}

template <int Rank, typename Context, typename SplitMethod>
int DimensionalSplitLevelIntegrator<
    Rank, Context, SplitMethod>::GetTotalRefineRatio(int fine_level,
                                                     int coarse_level) const {
  int refine_ratio = 1;
  for (int level = fine_level; level > coarse_level; --level) {
    refine_ratio *= Context::GetRatioToCoarserLevel(level).max();
  }
  return refine_ratio;
}

template <int R, typename Context, typename SplitMethod>
const Context&
DimensionalSplitLevelIntegrator<R, Context, SplitMethod>::GetContext() const
    noexcept {
  return *this;
}

template <int R, typename Context, typename SplitMethod>
Context& DimensionalSplitLevelIntegrator<R, Context,
                                         SplitMethod>::GetContext() noexcept {
  return *this;
}

template <int R, typename Context, typename SplitMethod>
const SplitMethod&
DimensionalSplitLevelIntegrator<R, Context, SplitMethod>::GetSplitMethod() const
    noexcept {
  return *this;
}

template <int Rank, typename Context, typename SplitMethod>
Duration
DimensionalSplitLevelIntegrator<Rank, Context, SplitMethod>::ComputeStableDt() {
  for (int level_num = 0; Context::LevelExists(level_num); ++level_num) {
    if (level_num > 0) {
      Context::FillGhostLayerTwoLevels(level_num, level_num - 1);
    } else {
      Context::FillGhostLayerSingleLevel(level_num);
    }
  }
  auto ComputeStableDt_Split = [this](Direction dir) -> Duration {
    int refine_ratio = 1;
    double coarse_dt = std::numeric_limits<double>::infinity();
    for (int level_num = 0; Context::LevelExists(level_num); ++level_num) {
      const double level_dt = Context::ComputeStableDt(level_num, dir).count();
      refine_ratio *= Context::GetRatioToCoarserLevel(level_num).max();
      coarse_dt = std::min(refine_ratio * level_dt, coarse_dt);
    }
    return Duration(coarse_dt);
  };
  Duration min_dt(std::numeric_limits<double>::infinity());
  for (int d = 0; d < Rank; ++d) {
    const Direction dir = static_cast<Direction>(d);
    min_dt = std::min(min_dt, ComputeStableDt_Split(dir));
  }
  MPI_Comm comm = Context::GetMpiCommunicator();
  const double local_dt = min_dt.count();
  double global_min_dt{0};
  MPI_Allreduce(&local_dt, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, comm);
  return Duration(global_min_dt);
}

template <int Rank, typename Context, typename SplitMethod>
Result<void, TimeStepTooLarge> DimensionalSplitLevelIntegrator<
    Rank, Context, SplitMethod>::AdvanceLevelNonRecursively(int this_level,
                                                            Duration dt,
                                                            int subcycle) {
  if (subcycle == 0 && this_level > 0) {
    Context::FillGhostLayerTwoLevels(this_level, this_level - 1);
  } else {
    Context::FillGhostLayerSingleLevel(this_level);
  }

  auto AdvanceLevel_Split = [&](Direction dir) {
    return [&, this_level, dir](
               Duration split_dt) -> Result<void, TimeStepTooLarge> {
      const int next_level = this_level + 1;
      const Duration level_dt = Context::ComputeStableDt(this_level, dir);
      if (level_dt < split_dt) {
        const int refine_ratio = GetTotalRefineRatio(this_level);
        return TimeStepTooLarge{refine_ratio * level_dt};
      }

      // Compute fluxes in the specified direction
      Context::ComputeNumericFluxes(this_level, split_dt, dir);

      if (Context::LevelExists(next_level)) {
        Context::ApplyFluxCorrection(next_level, this_level, split_dt);
      }

      if (this_level > 0) {
        Context::AccumulateCoarseFineFluxes(this_level,
                                            split_dt.count() / dt.count(), dir);
      }

      // Use the updated fluxes to update cons variables at the "SCRATCH"
      // context.
      Context::UpdateConservatively(this_level, split_dt, dir);

      // The conservative update and happened on conservative variables.
      // We have to reconstruct the missing variables in the complete state.
      Context::CompleteFromCons(this_level, split_dt);

      return boost::outcome_v2::success();
    };
  };
  return std::apply(
      [&](auto... directions) {
        return GetSplitMethod().Advance(dt, AdvanceLevel_Split(directions)...);
      },
      MakeSplitDirections<Rank>());
}

template <int Rank, typename Context, typename SplitMethod>
Result<void, TimeStepTooLarge>
DimensionalSplitLevelIntegrator<Rank, Context, SplitMethod>::AdvanceLevel(
    int this_level, Duration dt, int subcycle) {
  // PreAdvanceLevel might regrid all finer levels.
  // The Context must make sure that scratch data is allocated
  Context::PreAdvanceLevel(this_level, dt, subcycle);

  // If a finer level exists in the hierarchy, we subcycle that finer level
  // multiple times and use the fine fluxes on coarse-fine interfaces
  const int next_level = this_level + 1;
  if (Context::LevelExists(next_level)) {
    Context::ResetCoarseFineFluxes(next_level, this_level);
    const int refine_ratio = Context::GetRatioToCoarserLevel(next_level).max();
    for (int r = 0; r < refine_ratio; ++r) {
      auto result = AdvanceLevel(next_level, dt / refine_ratio, r);
      if (!result) {
        return result.as_failure();
      }
    }
  }

  Result<void, TimeStepTooLarge> result =
      AdvanceLevelNonRecursively(this_level, dt, subcycle);
  if (!result) {
    return result;
  }

  // Coarsen inner regions from next finer level to this level.
  if (Context::LevelExists(next_level)) {
    Context::CoarsenConservatively(next_level, this_level);

    // The conservative update and the coarsening happened on conservative
    // variables. We have to reconstruct the missing variables in the complete
    // state.
    Context::CompleteFromCons(this_level, dt);
  }

  // Apply any further context related work after advancing this level.
  // This function can also indicate if some error occured.
  // For example the context could detect unphysical states and return a
  // TooLargeTimeStep error condition.
  return Context::PostAdvanceLevel(this_level, dt, subcycle);
}

template <int Rank, typename Context, typename SplitMethod>
Result<void, TimeStepTooLarge>
DimensionalSplitLevelIntegrator<Rank, Context, SplitMethod>::AdvanceHierarchy(
    Duration dt) {
  return AdvanceLevel(0, dt, 0);
}

} // namespace fub

#endif
