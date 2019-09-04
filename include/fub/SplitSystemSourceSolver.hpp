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

#ifndef FUB_HYPERBOLIC_SPLIT_SYSTEM_SOURCE_SOLVER_HPP
#define FUB_HYPERBOLIC_SPLIT_SYSTEM_SOURCE_SOLVER_HPP

#include "fub/ext/outcome.hpp"
#include "fub/split_method/StrangSplitting.hpp"

#include <algorithm>

namespace fub {

template <typename Context, typename... Args>
using PreAdvanceHierarchy = decltype(
    std::declval<Context>().PreAdvanceHierarchy(std::declval<Args>()...));

template <typename Context, typename... Args>
using PostAdvanceHierarchy = decltype(
    std::declval<Context>().PostAdvanceHierarchy(std::declval<Args>()...));

////////////////////////////////////////////////////////////////////////////////
//                                                                  Decleration

template <typename SystemSolver, typename SourceTerm,
          typename SplittingMethod = StrangSplitting>
class DimensionalSplitSystemSourceSolver {
public:
  using GriddingAlgorithm = std::decay_t<decltype(
      *std::declval<SystemSolver&>().GetGriddingAlgorithm())>;

  static const int Rank = std::max(SystemSolver::Rank, SourceTerm::Rank);

  /// \brief Constructs a system source solver from given sub solvers.
  DimensionalSplitSystemSourceSolver(SystemSolver system_solver,
                                     SourceTerm source_term,
                                     SplittingMethod split = SplittingMethod());

  // Accessors

  auto& GetContext() noexcept;

  const auto& GetContext() const noexcept;

  /// \brief Returns the shared gridding algorithm of both subsolvers.
  const std::shared_ptr<GriddingAlgorithm>& GetGriddingAlgorithm() const
      noexcept {
    return system_solver_.GetGriddingAlgorithm();
  }

  /// \brief Returns the shared patch hierarchy of both subsolvers.
  const auto& GetPatchHierarchy() const;

  /// \brief Returns the global time point of the current data.
  Duration GetTimePoint() const;

  /// \brief Returns the number of cycles at the coarsest level.
  std::ptrdiff_t GetCycles() const;

  // Modifiers

  /// \brief Resets the hierarchy configuration of both sub solvers.
  template <typename... Args>
  void ResetHierarchyConfiguration(const Args&... args);

  /// \brief Invokes any "pre-advance" logic of both sub solvers.
  void PreAdvanceHierarchy();

  /// \brief Invokes any "post-advance" logic of both sub solvers.
  void PostAdvanceHierarchy();

  /// \brief Returns the minimum of stable time step sizes for both sub solvers.
  Duration ComputeStableDt();

  /// \brief Advances the hierarchy by time step size dt.
  Result<void, TimeStepTooLarge> AdvanceLevel(int level, Duration dt,
                                              int subcycle = 0);
  Result<void, TimeStepTooLarge> AdvanceHierarchy(Duration dt);

private:
  SystemSolver system_solver_;
  SourceTerm source_term_;
  SplittingMethod splitting_;
};

////////////////////////////////////////////////////////////////////////////////
//                                                               Implementation

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
DimensionalSplitSystemSourceSolver<SystemSolver, SourceTerm, SplittingMethod>::
    DimensionalSplitSystemSourceSolver(SystemSolver system_solver,
                                       SourceTerm source_term,
                                       SplittingMethod split)
    : system_solver_{std::move(system_solver)},
      source_term_{std::move(source_term)}, splitting_{split} {}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
template <typename... Args>
void DimensionalSplitSystemSourceSolver<
    SystemSolver, SourceTerm,
    SplittingMethod>::ResetHierarchyConfiguration(const Args&... args) {
  system_solver_.ResetHierarchyConfiguration(args...);
  source_term_.ResetHierarchyConfiguration(args...);
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
void DimensionalSplitSystemSourceSolver<
    SystemSolver, SourceTerm, SplittingMethod>::PreAdvanceHierarchy() {
  if constexpr (is_detected<::fub::PreAdvanceHierarchy, SystemSolver&>()) {
    system_solver_.PreAdvanceHierarchy();
  }
  if constexpr (is_detected<::fub::PreAdvanceHierarchy, SourceTerm&>()) {
    source_term_.PreAdvanceHierarchy();
  }
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
void DimensionalSplitSystemSourceSolver<
    SystemSolver, SourceTerm, SplittingMethod>::PostAdvanceHierarchy() {
  if constexpr (is_detected<::fub::PostAdvanceHierarchy, SourceTerm&>()) {
    source_term_.PostAdvanceHierarchy();
  }
  if constexpr (is_detected<::fub::PostAdvanceHierarchy, SystemSolver&>()) {
    system_solver_.PostAdvanceHierarchy();
  }
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
Duration
DimensionalSplitSystemSourceSolver<SystemSolver, SourceTerm,
                                   SplittingMethod>::ComputeStableDt() {
  return std::min(system_solver_.ComputeStableDt(),
                  source_term_.ComputeStableDt());
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
const auto&
DimensionalSplitSystemSourceSolver<SystemSolver, SourceTerm,
                                   SplittingMethod>::GetPatchHierarchy() const {
  return system_solver_.GetGriddingAlgorithm().GetPatchHierarchy();
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
Duration
DimensionalSplitSystemSourceSolver<SystemSolver, SourceTerm,
                                   SplittingMethod>::GetTimePoint() const {
  return system_solver_.GetTimePoint();
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
std::ptrdiff_t
DimensionalSplitSystemSourceSolver<SystemSolver, SourceTerm,
                                   SplittingMethod>::GetCycles() const {
  return system_solver_.GetCycles();
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
Result<void, TimeStepTooLarge> DimensionalSplitSystemSourceSolver<
    SystemSolver, SourceTerm, SplittingMethod>::AdvanceLevel(int this_level,
                                                             Duration dt,
                                                             int subcycle) {
  // PreAdvanceLevel might regrid this and all finer levels.
  // The Context must make sure that scratch data is allocated
  system_solver_.PreAdvanceLevel(this_level, dt, subcycle);

  // If a finer level exists in the hierarchy, we subcycle that finer level
  // multiple times and use the fine fluxes on coarse-fine interfaces
  const int next_level = this_level + 1;
  if (system_solver_.LevelExists(next_level)) {
    system_solver_.ResetCoarseFineFluxes(next_level, this_level);
    const int refine_ratio =
        system_solver_.GetRatioToCoarserLevel(next_level).max();
    for (int r = 0; r < refine_ratio; ++r) {
      auto result = AdvanceLevel(next_level, dt / refine_ratio, r);
      if (!result) {
        return result.as_failure();
      }
    }
  }

  auto AdvanceSystem = [&](Duration dt) -> Result<void, TimeStepTooLarge> {
    auto AdvanceLevel_Split = [&, dt](Direction dir) {
      return [&, dir, dt, split_cycle = 0](
                 Duration split_dt) mutable -> Result<void, TimeStepTooLarge> {
        if (dir == Direction::X && split_cycle == 0 && subcycle == 0 &&
            this_level > 0) {
          system_solver_.FillGhostLayerTwoLevels(this_level, this_level - 1);
        } else {
          system_solver_.FillGhostLayerSingleLevel(this_level);
        }
        const Duration level_dt =
            system_solver_.GetContext().ComputeStableDt(this_level, dir);
        if (level_dt < split_dt) {
          const int refine_ratio =
              system_solver_.GetTotalRefineRatio(this_level);
          return TimeStepTooLarge{refine_ratio * level_dt};
        }

        // Compute fluxes in the specified direction
        system_solver_.ComputeNumericFluxes(this_level, split_dt, dir);

        if (system_solver_.LevelExists(next_level)) {
          system_solver_.ApplyFluxCorrection(next_level, this_level, split_dt);
        }

        if (this_level > 0) {
          system_solver_.AccumulateCoarseFineFluxes(
              this_level, split_dt.count() / dt.count(), dir);
        }

        // Use the updated fluxes to update cons variables at the "SCRATCH"
        // system_solver_.
        system_solver_.UpdateConservatively(this_level, split_dt, dir);

        // The conservative update and happened on conservative variables.
        // We have to reconstruct the missing variables in the complete state.
        system_solver_.CompleteFromCons(this_level, split_dt);

        split_cycle += 1;
        return boost::outcome_v2::success();
      };
    };

    if (Result<void, TimeStepTooLarge> result = std::apply(
            [&](auto... directions) {
              return system_solver_.GetSplitMethod().Advance(
                  dt, AdvanceLevel_Split(directions)...);
            },
            MakeSplitDirections<Rank>());
        !result) {
      return result;
    }

    // Coarsen inner regions from next finer level to this level.
    if (system_solver_.LevelExists(next_level)) {
      system_solver_.CoarsenConservatively(next_level, this_level);

      // The conservative update and the coarsening happened on conservative
      // variables. We have to reconstruct the missing variables in the complete
      // state.
      system_solver_.CompleteFromCons(this_level, dt);
    }
    return boost::outcome_v2::success();
  };

  auto AdvanceSource = [&](Duration dt) -> Result<void, TimeStepTooLarge> {
    return source_term_.AdvanceLevel(this_level, dt);
  };

  if (Result<void, TimeStepTooLarge> result =
          splitting_.Advance(dt, AdvanceSource, AdvanceSystem);
      !result) {
    return result;
  }

  // Apply any further context related work after advancing this level.
  // This function can also indicate if some error occured.
  // For example the context could detect unphysical states and return a
  // TooLargeTimeStep error condition.
  return system_solver_.PostAdvanceLevel(this_level, dt, subcycle);
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
Result<void, TimeStepTooLarge> DimensionalSplitSystemSourceSolver<
    SystemSolver, SourceTerm, SplittingMethod>::AdvanceHierarchy(Duration dt) {
  return AdvanceLevel(0, dt);
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
auto& DimensionalSplitSystemSourceSolver<
    SystemSolver, SourceTerm, SplittingMethod>::GetContext() noexcept {
  return system_solver_.GetContext();
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
const auto& DimensionalSplitSystemSourceSolver<
    SystemSolver, SourceTerm, SplittingMethod>::GetContext() const noexcept {
  return system_solver_.GetContext();
}

} // namespace fub

#endif
