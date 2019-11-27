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

#include "fub/Meta.hpp"
#include "fub/ext/outcome.hpp"
#include "fub/split_method/StrangSplitting.hpp"

#include <algorithm>

namespace fub {

////////////////////////////////////////////////////////////////////////////////
//                                                                  Decleration

template <typename SystemSolver, typename SourceTerm,
          typename SplittingMethod = StrangSplitting>
class SplitSystemSourceLevelIntegrator : private SystemSolver {
public:
  using GriddingAlgorithm = meta::GriddingAlgorithm<SystemSolver&>;

  static const int Rank = std::max(SystemSolver::Rank, SourceTerm::Rank);

  /// \brief Constructs a system source solver from given sub solvers.
  SplitSystemSourceLevelIntegrator(SystemSolver system_solver,
                                   SourceTerm source_term,
                                   SplittingMethod split = SplittingMethod());

  // Accessors

  SystemSolver& GetSystem() noexcept { return *this; }
  const SystemSolver& GetSystem() const noexcept { return *this; }

  SourceTerm& GetSource() noexcept { return source_term_; }
  const SourceTerm& GetSource() const noexcept { return source_term_; }

  using SystemSolver::GetContext;
  using SystemSolver::GetCycles;
  using SystemSolver::GetGriddingAlgorithm;
  using SystemSolver::GetMpiCommunicator;
  using SystemSolver::GetTimePoint;

  using SystemSolver::FillGhostLayerSingleLevel;
  using SystemSolver::FillGhostLayerTwoLevels;

  using SystemSolver::GetRatioToCoarserLevel;
  using SystemSolver::LevelExists;

  using SystemSolver::CoarsenConservatively;
  using SystemSolver::CompleteFromCons;
  using SystemSolver::ComputeNumericFluxes;
  using SystemSolver::UpdateConservatively;

  using SystemSolver::AccumulateCoarseFineFluxes;
  using SystemSolver::ApplyFluxCorrection;
  using SystemSolver::ResetCoarseFineFluxes;

  /// \brief Invokes any "pre-advance" logic of both sub solvers.
  void PreAdvanceHierarchy();

  /// \brief Invokes any "post-advance" logic of both sub solvers.
  void PostAdvanceHierarchy(Duration dt);

  void PreAdvanceLevel(int level, Duration dt, std::pair<int, int> subcycle);

  [[nodiscard]] Result<void, TimeStepTooLarge>
  PostAdvanceLevel(int level, Duration dt, std::pair<int, int> subcycle);

  // Modifiers

  /// \brief Resets the hierarchy configuration of both sub solvers.
  template <typename... Args>
  void ResetHierarchyConfiguration(const Args&... args);

  /// \brief Returns the minimum of stable time step sizes for both sub solvers.
  [[nodiscard]] Duration ComputeStableDt(int level);

  [[nodiscard]] Result<void, TimeStepTooLarge>
  AdvanceLevelNonRecursively(int level, Duration dt,
                             std::pair<int, int> subcycle = {0, 1});

private:
  SourceTerm source_term_;
  SplittingMethod splitting_;
};

////////////////////////////////////////////////////////////////////////////////
//                                                               Implementation

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
SplitSystemSourceLevelIntegrator<SystemSolver, SourceTerm, SplittingMethod>::
    SplitSystemSourceLevelIntegrator(SystemSolver system_solver,
                                     SourceTerm source_term,
                                     SplittingMethod split)
    : SystemSolver(std::move(system_solver)),
      source_term_{std::move(source_term)}, splitting_{split} {}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
template <typename... Args>
void SplitSystemSourceLevelIntegrator<
    SystemSolver, SourceTerm,
    SplittingMethod>::ResetHierarchyConfiguration(const Args&... args) {
  SystemSolver::ResetHierarchyConfiguration(args...);
  source_term_.ResetHierarchyConfiguration(args...);
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
void SplitSystemSourceLevelIntegrator<SystemSolver, SourceTerm,
                                      SplittingMethod>::PreAdvanceHierarchy() {
  if constexpr (is_detected<meta::PreAdvanceHierarchy, SystemSolver&>()) {
    SystemSolver::PreAdvanceHierarchy();
  }
  if constexpr (is_detected<meta::PreAdvanceHierarchy, SourceTerm&>()) {
    source_term_.PreAdvanceHierarchy();
  }
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
void SplitSystemSourceLevelIntegrator<
    SystemSolver, SourceTerm,
    SplittingMethod>::PostAdvanceHierarchy([[maybe_unused]] Duration dt) {
  if constexpr (is_detected<meta::PostAdvanceHierarchy, SourceTerm&,
                            Duration>()) {
    source_term_.PostAdvanceHierarchy(dt);
  } else if constexpr (is_detected<meta::PostAdvanceHierarchy, SourceTerm&>()) {
    source_term_.PostAdvanceHierarchy();
  }

  if constexpr (is_detected<meta::PostAdvanceHierarchy, SystemSolver&,
                            Duration>()) {
    SystemSolver::PostAdvanceHierarchy(dt);
  } else if constexpr (is_detected<meta::PostAdvanceHierarchy,
                                   SystemSolver&>()) {
    SystemSolver::PostAdvanceHierarchy();
  }
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
Duration
SplitSystemSourceLevelIntegrator<SystemSolver, SourceTerm,
                                 SplittingMethod>::ComputeStableDt(int level) {
  return std::min(SystemSolver::ComputeStableDt(level),
                  source_term_.ComputeStableDt(level));
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
void SplitSystemSourceLevelIntegrator<
    SystemSolver, SourceTerm,
    SplittingMethod>::PreAdvanceLevel(int this_level, Duration dt,
                                      std::pair<int, int> subcycle) {
  SystemSolver::PreAdvanceLevel(this_level, dt, subcycle);
  if constexpr (is_detected<meta::PreAdvanceLevel, SourceTerm&, int, Duration,
                            std::pair<int, int>>()) {
    source_term_.PreAdvanceLevel(this_level, dt, subcycle);
  }
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
Result<void, TimeStepTooLarge>
SplitSystemSourceLevelIntegrator<SystemSolver, SourceTerm, SplittingMethod>::
    PostAdvanceLevel(int this_level, Duration dt,
                     std::pair<int, int> subcycle) {
  if constexpr (is_detected<meta::PostAdvanceLevel, SourceTerm&, int, Duration,
                            std::pair<int, int>>()) {
    source_term_.PostAdvanceLevel(this_level, dt, subcycle);
  }
  return SystemSolver::PostAdvanceLevel(this_level, dt, subcycle);
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
Result<void, TimeStepTooLarge>
SplitSystemSourceLevelIntegrator<SystemSolver, SourceTerm, SplittingMethod>::
    AdvanceLevelNonRecursively(int this_level, Duration dt,
                               std::pair<int, int> subcycle) {
  auto AdvanceSystem = [&](Duration dt) -> Result<void, TimeStepTooLarge> {
    Result<void, TimeStepTooLarge> result =
        SystemSolver::AdvanceLevelNonRecursively(this_level, dt, subcycle);
    if (!result) {
      return result;
    }
    const int next_level = this_level + 1;
    // Coarsen inner regions from next finer level to this level.
    if (SystemSolver::LevelExists(next_level)) {
      CoarsenConservatively(next_level, this_level);

      // The conservative update and the coarsening happened on conservative
      // variables. We have to reconstruct the missing variables in the complete
      // state.
      CompleteFromCons(this_level, dt);
    }
    return result;
  };

  auto AdvanceSource = [&](Duration dt) -> Result<void, TimeStepTooLarge> {
    return source_term_.AdvanceLevel(this_level, dt);
  };

  return splitting_.Advance(dt, AdvanceSource, AdvanceSystem);
}

} // namespace fub

#endif
