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
class SplitSystemSourceSolver {
public:
  using GriddingAlgorithm = typename SystemSolver::GriddingAlgorithm;
  using Equation = typename SystemSolver::Equation;

  /// \brief Constructs a system source solver from given sub solvers.
  SplitSystemSourceSolver(SystemSolver system_solver, SourceTerm source_term,
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

  /// \brief Returns the underlying equation object of system solver.
  const Equation& GetEquation() const;

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
  Result<void, TimeStepTooLarge> AdvanceHierarchy(Duration dt);

private:
  SystemSolver system_solver_;
  SourceTerm source_term_;
  SplittingMethod splitting_;
};

////////////////////////////////////////////////////////////////////////////////
//                                                               Implementation

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
SplitSystemSourceSolver<SystemSolver, SourceTerm, SplittingMethod>::
    SplitSystemSourceSolver(SystemSolver system_solver, SourceTerm source_term,
                            SplittingMethod split)
    : system_solver_{std::move(system_solver)},
      source_term_{std::move(source_term)}, splitting_{split} {}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
template <typename... Args>
void SplitSystemSourceSolver<SystemSolver, SourceTerm, SplittingMethod>::
    ResetHierarchyConfiguration(const Args&... args) {
  system_solver_.ResetHierarchyConfiguration(args...);
  source_term_.ResetHierarchyConfiguration(args...);
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
void SplitSystemSourceSolver<SystemSolver, SourceTerm,
                             SplittingMethod>::PreAdvanceHierarchy() {
  if constexpr (is_detected<::fub::PreAdvanceHierarchy, SystemSolver&>()) {
    system_solver_.PreAdvanceHierarchy();
  }
  if constexpr (is_detected<::fub::PreAdvanceHierarchy, SourceTerm&>()) {
    source_term_.PreAdvanceHierarchy();
  }
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
void SplitSystemSourceSolver<SystemSolver, SourceTerm,
                             SplittingMethod>::PostAdvanceHierarchy() {
  if constexpr (is_detected<::fub::PostAdvanceHierarchy, SourceTerm&>()) {
    source_term_.PostAdvanceHierarchy();
  }
  if constexpr (is_detected<::fub::PostAdvanceHierarchy, SystemSolver&>()) {
    system_solver_.PostAdvanceHierarchy();
  }
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
Duration SplitSystemSourceSolver<SystemSolver, SourceTerm,
                                 SplittingMethod>::ComputeStableDt() {
  return std::min(system_solver_.ComputeStableDt(),
                  source_term_.ComputeStableDt());
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
const auto&
SplitSystemSourceSolver<SystemSolver, SourceTerm,
                        SplittingMethod>::GetPatchHierarchy() const {
  return system_solver_.GetGriddingAlgorithm().GetPatchHierarchy();
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
Duration SplitSystemSourceSolver<SystemSolver, SourceTerm,
                                 SplittingMethod>::GetTimePoint() const {
  return system_solver_.GetTimePoint();
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
std::ptrdiff_t
SplitSystemSourceSolver<SystemSolver, SourceTerm, SplittingMethod>::GetCycles()
    const {
  return system_solver_.GetCycles();
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
const typename SplitSystemSourceSolver<SystemSolver, SourceTerm,
                                       SplittingMethod>::Equation&
SplitSystemSourceSolver<SystemSolver, SourceTerm,
                        SplittingMethod>::GetEquation() const {
  return system_solver_.GetEquation();
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
Result<void, TimeStepTooLarge>
SplitSystemSourceSolver<SystemSolver, SourceTerm,
                        SplittingMethod>::AdvanceHierarchy(Duration dt) {
  auto system_solver = [&](Duration dt) {
    return system_solver_.AdvanceHierarchy(dt);
  };
  auto source_term = [&](Duration dt) {
    return source_term_.AdvanceHierarchy(dt);
  };
  return splitting_.Advance(dt, source_term, system_solver);
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
  auto& SplitSystemSourceSolver<SystemSolver, SourceTerm,
  SplittingMethod>::GetContext() noexcept {
  return system_solver_.GetContext();
}

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod>
  const auto& SplitSystemSourceSolver<SystemSolver, SourceTerm,
  SplittingMethod>::GetContext() const noexcept {
  return system_solver_.GetContext();
}

} // namespace fub

#endif
