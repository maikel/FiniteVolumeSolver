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

#include "fub/split_method/GodunovSplitting.hpp"
#include "fub/ext/outcome.hpp"

namespace fub {

template <typename LevelIntegrator>
using IntegratorContext = std::decay_t<decltype(
    std::declval<LevelIntegrator>().GetIntegratorContext())>;

template <typename Context, typename... Args>
using PreAdvanceHierarchy = decltype(
    std::declval<Context>().PreAdvanceHierarchy(std::declval<Args>()...));

template <typename Context, typename... Args>
using PostAdvanceHierarchy = decltype(
    std::declval<Context>().PostAdvanceHierarchy(std::declval<Args>()...));

template <typename SystemSolver, typename SourceTerm, typename SplittingMethod = GodunovSplitting>
struct HyperbolicSplitSystemSourceSolver {
  HyperbolicSplitSystemSourceSolver(SystemSolver system_solver, SourceTerm source_term,
                              SplittingMethod split = SplittingMethod())
      : system_solver_{std::move(system_solver)}, source_term_{std::move(source_term)}, splitting{split} {}

  void ResetHierarchyConfiguration() {
    system_solver_.ResetHierarchyConfiguration();
  }

  void PreAdvanceHierarchy() {
    if constexpr (is_detected<PreAdvanceHierarchy, SystemSolver&>()) {
      system_solver_.PreAdvanceHierarchy();
    }
    if constexpr (is_detected<PreAdvanceHierarchy, SourceTerm&>()) {
      source_term_.PreAdvanceHierarchy();
    }
  }

  void PostAdvanceHierarchy() {
    if constexpr (is_detected<PostAdvanceHierarchy, Context&>()) {
      system_solver_.PostAdvanceHierarchy();
    }
    if constexpr (is_detected<PostAdvanceHierarchy, Context&>()) {
      source_term_.PostAdvanceHierarchy();
    }
  }

  Duration ComputeStableDt() {
    return std::min(system_solver_.ComputeStableDt(), source_term_.ComputeStableDt());
  }

  const auto& GetPatchHierarchy() const {
    return system_solver_.GetPatchHierarchy();
  }

  Duration GetTimePoint() const {
    return system_solver_.GetTimePoint();
  }

  std::ptrdiff_t GetCycles() const {
    return system_solver_.GetCycles();
  }

  Result<void, TimeStepTooLarge>
  AdvanceHierarchy(std::chrono::duration<double> dt) {
    auto system_solver = [&](std::chrono::duration<double> dt) {
      return system_solver_.AdvanceLevel(dt);
    };
    auto source_term = [&](std::chrono::duration<double> dt) {
      return source_term_.AdvanceLevel(dt);
    };
    return splitting.Advance(dt, system_solver, soruce_term);
  }

  SystemSolver system_solver_;
  SourceTerm source_term_;
  SplittingMethod splitting_;
};

} // namespace fub

#endif
