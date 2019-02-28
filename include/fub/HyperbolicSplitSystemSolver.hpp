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

#ifndef FUB_HYPERBOLIC_SPLIT_SYSTEM_SOLVER_HPP
#define FUB_HYPERBOLIC_SPLIT_SYSTEM_SOLVER_HPP

#include "fub/split_method/GodunovSplitting.hpp"

namespace fub {

template <typename LevelIntegrator, typename SplittingMethod = GodunovSplitting>
struct HyperbolicSplitSystemSolver {
  HyperbolicSplitSystemSolver(LevelIntegrator level_integrator,
                              SplittingMethod split = SplittingMethod())
      : integrator{std::move(level_integrator)}, splitting{split} {}

  void ResetHierarchyConfiguration() {
    integrator.ResetHierarchyConfiguration();
  }

  Duration ComputeStableDt() {
    Duration dt(std::numeric_limits<double>::infinity());
    for (int d = 0; d < LevelIntegrator::Rank(); ++d) {
      Direction dir = static_cast<Direction>(d);
      dt = std::min(dt, integrator.ComputeStableDt(dir));
    }
    return Duration(dt);
  }

  std::array<Direction, LevelIntegrator::Rank()> MakeSplitDirections() {
    std::array<Direction, LevelIntegrator::Rank()> dirs;
    for (int i = 0; i < LevelIntegrator::Rank(); ++i) {
      dirs[i] = Direction(i);
    }
    return dirs;
  }

  const auto& GetPatchHierarchy() const {
    return integrator.GetPatchHierarchy();
  }

  Duration GetTimePoint() const {
    return integrator.GetTimePoint(0, Direction::X);
  }

  std::ptrdiff_t GetCycles() const {
    return integrator.GetCycles(0, Direction::X);
  }

  void AdvanceHierarchy(std::chrono::duration<double> dt) {
    // This transforms a direction into a function which statisfies
    // is_invokable<void, Duration>.
    auto MakeAdvanceFunction = [&](Direction dir) {
      return [&, dir](std::chrono::duration<double> dt) {
        integrator.AdvanceLevel(0, dir, dt);
      };
    };
    // We construct for each split direction a function which will be passed to
    // our splitting method.
    std::apply(
        [&](auto... dir) {
          splitting.Advance(dt, MakeAdvanceFunction(dir)...);
        },
        // TODO: Maybe randomize the directions array?
        MakeSplitDirections());
  }

  LevelIntegrator integrator;
  SplittingMethod splitting;
};

} // namespace fub

#endif