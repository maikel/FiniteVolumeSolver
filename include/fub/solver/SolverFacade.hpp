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

#ifndef FUB_SOLVER_FACADE_HPP
#define FUB_SOLVER_FACADE_HPP

#include <type_traits>

namespace fub {

template <typename LevelIntegrator>
class SolverFacade : private LevelIntegrator {
public:
  template <typename... Args>
  SolverFacade(Args&&... args) : LevelIntegrator(std::forward<Args>(args)...) {}

  [[nodiscard]] LevelIntegrator& GetLevelIntegrator() noexcept { return *this; }
  [[nodiscard]] const LevelIntegrator& GetLevelIntegrator() const noexcept {
    return *this;
  }

  using LevelIntegrator::AdvanceLevelNonRecursively;
  using LevelIntegrator::ApplyFluxCorrection;
  using LevelIntegrator::CoarsenConservatively;
  using LevelIntegrator::CompleteFromCons;
  using LevelIntegrator::ComputeStableDt;
  using LevelIntegrator::CopyDataToScratch;
  using LevelIntegrator::CopyScratchToData;
  using LevelIntegrator::GetContext;
  using LevelIntegrator::GetCycles;
  using LevelIntegrator::GetGriddingAlgorithm;
  using LevelIntegrator::GetMpiCommunicator;
  using LevelIntegrator::GetRatioToCoarserLevel;
  using LevelIntegrator::GetTimePoint;
  using LevelIntegrator::LevelExists;
  using LevelIntegrator::PostAdvanceHierarchy;
  using LevelIntegrator::PostAdvanceLevel;
  using LevelIntegrator::PreAdvanceHierarchy;
  using LevelIntegrator::PreAdvanceLevel;
  using LevelIntegrator::ResetCoarseFineFluxes;
  using LevelIntegrator::ResetHierarchyConfiguration;
};

} // namespace fub

#endif