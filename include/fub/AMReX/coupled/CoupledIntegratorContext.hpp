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

#ifndef FUB_AMREX_COUPLED_INTEGRATOR_CONTEXT_HPP
#define FUB_AMREX_COUPLED_INTEGRATOR_CONTEXT_HPP

#include "fub/AMReX/HyperbolicSplitIntegratorContext.hpp"
#include "fub/AMReX/cutcell/HyperbolicSplitIntegratorContext.hpp"

#include <vector>

namespace fub::amrex {

class CoupledIntegratorContext;
using HyperbolicMethod = ::fub::HyperbolicMethod<CoupledIntegratorContext>;

class CoupledIntegratorContext {
private:
  std::vector<HyperbolicSplitIntegratorContext> tubes_;
  std::vector<cutcell::HyperbolicSplitIntegratorContext> plena_;
  std::shared_ptr<CoupledGriddingAlgorithm> gridding_;
  HyperbolicMethod method_;

public:
  CoupledIntegratorContext(std::shared_ptr<CoupledGriddingAlgorithm> gridding,
                    HyperbolicMethod method) = default;

  span<HyperbolicSplitIntegratorContext> Tubes() noexcept;
  span<const HyperbolicSplitIntegratorContext> Tubes() const noexcept;

  span<cutcell::HyperbolicSplitIntegratorContext> Plena() noexcept;
  span<const cutcell::HyperbolicSplitIntegratorContext> Plena() const noexcept;

  const std::shared_ptr<CoupledGriddingAlgorithm>& GetGriddingAlgorithm() const noexcept;

  /// \brief Returns the current boundary condition for the specified level.
  CoupledBoundaryCondition& GetBoundaryCondition(int level);

  /// @{
  /// \name Observers

  /// \brief Returns true if the data exists for the specified level number.
  bool LevelExists(int level) const noexcept;

  /// \brief Returns the refinement ratio in the specified direction.
  int GetRatioToCoarserLevel(int level, Direction dir) const noexcept;

  /// \brief Returns the refinement ratio for all directions.
  ::amrex::IntVect GetRatioToCoarserLevel(int level) const noexcept;
  /// @}

  /// @{
  /// \name Modifiers

  /// \brief Replaces the underlying gridding algorithm with the specified one.
  void ResetHierarchyConfiguration(std::shared_ptr<CoupledGriddingAlgorithm> gridding);

  /// \brief Whenever the gridding algorithm changes the data hierarchy this
  /// function will regrid all distributed helper variables managed by the
  /// context.
  ///
  /// \param[in] level  The level number of the coarsest level which changed its
  /// shape. Regrid all levels finer than level.
  void ResetHierarchyConfiguration(int level = 0);

  /// \brief Sets the cycle count for a specific level number and direction.
  void SetCycles(std::ptrdiff_t cycle, int level, Direction dir);

  /// \brief Sets the time point for a specific level number and direction.
  void SetTimePoint(Duration t, int level, Direction dir);
  /// @}

  /// @{
  /// \name Member functions relevant for the level integrator algorithm.

  /// \brief On each first subcycle this will regrid the data if neccessary.
  void PreAdvanceLevel(int level_num, Direction dir, Duration dt, int subcycle);

  /// \brief Increases the internal time stamps and cycle counters for the
  /// specified level number and direction.
  Result<void, TimeStepTooLarge> PostAdvanceLevel(int level_num, Direction dir,
                                                  Duration dt, int subcycle);

  /// \brief Fills the ghost layer of the scratch data and interpolates in the
  /// coarse fine layer.
  void FillGhostLayerTwoLevels(int level, int coarse, Direction direction);

  /// \brief Fills the ghost layer of the scratch data and does nothing in the
  /// coarse fine layer.
  void FillGhostLayerSingleLevel(int level, Direction direction);

  /// \brief Returns a estimate for a stable time step size which can be taken
  /// for specified level number in direction dir.
  Duration ComputeStableDt(int level, Direction dir);

  /// \brief Fill the flux MultiFab with numeric fluxes based on current states
  /// in scratch.
  void ComputeNumericFluxes(int level, Duration dt, Direction dir);

  /// \brief Apply a conservative time update for each conservative variable on
  /// the specified level number and direction.
  void UpdateConservatively(int level, Duration dt, Direction dir);

  /// \brief Reconstruct complete state variables from conservative ones.
  void CompleteFromCons(int level, Duration dt, Direction dir);

  /// \brief Accumulate fluxes on the coarse fine interfaces for a specified
  /// fine level number.
  void AccumulateCoarseFineFluxes(int level, Duration dt, Direction dir);

  /// \brief Replace the coarse fluxes by accumulated fine fluxes on the coarse
  /// fine interfaces.
  void ApplyFluxCorrection(int fine, int coarse, Duration dt, Direction dir);

  /// \brief Resets all accumulates fluxes to zero.
  void ResetCoarseFineFluxes(int fine, int coarse, Direction dir);

  /// \brief Coarsen scratch data from a fine level number to a coarse level
  /// number.
  void CoarsenConservatively(int fine, int coarse, Direction dir);
  ///@}
};

} // namespace fub::amrex