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

#ifndef FUB_AMREX_MultiBlock_INTEGRATOR_CONTEXT_HPP
#define FUB_AMREX_MultiBlock_INTEGRATOR_CONTEXT_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/AMReX/cutcell/IntegratorContext.hpp"

#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp"

#include <memory>
#include <vector>

namespace fub::amrex {

/// \ingroup IntegratorContext
class MultiBlockIntegratorContext {
public:
  MultiBlockIntegratorContext(FlameMasterReactor reactor,
                              std::vector<IntegratorContext> tubes,
                              std::vector<cutcell::IntegratorContext> plena,
                              std::vector<BlockConnection> connectivity);

  MultiBlockIntegratorContext(
      FlameMasterReactor reactor, std::vector<IntegratorContext> tubes,
      std::vector<cutcell::IntegratorContext> plena,
      std::vector<BlockConnection> connectivity,
      std::vector<std::shared_ptr<PressureValve>> valves);

  /// @{
  /// \name Member Accessors
  [[nodiscard]] span<IntegratorContext> Tubes() noexcept;
  [[nodiscard]] span<const IntegratorContext> Tubes() const noexcept;

  [[nodiscard]] span<cutcell::IntegratorContext> Plena() noexcept;
  [[nodiscard]] span<const cutcell::IntegratorContext> Plena() const noexcept;

  [[nodiscard]] const std::shared_ptr<MultiBlockGriddingAlgorithm>&
  GetGriddingAlgorithm() const noexcept;

  /// \brief Returns the current time level for data at the specified refinement
  /// level and direction.
  [[nodiscard]] Duration GetTimePoint(int level = 0) const;

  /// \brief Returns the current number of cycles for data at the specified
  /// refinement level and direction.
  [[nodiscard]] std::ptrdiff_t GetCycles(int level = 0) const;
  /// @}

  /// \brief Returns the current boundary condition for the specified level.
  [[nodiscard]] MultiBlockBoundary& GetBoundaryCondition(int level);

  [[nodiscard]] MPI_Comm GetMpiCommunicator() const noexcept;

  /// \brief Returns a shared pointer to the counter registry.
  [[nodiscard]] const std::shared_ptr<CounterRegistry>&
  GetCounterRegistry() const noexcept;

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

  void CopyDataToScratch(int level);
  void CopyScratchToData(int level);

  /// \brief Replaces the underlying gridding algorithm with the specified one.
  void ResetHierarchyConfiguration(
      std::shared_ptr<MultiBlockGriddingAlgorithm> gridding);

  /// \brief Whenever the gridding algorithm changes the data hierarchy this
  /// function will regrid all distributed helper variables managed by the
  /// context.
  ///
  /// \param[in] level  The level number of the coarsest level which changed its
  /// shape. Regrid all levels finer than level.
  void ResetHierarchyConfiguration(int level = 0);

  /// \brief Sets the cycle count for a specific level number and direction.
  void SetCycles(std::ptrdiff_t cycle, int level);

  /// \brief Sets the time point for a specific level number and direction.
  void SetTimePoint(Duration t, int level);
  /// @}

  /// @{
  /// \name Member functions relevant for the level integrator algorithm.

  void PreAdvanceHierarchy();

  void PostAdvanceHierarchy();

  /// \brief On each first subcycle this will regrid the data if neccessary.
  int PreAdvanceLevel(int level_num, Duration dt, std::pair<int, int> subcycle);

  /// \brief Increases the internal time stamps and cycle counters for the
  /// specified level number and direction.
  Result<void, TimeStepTooLarge> PostAdvanceLevel(int level_num, Duration dt,
                                                  std::pair<int, int> subcycle);

  /// \brief Applies the boundary condition for the scratch space on level
  /// `level` in direcition `dir`.
  ///
  /// \param level  The refinement level on which the boundary condition shall
  /// be used.
  /// \param dir The dimensional split direction which will be used.
  void ApplyBoundaryCondition(int level, Direction dir);

  /// \brief Fills the ghost layer of the scratch data and interpolates in the
  /// coarse fine layer.
  void FillGhostLayerTwoLevels(int level, int coarse);

  /// \brief Fills the ghost layer of the scratch data and does nothing in the
  /// coarse fine layer.
  void FillGhostLayerSingleLevel(int level);

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
  void CompleteFromCons(int level, Duration dt);

  /// \brief Accumulate fluxes on the coarse fine interfaces for a specified
  /// fine level number.
  void AccumulateCoarseFineFluxes(int level, double time_scale, Direction dir);

  /// \brief Replace the coarse fluxes by accumulated fine fluxes on the coarse
  /// fine interfaces.
  void ApplyFluxCorrection(int fine, int coarse, Duration dt);

  /// \brief Resets all accumulates fluxes to zero.
  void ResetCoarseFineFluxes(int fine, int coarse);

  /// \brief Coarsen scratch data from a fine level number to a coarse level
  /// number.
  void CoarsenConservatively(int fine, int coarse);
  ///@}

private:
  std::vector<IntegratorContext> tubes_;
  std::vector<cutcell::IntegratorContext> plena_;
  std::shared_ptr<MultiBlockGriddingAlgorithm> gridding_;
};

} // namespace fub::amrex

#endif
