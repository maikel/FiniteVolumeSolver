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

#ifndef FUB_AMREX_MultiBlock_INTEGRATOR_CONTEXT2_HPP
#define FUB_AMREX_MultiBlock_INTEGRATOR_CONTEXT2_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/AMReX/cutcell/IntegratorContext.hpp"

#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm2.hpp"

#include <memory>
#include <vector>

namespace fub::amrex {

/// \ingroup IntegratorContext
class MultiBlockIntegratorContext2 {
public:
  using FeedbackFn =
      std::function<void(MultiBlockIntegratorContext2&, Duration)>;

  template <typename TubeEquation, typename PlenumEquation>
  MultiBlockIntegratorContext2(const TubeEquation& tube_equation,
                               const PlenumEquation& plenum_equation,
                               std::vector<IntegratorContext> tubes,
                               std::vector<cutcell::IntegratorContext> plena,
                               std::vector<BlockConnection> connectivity);
  /// @{
  /// \name Member Accessors
  [[nodiscard]] span<IntegratorContext> Tubes() noexcept;
  [[nodiscard]] span<const IntegratorContext> Tubes() const noexcept;

  [[nodiscard]] span<cutcell::IntegratorContext> Plena() noexcept;
  [[nodiscard]] span<const cutcell::IntegratorContext> Plena() const noexcept;

  [[nodiscard]] const std::shared_ptr<MultiBlockGriddingAlgorithm2>&
  GetGriddingAlgorithm() const noexcept;

  const IntegratorContextOptions& GetOptions() const noexcept {
    return plena_[0].GetOptions();
  }

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

  /// \brief Synchronize internal ghost cell data for each inter domain connection.
  void ComputeMultiBlockBoundaryData(int level);

  void CopyDataToScratch(int level);
  void CopyScratchToData(int level);

  /// \brief Replaces the underlying gridding algorithm with the specified one.
  void ResetHierarchyConfiguration(
      std::shared_ptr<MultiBlockGriddingAlgorithm2> gridding);

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

  /// \brief Sets a feedback function that will be called in
  /// PostAdvanceHierarchy.
  ///
  /// This is used to inject user defined behaviour after each time step
  void SetPostAdvanceHierarchyFeedback(FeedbackFn feedback);
  /// @}

  /// @{
  /// \name Member functions relevant for the level integrator algorithm.

  void PreAdvanceHierarchy();

  void PostAdvanceHierarchy(Duration dt);

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
  std::shared_ptr<MultiBlockGriddingAlgorithm2> gridding_;
  FeedbackFn post_advance_hierarchy_feedback_;
};

template <typename TubeEquation, typename PlenumEquation>
MultiBlockIntegratorContext2::MultiBlockIntegratorContext2(
    const TubeEquation& tube_equation, const PlenumEquation& plenum_equation,
    std::vector<IntegratorContext> tubes,
    std::vector<cutcell::IntegratorContext> plena,
    std::vector<BlockConnection> connectivity)
    : tubes_{std::move(tubes)}, plena_{std::move(plena)} {
  std::vector<std::shared_ptr<GriddingAlgorithm>> tube_grids(tubes_.size());
  std::transform(
      tubes_.begin(), tubes_.end(), tube_grids.begin(),
      [](IntegratorContext& ctx) { return ctx.GetGriddingAlgorithm(); });
  std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plenum_grids(
      plena_.size());
  std::transform(plena_.begin(), plena_.end(), plenum_grids.begin(),
                 [](cutcell::IntegratorContext& ctx) {
                   return ctx.GetGriddingAlgorithm();
                 });
  gridding_ = std::make_shared<MultiBlockGriddingAlgorithm2>(
      tube_equation, plenum_equation, std::move(tube_grids),
      std::move(plenum_grids), std::move(connectivity));
  const int nlevel = plena_[0].GetPatchHierarchy().GetNumberOfLevels();
  for (int level = 0; level < nlevel; ++level) {
    ComputeMultiBlockBoundaryData(level);
    ApplyBoundaryCondition(level, Direction::X);
  }
}

} // namespace fub::amrex

#endif
