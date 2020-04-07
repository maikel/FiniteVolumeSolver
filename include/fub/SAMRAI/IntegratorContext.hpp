// Copyright (c) 2019 Maikel Nadolski
// Copyright (c) 2019 Patrick Denzler
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

#ifndef FUB_SAMRAI_INTEGRATOR_CONTEXT_HPP
#define FUB_SAMRAI_INTEGRATOR_CONTEXT_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/HyperbolicMethod.hpp"
#include "fub/TimeStepError.hpp" 
#include "fub/ext/outcome.hpp"
#include "fub/counter/CounterRegistry.hpp"

#include "fub/SAMRAI/GriddingAlgorithm.hpp"

#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>

#include <array>
#include <memory>
#include <vector>

namespace fub {
namespace samrai {

class IntegratorContext;

using HyperbolicMethod = ::fub::HyperbolicMethod<IntegratorContext>;

/// \brief This class is used by the HypebrolicSplitLevelIntegrator and
/// delegates AMR related tasks to the AMReX library.
class IntegratorContext {
public:
  struct AuxialiaryDataDescription {
    std::vector<int> scratch_ids;
    std::vector<int> flux_ids;
    std::vector<int> coarse_fine_ids;
  };

  static AuxialiaryDataDescription
  RegisterVariables(const DataDescription& desc, int scratch_ghost_cell_width,
                    int flux_ghost_cell_width);

  template <typename Method>
  static AuxialiaryDataDescription
  RegisterVariables(const DataDescription& desc, const Method& method) {
    return IntegratorContext::RegisterVariables(desc, method.GetStencilWidth(),
                                                1);
  }

  /// @{
  /// \name Constructors and Assignments

  /// \brief Constructs a context object from given a gridding algorithm and a
  /// numerical method.
  IntegratorContext(std::shared_ptr<GriddingAlgorithm> gridding,
                    HyperbolicMethod method,
                    AuxialiaryDataDescription aux_desc);

  /// \brief Deeply copies a context and all its distributed data for all MPI
  /// ranks.
  IntegratorContext(const IntegratorContext&) = default;

  /// \brief Deeply copies a context and all its distributed data for all MPI
  /// ranks.
  IntegratorContext& operator=(const IntegratorContext&) = default;

  IntegratorContext(IntegratorContext&&) noexcept = default;

  IntegratorContext& operator=(IntegratorContext&&) noexcept = default;

  ~IntegratorContext() = default;
  /// @}

  /// @{
  /// \name Member Accessors

  /// \brief Returns the current boundary condition for the specified level.
  [[nodiscard]] const AnyBoundaryCondition& GetBoundaryCondition() const;
  [[nodiscard]] AnyBoundaryCondition& GetBoundaryCondition();

  /// \brief Returns a shared pointer to the underlying GriddingAlgorithm which
  /// owns the simulation.
  [[nodiscard]] const std::shared_ptr<GriddingAlgorithm>&
  GetGriddingAlgorithm() const noexcept;

  /// \brief Returns a reference to const PatchHierarchy which is a member of
  /// the GriddingAlgorithm.
  [[nodiscard]] const PatchHierarchy& GetPatchHierarchy() const noexcept;

  /// \brief Returns a reference to PatchHierarchy which is a member of the
  /// GriddingAlgorithm.
  [[nodiscard]] PatchHierarchy& GetPatchHierarchy() noexcept;

  /// \brief Returns the MPI communicator which is associated with this context.
  [[nodiscard]] MPI_Comm GetMpiCommunicator() const noexcept;
  /// @}

  /// @{
  /// \name Access Level-specific data

  /// \brief Returns the MultiFab associated with level data on the specifed
  /// level number.
  [[nodiscard]] SAMRAI::hier::PatchLevel& GetPatchLevel(int level);
  [[nodiscard]] const SAMRAI::hier::PatchLevel& GetPatchLevel(int level) const;

  [[nodiscard]] span<const int> GetDataIds() const;
  [[nodiscard]] span<const int> GetScratchIds() const;
  [[nodiscard]] span<const int> GetFluxIds() const;

  /// \brief Returns the current time level for data at the specified refinement
  /// level and direction.
  [[nodiscard]] Duration GetTimePoint(int level = 0) const;

  /// \brief Returns the current number of cycles for data at the specified
  /// refinement level and direction.
  [[nodiscard]] std::ptrdiff_t GetCycles(int level = 0) const;

  /// \brief Returns the geometry object for the specified refinement level.
  [[nodiscard]] const SAMRAI::geom::CartesianGridGeometry&
  GetGeometry(int level) const;
  /// @}

  /// @{
  /// \name Observers

  /// \brief Returns true if the data exists for the specified level number.
  [[nodiscard]] bool LevelExists(int level) const noexcept;

  /// \brief Returns the refinement ratio in the specified direction.
  [[nodiscard]] int GetRatioToCoarserLevel(int level, Direction dir) const
      noexcept;

  /// \brief Returns the refinement ratio for all directions.
  [[nodiscard]] SAMRAI::hier::IntVector GetRatioToCoarserLevel(int level) const
      noexcept;
  /// @}

  /// @{
  /// \name Modifiers

  /// \brief Replaces the underlying gridding algorithm with the specified one.
  void ResetHierarchyConfiguration(std::shared_ptr<GriddingAlgorithm> gridding);

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

  /// \brief On each first subcycle this will regrid the data if neccessary.
  int PreAdvanceLevel(int level_num, Duration dt, std::pair<int, int> subcycle);

  /// \brief Increases the internal time stamps and cycle counters for the
  /// specified level number and direction.
  [[nodiscard]] Result<void, TimeStepTooLarge>
  PostAdvanceLevel(int level_num, Duration dt, int subcycle);

  void ApplyBoundaryCondition(int level, Direction dir);

  /// \brief Fills the ghost layer of the scratch data and interpolates in the
  /// coarse fine layer.
  void FillGhostLayerTwoLevels(int level, int coarse);

  /// \brief Fills the ghost layer of the scratch data and does nothing in the
  /// coarse fine layer.
  void FillGhostLayerSingleLevel(int level);

  /// \brief Returns a estimate for a stable time step size which can be taken
  /// for specified level number in direction dir.
  [[nodiscard]] Duration ComputeStableDt(int level, Direction dir);

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

  void CopyDataToScratch(int level);

  void CopyScratchToData(int level);

  [[nodiscard]] const std::shared_ptr<CounterRegistry>& GetCounterRegistry() const noexcept;


private:
  int ghost_cell_width_;
  std::shared_ptr<GriddingAlgorithm> gridding_;
  HyperbolicMethod method_;
  AuxialiaryDataDescription aux_desc_;

  std::vector<Duration> time_points_;
  std::vector<Duration> regrid_time_points_;
  std::vector<std::ptrdiff_t> cycles_;

  using RefineAlgorithm = SAMRAI::xfer::RefineAlgorithm;
  using RefineSchedule = SAMRAI::xfer::RefineSchedule;
  using CoarsenAlgorithm = SAMRAI::xfer::CoarsenAlgorithm;
  using CoarsenSchedule = SAMRAI::xfer::CoarsenSchedule;

  std::shared_ptr<RefineAlgorithm> fill_scratch_;
  std::shared_ptr<CoarsenAlgorithm> coarsen_fluxes_;
  std::shared_ptr<CoarsenAlgorithm> coarsen_scratch_;

  std::vector<std::shared_ptr<RefineSchedule>> fill_scratch_two_level_schedule_;
  std::vector<std::shared_ptr<RefineSchedule>> fill_scratch_one_level_schedule_;
  std::vector<std::shared_ptr<CoarsenSchedule>> coarsen_scratch_schedule_;
  std::vector<std::shared_ptr<CoarsenSchedule>> coarsen_fluxes_schedule_;

  std::vector<std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy>> boundaries_;
};

} // namespace samrai
} // namespace fub

#endif
