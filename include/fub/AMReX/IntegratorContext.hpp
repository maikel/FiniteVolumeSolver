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

#ifndef FUB_AMREX_INTEGRATOR_CONTEXT_HPP
#define FUB_AMREX_INTEGRATOR_CONTEXT_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/HyperbolicMethod.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/ext/outcome.hpp"

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/PatchHierarchy.hpp"

#include <AMReX_FluxRegister.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <memory>
#include <vector>

namespace fub {
namespace amrex {

class IntegratorContext;

using HyperbolicMethod = ::fub::HyperbolicMethod<IntegratorContext>;

/// \brief This class is used by the HypebrolicSplitLevelIntegrator and
/// delegates AMR related tasks to the AMReX library.
class IntegratorContext {
public:
  /// @{
  /// \name Constructors and Assignments

  /// \brief Constructs a context object from given a gridding algorithm and a
  /// numerical method.
  IntegratorContext(std::shared_ptr<GriddingAlgorithm> gridding,
                    HyperbolicMethod method);

  /// \brief Deeply copies a context and all its distributed data for all MPI
  /// ranks.
  IntegratorContext(const IntegratorContext&);

  /// \brief Deeply copies a context and all its distributed data for all MPI
  /// ranks.
  IntegratorContext& operator=(const IntegratorContext&);

  IntegratorContext(IntegratorContext&&) noexcept = default;

  IntegratorContext& operator=(IntegratorContext&&) noexcept = default;

  ~IntegratorContext() = default;
  /// @}

  /// @{
  /// \name Member Accessors

  /// \brief Returns the current boundary condition for the specified level.
  const BoundaryCondition& GetBoundaryCondition(int level) const;
  BoundaryCondition& GetBoundaryCondition(int level);

  /// \brief Returns a shared pointer to the underlying GriddingAlgorithm which
  /// owns the simulation.
  const std::shared_ptr<GriddingAlgorithm>& GetGriddingAlgorithm() const
      noexcept;

  /// \brief Returns a reference to const PatchHierarchy which is a member of
  /// the GriddingAlgorithm.
  const PatchHierarchy& GetPatchHierarchy() const noexcept;

  /// \brief Returns a reference to PatchHierarchy which is a member of the
  /// GriddingAlgorithm.
  PatchHierarchy& GetPatchHierarchy() noexcept;

  /// \brief Returns the MPI communicator which is associated with this context.
  MPI_Comm GetMpiCommunicator() const noexcept;
  /// @}

  /// @{
  /// \name Access Level-specific data

  /// \brief Returns the MultiFab associated with level data on the specifed
  /// level number.
  ::amrex::MultiFab& GetData(int level);

  /// \brief Returns the MultiFab associated with level data with ghost cells on
  /// the specifed level number and direction.
  ::amrex::MultiFab& GetScratch(int level, Direction dir);

  /// \brief Returns the MultiFab associated with flux data on the specifed
  /// level number and direction.
  ::amrex::MultiFab& GetFluxes(int level, Direction dir);

  /// \brief Returns the current time level for data at the specified refinement
  /// level and direction.
  Duration GetTimePoint(int level, Direction dir) const;

  /// \brief Returns the current number of cycles for data at the specified
  /// refinement level and direction.
  std::ptrdiff_t GetCycles(int level, Direction dir) const;

  /// \brief Returns the geometry object for the specified refinement level.
  const ::amrex::Geometry& GetGeometry(int level) const;

  /// @}

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
  void ResetHierarchyConfiguration(std::shared_ptr<GriddingAlgorithm> gridding);

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
  void FillGhostLayerTwoLevels(int level, BoundaryCondition& fbc, int coarse, BoundaryCondition& cbc, Direction direction);
  void FillGhostLayerTwoLevels(int level, int coarse, Direction direction);

  /// \brief Fills the ghost layer of the scratch data and does nothing in the
  /// coarse fine layer.
  void FillGhostLayerSingleLevel(int level, BoundaryCondition& bc, Direction direction);
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

private:
  /// \brief This class holds auxiliary data on each refinement level.
  struct LevelData {
    LevelData() = default;
    LevelData(const LevelData& other) = delete;
    LevelData& operator=(const LevelData& other) = delete;
    LevelData(LevelData&&) noexcept = default;
    LevelData& operator=(LevelData&&) noexcept;
    ~LevelData() noexcept = default;

    /// Scratch space with ghost cell widths
    std::array<::amrex::MultiFab, 3> scratch;

    /// These arrays will store the fluxes for each patch level which is present
    /// in the patch hierarchy. These will need to be rebuilt if the
    /// PatchHierarchy changes.
    std::array<::amrex::MultiFab, 3> fluxes;

    /// FluxRegister accumulate fluxes on coarse fine interfaces between
    /// refinement level. These will need to be rebuilt whenever the hierarchy
    /// changes.
    ::amrex::FluxRegister coarse_fine;

    std::array<Duration, 3> time_point;
    std::array<Duration, 3> regrid_time_point;
    std::array<std::ptrdiff_t, 3> cycles;
  };

  int ghost_cell_width_;
  std::shared_ptr<GriddingAlgorithm> gridding_;
  std::vector<LevelData> data_;
  HyperbolicMethod method_;
};

} // namespace amrex
} // namespace fub

#endif
