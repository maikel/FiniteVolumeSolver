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

#ifndef FUB_AMREX_CUTCELL_INTEGRATOR_CONTEXT_HPP
#define FUB_AMREX_CUTCELL_INTEGRATOR_CONTEXT_HPP

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/HyperbolicMethod.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/ext/outcome.hpp"

#include <AMReX_FluxRegister.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <memory>
#include <vector>

namespace fub::amrex::cutcell {

class IntegratorContext;

using HyperbolicMethod = ::fub::HyperbolicMethod<IntegratorContext>;

/// \brief This class manages data and the numerical method to perform cut cell
/// simulations with the AMReX library.
class IntegratorContext {
public:
  /// @{
  /// \name Constructors, Assignment Operators and Desctructor

  IntegratorContext(std::shared_ptr<GriddingAlgorithm> gridding,
                    HyperbolicMethod method);

  IntegratorContext(std::shared_ptr<GriddingAlgorithm> gridding,
                    HyperbolicMethod method, int cell_gcw, int face_gcw);

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

  /// \brief Returns the hyperbolic method member object.
  [[nodiscard]] const HyperbolicMethod& GetHyperbolicMethod() const noexcept;

  /// \brief Returns a shared pointer to the underlying GriddingAlgorithm which
  /// owns the simulation.
  [[nodiscard]] const std::shared_ptr<GriddingAlgorithm>&
  GetGriddingAlgorithm() const noexcept;

  /// \brief Returns a shared pointer to the counter registry.
  [[nodiscard]] const std::shared_ptr<CounterRegistry>&
  GetCounterRegistry() const noexcept;

  /// \brief Returns a reference to const PatchHierarchy which is a member of
  /// the GriddingAlgorithm.
  [[nodiscard]] const PatchHierarchy& GetPatchHierarchy() const noexcept;
  [[nodiscard]] PatchHierarchy& GetPatchHierarchy() noexcept;

  /// \brief Returns the MPI communicator which is associated with this context.
  [[nodiscard]] MPI_Comm GetMpiCommunicator() const noexcept;
  /// @}

  /// @{
  /// \name Access Level-specific data

  /// \brief Returns the current boundary condition for the specified level.
  [[nodiscard]] const AnyBoundaryCondition&
  GetBoundaryCondition(int level) const;
  [[nodiscard]] AnyBoundaryCondition& GetBoundaryCondition(int level);

  [[nodiscard]] ::amrex::EBFArrayBoxFactory& GetEmbeddedBoundary(int level);
  [[nodiscard]] const ::amrex::EBFArrayBoxFactory&
  GetEmbeddedBoundary(int level) const;

  /// \brief Returns the MultiFab associated with level data on the specifed
  /// level number.
  [[nodiscard]] ::amrex::MultiFab& GetData(int level);

  /// \brief Returns the MultiFab associated with flux data on the specifed
  /// level number and direction.
  [[nodiscard]] ::amrex::MultiFab& GetReferenceStates(int level);

  /// \brief Returns the MultiFab associated with level data with ghost cells on
  /// the specifed level number and direction.
  [[nodiscard]] ::amrex::MultiFab& GetScratch(int level);

  /// \brief Returns the MultiFab associated with flux data on the specifed
  /// level number and direction.
  [[nodiscard]] ::amrex::MultiCutFab& GetBoundaryFluxes(int level);

  /// \brief Returns the MultiFab associated with flux data on the specifed
  /// level number and direction.
  [[nodiscard]] ::amrex::MultiFab& GetFluxes(int level, Direction dir);

  /// \brief Returns the MultiFab associated with flux data on the specifed
  /// level number and direction.
  [[nodiscard]] ::amrex::MultiFab& GetStabilizedFluxes(int level,
                                                       Direction dir);

  /// \brief Returns the MultiFab associated with flux data on the specifed
  /// level number and direction.
  [[nodiscard]] ::amrex::MultiFab& GetShieldedFromLeftFluxes(int level,
                                                             Direction dir);

  /// \brief Returns the MultiFab associated with flux data on the specifed
  /// level number and direction.
  [[nodiscard]] ::amrex::MultiFab& GetShieldedFromRightFluxes(int level,
                                                              Direction dir);

  /// \brief Returns the MultiFab associated with flux data on the specifed
  /// level number and direction.
  [[nodiscard]] ::amrex::MultiFab& GetDoublyShieldedFluxes(int level,
                                                           Direction dir);

  /// \brief Returns the current time level for data at the specified refinement
  /// level and direction.
  [[nodiscard]] Duration GetTimePoint(int level = 0) const;

  /// \brief Returns the current number of cycles for data at the specified
  /// refinement level and direction.
  [[nodiscard]] std::ptrdiff_t GetCycles(int level = 0) const;

  /// \brief Returns the geometry object for the specified refinement level.
  [[nodiscard]] const ::amrex::Geometry& GetGeometry(int level) const;

  /// @}

  /// @{
  /// \name Observers
  /// \brief Returns true if the data exists for the specified level number.
  [[nodiscard]] bool LevelExists(int level) const noexcept;

  /// \brief Returns the refinement ratio in the specified direction.
  [[nodiscard]] int GetRatioToCoarserLevel(int level, Direction dir) const
      noexcept;

  /// \brief Returns the refinement ratio for all directions.
  [[nodiscard]] ::amrex::IntVect GetRatioToCoarserLevel(int level) const
      noexcept;

  [[nodiscard]] double GetDx(int level, Direction dir) const noexcept;

  [[nodiscard]] ::amrex::FabType GetFabType(int level,
                                            const ::amrex::MFIter& mfi) const
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

  /// \brief Compute and store reference states for cut-cells over all split
  /// cycles.
  ///
  /// \note This assumes that cut cells are always refined to the finest level
  /// of the hierarchy.
  void PreAdvanceHierarchy();
  void PostAdvanceHierarchy();

  /// \brief On each first subcycle this will regrid the data if neccessary.
  int PreAdvanceLevel(int level_num, Duration dt, std::pair<int, int> subcycle);

  /// \brief Increases the internal time stamps and cycle counters for the
  /// specified level number and direction.
  [[nodiscard]] Result<void, TimeStepTooLarge>
  PostAdvanceLevel(int level_num, Duration dt, std::pair<int, int> subcycle);

  void CopyDataToScratch(int level);
  void CopyScratchToData(int level);

  /// \brief Applies the boundary condition for the scratch space on level
  /// `level` in direcition `dir`.
  ///
  /// \param level  The refinement level on which the boundary condition shall
  /// be used.
  void ApplyBoundaryCondition(int level, Direction dir);
  void ApplyBoundaryCondition(int level, Direction dir,
                              AnyBoundaryCondition& bc);

  /// \brief Fills the ghost layer of the scratch data and interpolates in the
  /// coarse fine layer.
  void FillGhostLayerTwoLevels(int level, AnyBoundaryCondition& fbc, int coarse,
                               AnyBoundaryCondition& cbc);
  void FillGhostLayerTwoLevels(int level, int coarse);

  /// \brief Fills the ghost layer of the scratch data and does nothing in the
  /// coarse fine layer.
  void FillGhostLayerSingleLevel(int level, AnyBoundaryCondition& bc);
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

  std::shared_ptr<CounterRegistry> registry_;
  bool count_per_level{false};

private:
  struct LevelData {
    LevelData() = default;
    LevelData(const LevelData& other) = delete;
    LevelData& operator=(const LevelData& other) = delete;
    LevelData(LevelData&&) noexcept = default;
    LevelData& operator=(LevelData&&) noexcept;
    ~LevelData() noexcept = default;

    /// This eb_factory is shared with the underlying patch hierarchy.
    std::shared_ptr<::amrex::EBFArrayBoxFactory> eb_factory{};

    ///////////////////////////////////////////////////////////////////////////
    // [cell-centered]

    /// reference states which are used to compute embedded boundary fluxes
    std::optional<::amrex::MultiFab> reference_states{};

    /// scratch space filled with data in ghost cells
    ::amrex::MultiFab scratch{};

    /// fluxes for the embedded boundary
    std::unique_ptr<::amrex::MultiCutFab> boundary_fluxes{};

    ///////////////////////////////////////////////////////////////////////////
    // [face-centered]

    /// @{
    /// various flux types needed by the numerical scheme
    std::array<::amrex::MultiFab, AMREX_SPACEDIM> fluxes{};
    std::array<::amrex::MultiFab, AMREX_SPACEDIM> stabilized_fluxes{};
    std::array<::amrex::MultiFab, AMREX_SPACEDIM> shielded_left_fluxes{};
    std::array<::amrex::MultiFab, AMREX_SPACEDIM> shielded_right_fluxes{};
    std::array<::amrex::MultiFab, AMREX_SPACEDIM> doubly_shielded_fluxes{};
    /// @}

    ///////////////////////////////////////////////////////////////////////////
    // [misc]

    /// FluxRegister accumulate fluxes on coarse fine interfaces between
    /// refinement level. These will need to be rebuilt whenever the hierarchy
    /// changes.
    ::amrex::FluxRegister coarse_fine{};

    Duration time_point{};
    Duration regrid_time_point{};
    std::ptrdiff_t cycles{};
  };

  int scratch_ghost_cell_width_;
  int flux_ghost_cell_width_{scratch_ghost_cell_width_};
  std::shared_ptr<GriddingAlgorithm> gridding_;
  std::vector<LevelData> data_;
  HyperbolicMethod method_;
};

} // namespace fub::amrex::cutcell

#endif
