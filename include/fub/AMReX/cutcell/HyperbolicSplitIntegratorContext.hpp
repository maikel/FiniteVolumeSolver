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

#ifndef FUB_AMREX_CUTCELL_HYPERBOLIC_SPLIT_LEVEL_INTEGRATOR_HPP
#define FUB_AMREX_CUTCELL_HYPERBOLIC_SPLIT_LEVEL_INTEGRATOR_HPP

#include "fub/CartesianCoordinates.hpp"
#include "fub/CutCellData.hpp"
#include "fub/Duration.hpp"
#include "fub/Execution.hpp"
#include "fub/core/function_ref.hpp"
#include "fub/core/mdspan.hpp"

#include "fub/AMReX/PatchHandle.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/cutcell/BoundaryCondition.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"

#include <memory>

namespace fub {
namespace amrex {
namespace cutcell {

class HyperbolicSplitIntegratorContext {
public:
  using PatchHandle = ::fub::amrex::PatchHandle;
  using GriddingAlgorithm = ::fub::amrex::cutcell::GriddingAlgorithm;
  template <typename V> using MakeView = MakeViewImpl<V>;

  static constexpr int Rank = AMREX_SPACEDIM;

  HyperbolicSplitIntegratorContext(std::shared_ptr<GriddingAlgorithm> gridding,
                                   int gcw);

  HyperbolicSplitIntegratorContext(const HyperbolicSplitIntegratorContext&);

  HyperbolicSplitIntegratorContext&
  operator=(const HyperbolicSplitIntegratorContext&);

  HyperbolicSplitIntegratorContext(
      HyperbolicSplitIntegratorContext&&) noexcept = default;

  HyperbolicSplitIntegratorContext&
  operator=(HyperbolicSplitIntegratorContext&&) noexcept = default;

  ~HyperbolicSplitIntegratorContext() = default;

  void ResetHierarchyConfiguration(std::shared_ptr<GriddingAlgorithm> gridding);
  void ResetHierarchyConfiguration(int level = 0);

  MPI_Comm GetMpiCommunicator() const noexcept;

  double GetDx(PatchHandle patch, Direction dir) const;

  CartesianCoordinates GetCartesianCoordinates(PatchHandle patch) const;

  bool LevelExists(int level) const noexcept;

  int GetRatioToCoarserLevel(int level, Direction dir) const noexcept;
  ::amrex::IntVect GetRatioToCoarserLevel(int level) const noexcept;

  int GetGhostCellWidth(PatchHandle, Direction);

  void FillGhostLayer(::amrex::MultiFab& dest, int level);

  void FillGhostLayerTwoLevels(int level, int coarse, Direction direction);

  void FillGhostLayerSingleLevel(int level, Direction direction);

  void AccumulateCoarseFineFluxes(int level, Direction dir, Duration dt);
  void ApplyFluxCorrection(int fine, int coarse, Duration dt, Direction dir);
  void ResetCoarseFineFluxes(int fine, int coarse, Direction dir);
  void CoarsenConservatively(int fine, int coarse, Direction dir);

  const std::shared_ptr<GriddingAlgorithm>& GetGriddingAlgorithm() const
      noexcept;

  PatchHierarchy& GetPatchHierarchy() noexcept;
  const PatchHierarchy& GetPatchHierarchy() const noexcept;

  ::amrex::MultiFab& GetData(int level);
  PatchDataView<double, Rank + 1> GetData(PatchHandle patch);

  ::amrex::MultiFab& GetScratch(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetScratch(PatchHandle patch, Direction dir);

  ::amrex::MultiFab& GetFluxes(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetFluxes(PatchHandle patch, Direction dir);

  ::amrex::MultiCutFab& GetBoundaryFluxes(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetBoundaryFluxes(PatchHandle patch,
                                                    Direction dir);

  ::amrex::MultiFab& GetReferenceStates(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetReferenceStates(PatchHandle patch);

  ::amrex::MultiFab& GetStabilizedFluxes(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetStabilizedFluxes(PatchHandle patch,
                                                      Direction dir);

  ::amrex::MultiFab& GetShieldedLeftFluxes(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetShieldedLeftFluxes(PatchHandle patch,
                                                        Direction dir);

  ::amrex::MultiFab& GetShieldedRightFluxes(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetShieldedRightFluxes(PatchHandle patch,
                                                         Direction dir);

  ::amrex::MultiFab& GetDoublyShieldedFluxes(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetDoublyShieldedFluxes(PatchHandle patch,
                                                          Direction dir);

  CutCellData<Rank> GetCutCellData(PatchHandle patch, Direction dir);

  Duration GetTimePoint(int level, Direction dir) const;
  void SetTimePoint(Duration t, int level, Direction dir);

  std::ptrdiff_t GetCycles(int level, Direction dir) const;
  void SetCycles(std::ptrdiff_t cycle, int level, Direction dir);

  const ::amrex::Geometry& GetGeometry(int level) const;

  void PreAdvanceHierarchy();
  void PostAdvanceHierarchy();

  void PreAdvanceLevel(int level_num, Direction dir, Duration dt, int subcycle);

  void PostAdvanceLevel(int level_num, Direction dir, Duration dt,
                        int subcycle);

  template <typename F> F ForEachPatch(int level, F function) {
    return GetPatchHierarchy().ForEachPatch(level, function);
  }

  template <typename F> double Minimum(int level, F function) {
    return GetPatchHierarchy().Minimum(level, function);
  }

  template <typename F>
  F ForEachPatch(execution::OpenMpTag, int level, F function) {
    return GetPatchHierarchy().ForEachPatch(execution::openmp, level, function);
  }

  template <typename F>
  double Minimum(execution::OpenMpTag, int level, F function) {
    return GetPatchHierarchy().Minimum(execution::openmp, level, function);
  }

  ::amrex::FabType GetCutCellPatchType(PatchHandle handle, int gcw = 0) const;

  BoundaryCondition GetBoundaryCondition(int level) const;

  struct LevelData {
    LevelData() = default;
    LevelData(const LevelData& other) = delete;
    LevelData& operator=(const LevelData& other) = delete;
    LevelData(LevelData&&) noexcept = default;
    LevelData& operator=(LevelData&&) noexcept;
    ~LevelData() noexcept = default;

    /// This eb_factory is shared with the underlying patch hierarchy.
    std::shared_ptr<::amrex::EBFArrayBoxFactory> eb_factory;

    ///////////////////////////////////////////////////////////////////////////
    // [cell-centered]

    /// reference states which are used to compute embedded boundary fluxes
    std::unique_ptr<::amrex::MultiFab> reference_states;

    /// scratch space filled with data in ghost cells
    std::array<::amrex::MultiFab, Rank> scratch;

    /// fluxes for the embedded boundary
    ::amrex::MultiCutFab boundary_fluxes;

    ///////////////////////////////////////////////////////////////////////////
    // [face-centered]

    /// @{
    /// various flux types needed by the numerical scheme
    std::array<::amrex::MultiFab, Rank> fluxes;
    std::array<::amrex::MultiFab, Rank> stabilized_fluxes;
    std::array<::amrex::MultiFab, Rank> shielded_left_fluxes;
    std::array<::amrex::MultiFab, Rank> shielded_right_fluxes;
    std::array<::amrex::MultiFab, Rank> doubly_shielded_fluxes;
    /// @}

    ///////////////////////////////////////////////////////////////////////////
    // [misc]

    /// FluxRegister accumulate fluxes on coarse fine interfaces between
    /// refinement level. These will need to be rebuilt whenever the hierarchy
    /// changes.
    ::amrex::FluxRegister coarse_fine;

    std::array<Duration, Rank> time_point;
    std::array<Duration, Rank> regrid_time_point;
    std::array<std::ptrdiff_t, Rank> cycles;
  };

private:
  int ghost_cell_width_;
  std::shared_ptr<GriddingAlgorithm> gridding_;
  std::vector<LevelData> data_;
};

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
