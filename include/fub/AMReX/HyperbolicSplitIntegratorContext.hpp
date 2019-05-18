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

#ifndef FUB_AMREX_HYPERBOLIC_SPLIT_INTEGRATOR_CONTEXT_HPP
#define FUB_AMREX_HYPERBOLIC_SPLIT_INTEGRATOR_CONTEXT_HPP

#include "fub/CartesianCoordinates.hpp"
#include "fub/Duration.hpp"
#include "fub/PatchDataView.hpp"
#include "fub/core/function_ref.hpp"
#include "fub/core/mdspan.hpp"

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/PatchHandle.hpp"
#include "fub/AMReX/PatchHierarchy.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"

#include <memory>

namespace fub {
namespace amrex {

/// \brief This class is used by the HypebrolicSplitLevelIntegrator and
/// delegates AMR related tasks to the AMReX library.
class HyperbolicSplitIntegratorContext {
public:
  using PatchHandle = ::fub::amrex::PatchHandle;
  using GriddingAlgorithm = ::fub::amrex::GriddingAlgorithm;
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

  /// \brief Whenever the gridding algorithm changes the data hierarchy this
  /// function will regrid all distributed helper variables managed by the
  /// context.
  ///
  /// \param[in] level  The level number of the coarsest level which changed its
  /// shape. Regrid all levels finer than level.
  void ResetHierarchyConfiguration(int level = 0);

  template <typename F> F ForEachPatch(int level, F function) {
    return GetPatchHierarchy().ForEachPatch(level, function);
  }

  template <typename Feedback> double Minimum(int level, Feedback feedback) {
    return GetPatchHierarchy().Minimum(level, feedback);
  }

  MPI_Comm GetMpiCommunicator() const noexcept;

  double GetDx(PatchHandle patch, Direction dir) const;

  CartesianCoordinates GetCartesianCoordinates(PatchHandle patch) const;

  bool LevelExists(int level) const noexcept {
    return level < GetPatchHierarchy().GetNumberOfLevels();
  }

  int GetRatioToCoarserLevel(int level, Direction dir) const noexcept;
  ::amrex::IntVect GetRatioToCoarserLevel(int level) const noexcept;

  int GetGhostCellWidth(PatchHandle, Direction) { return ghost_cell_width_; }

  void FillGhostLayerTwoLevels(int level, int coarse, Direction direction);

  void FillGhostLayerSingleLevel(int level, Direction direction);

  void AccumulateCoarseFineFluxes(int level, Direction dir, Duration dt);
  void ApplyFluxCorrection(int fine, int coarse, Duration dt, Direction dir);
  void ResetCoarseFineFluxes(int fine, int coarse, Direction dir);
  void CoarsenConservatively(int fine, int coarse, Direction dir);

  const std::shared_ptr<GriddingAlgorithm>& GetGriddingAlgorithm() const
      noexcept;

  const PatchHierarchy& GetPatchHierarchy() const noexcept;
  PatchHierarchy& GetPatchHierarchy() noexcept;

  ::amrex::MultiFab& GetData(int level);
  PatchDataView<double, Rank + 1> GetData(PatchHandle patch);

  ::amrex::MultiFab& GetScratch(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetScratch(PatchHandle patch, Direction dir);

  ::amrex::MultiFab& GetFluxes(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetFluxes(PatchHandle patch, Direction dir);

  Duration GetTimePoint(int level, Direction dir) const;
  void SetTimePoint(Duration t, int level, Direction dir);
  std::ptrdiff_t GetCycles(int level, Direction dir) const;
  void SetCycles(std::ptrdiff_t cycle, int level, Direction dir);
  const ::amrex::Geometry& GetGeometry(int level) const;

  void PreAdvanceLevel(int level_num, Direction dir, Duration dt, int subcycle);

  void PostAdvanceLevel(int level_num, Direction dir, Duration dt,
                        int subcycle);

  BoundaryCondition GetBoundaryCondition(int level) const;

private:
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
};

} // namespace amrex
} // namespace fub

#endif
