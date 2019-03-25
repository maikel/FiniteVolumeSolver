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
#include "fub/core/function_ref.hpp"
#include "fub/core/mdspan.hpp"

#include "fub/grid/AMReX/PatchHandle.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"
#include "fub/grid/AMReX/cutcell/GriddingAlgorithm.hpp"

#include <memory>

namespace fub {
namespace amrex {
namespace cutcell {

class HyperbolicSplitIntegratorContext {
public:
  using PatchHandle = ::fub::amrex::PatchHandle;
  static constexpr int Rank = AMREX_SPACEDIM;

  HyperbolicSplitIntegratorContext(std::shared_ptr<GriddingAlgorithm> gridding,
                                   int gcw);

  void ResetHierarchyConfiguration(int level = 0);

  MPI_Comm GetMpiCommunicator() const noexcept;

  double GetDx(PatchHandle patch, Direction dir) const;

  CartesianCoordinates GetCartesianCoordinates(PatchHandle patch) const;

  bool LevelExists(int level) const noexcept;

  int GetRatioToCoarserLevel(int level) const noexcept;

  int GetGhostCellWidth(PatchHandle, Direction);

  void FillGhostLayerTwoLevels(int level, int coarse, Direction direction);

  void FillGhostLayerSingleLevel(int level, Direction direction);

  void AccumulateCoarseFineFluxes(int level, Direction dir, Duration dt);
  void ApplyFluxCorrection(int fine, int coarse, Duration dt, Direction dir);
  void ResetCoarseFineFluxes(int fine, int coarse, Direction dir);
  void CoarsenConservatively(int fine, int coarse, Direction dir);

  const std::shared_ptr<PatchHierarchy>& GetPatchHierarchy() const noexcept;

  ::amrex::MultiFab& GetData(int level);
  PatchDataView<double, Rank + 1> GetData(PatchHandle patch);

  ::amrex::MultiFab& GetScratch(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetScratch(PatchHandle patch, Direction dir);

  ::amrex::MultiFab& GetFluxes(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetFluxes(PatchHandle patch, Direction dir);

  ::amrex::MultiFab& GetBoundaryFluxes(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetBoundaryFluxes(PatchHandle patch,
                                                    Direction dir);

  ::amrex::MultiFab& GetReferenceStates(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetReferenceStates(PatchHandle patch,
                                                     Direction dir);

  ::amrex::MultiFab& GetStabilizedFluxes(int level, Direction dir);
  PatchDataView<double, Rank + 1> GetStabilizedFluxes(PatchHandle patch,
                                                      Direction dir);

  CutCellData<Rank> GetCutCellData(PatchHandle patch, Direction dir);

  Duration GetTimePoint(int level, Direction dir) const;
  void SetTimePoint(Duration t, int level, Direction dir);

  std::ptrdiff_t GetCycles(int level, Direction dir) const;
  void SetCycles(std::ptrdiff_t cycle, int level, Direction dir);

  const ::amrex::Geometry& GetGeometry(int level) const;

  void PreAdvanceLevel(int level_num, Direction dir, Duration dt, int subcycle);

  void PostAdvanceLevel(int level_num, Direction dir, Duration dt,
                        int subcycle);

  void PostProcessNumericFluxes(int level_num, Direction dir, Duration dt);

  void Advance(int level, Direction dir, Duration dt, int subcycle);

  void UpdateConservatively(int level, Direction dir, Duration dt);

  template <typename F> F ForEachPatch(int level, F function) {
    return GetPatchHierarchy()->ForEachPatch(level, function);
  }

  ::amrex::FabType GetCutCellPatchType(PatchHandle handle, int gcw = 0) const;

  BoundaryCondition GetBoundaryCondition(int level) const;

private:
  struct LevelData {
    /// Scratch space with ghost cell widths
    std::array<::amrex::MultiFab, 3> scratch;
    ::amrex::MultiFab reference_states;

    /// These arrays will store the fluxes for each patch level which is present
    /// in the patch hierarchy. These will need to be rebuilt if the
    /// PatchHierarchy changes.
    std::array<::amrex::MultiFab, 3> fluxes;
    std::array<::amrex::MultiFab, 3> stabilized_fluxes;
    std::array<::amrex::MultiFab, 3> boundary_fluxes;

    /// FluxRegister accumulate fluxes on coarse fine interfaces between
    /// refinement level. These will need to be rebuilt whenever the hierarchy
    /// changes.
    ::amrex::FluxRegister coarse_fine;

    std::array<Duration, 3> time_point;
    std::array<Duration, 3> regrid_time_point;
    std::array<std::ptrdiff_t, 3> cycles;
  };

  int ghost_cell_width_;
  std::shared_ptr<cutcell::GriddingAlgorithm> gridding_;
  std::vector<LevelData> data_;
};

template <typename State, typename T, typename Equation>
State MakeViewImpl(const HyperbolicSplitIntegratorContext&,
                   boost::hana::basic_type<State>,
                   mdspan<T, AMREX_SPACEDIM + 1> fab,
                   const Equation& equation) {
  return MakeView(boost::hana::type_c<State>, fab, equation);
}

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
