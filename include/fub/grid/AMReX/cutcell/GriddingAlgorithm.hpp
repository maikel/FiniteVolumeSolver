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

#ifndef FUB_AMREX_CUT_CELL_GRIDDING_ALGORITHM_HPP
#define FUB_AMREX_CUT_CELL_GRIDDING_ALGORITHM_HPP

#include "fub/grid/AMReX/BoundaryCondition.hpp"
#include "fub/grid/AMReX/InitialData.hpp"
#include "fub/grid/AMReX/Tagging.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"
#include "fub/grid/AMReX/cutcell/PatchHierarchy.hpp"

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>

#include <memory>

namespace fub {
namespace amrex {
namespace cutcell {

class GriddingAlgorithm : private ::amrex::AmrCore {
public:
  using BoundaryCondition =
      std::function<void(const PatchDataView<double, AMREX_SPACEDIM + 1>&,
                         PatchHandle, Location, int, Duration)>;
  static constexpr int Rank = AMREX_SPACEDIM;

  GriddingAlgorithm(std::shared_ptr<PatchHierarchy> hier, InitialData data,
                    Tagging tagging, BoundaryCondition boundary);

  const std::shared_ptr<PatchHierarchy>& GetPatchHierarchy() const noexcept;

  bool RegridAllFinerlevels(int which_level);
  void InitializeHierarchy(double level_time);

  const BoundaryCondition& GetBoundaryCondition() const noexcept;

private:
  void FillMultiFabFromLevel(::amrex::MultiFab& mf, int level_number);

  void ErrorEst(int level, ::amrex::TagBoxArray& tags, double time_point,
                int /* ngrow */) override;

  void MakeNewLevelFromScratch(
      int level, double time_point, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override;

  void MakeNewLevelFromCoarse(
      int level, double time_point, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override;

  void RemakeLevel(
      int level, double time_point, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override;

  void ClearLevel(int level) override;

  std::shared_ptr<PatchHierarchy> hierarchy_;
  InitialData initial_condition_;
  BoundaryCondition boundary_condition_;
  Tagging tagging_;
};

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
