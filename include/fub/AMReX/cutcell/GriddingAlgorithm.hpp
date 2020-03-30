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

#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/cutcell/BoundaryCondition.hpp"
#include "fub/AMReX/cutcell/PatchHierarchy.hpp"
#include "fub/AMReX/cutcell/Tagging.hpp"
#include "fub/AMReX/InitialData.hpp"


#include <AMReX_AmrCore.H>
// #include <AMReX_MultiFabUtil.H>

#include <memory>

namespace fub::amrex::cutcell {

/// \ingroup GriddingAlgorithm
class GriddingAlgorithm : private ::amrex::AmrCore {
public:
  static constexpr int Rank = AMREX_SPACEDIM;

  GriddingAlgorithm() = delete;

  GriddingAlgorithm(const GriddingAlgorithm&);
  GriddingAlgorithm& operator=(const GriddingAlgorithm&);

  GriddingAlgorithm(GriddingAlgorithm&&) noexcept;
  GriddingAlgorithm& operator=(GriddingAlgorithm&&) noexcept;

  ~GriddingAlgorithm() noexcept override = default;

  GriddingAlgorithm(PatchHierarchy hier, AnyInitialData<GriddingAlgorithm> data, Tagging tagging,
                    AnyBoundaryCondition boundary);

  [[nodiscard]] const PatchHierarchy& GetPatchHierarchy() const noexcept;
  PatchHierarchy& GetPatchHierarchy() noexcept;

  int RegridAllFinerlevels(int which_level);
  void InitializeHierarchy(double level_time);

  void SetBoundaryCondition(int level, AnyBoundaryCondition&& condition);
  void SetBoundaryCondition(int level, const AnyBoundaryCondition& condition);

  [[nodiscard]] std::ptrdiff_t GetCycles() const noexcept {
    return hierarchy_.GetCycles();
  }

  [[nodiscard]] Duration GetTimePoint() const noexcept {
    return hierarchy_.GetTimePoint();
  }

  [[nodiscard]] const AnyBoundaryCondition&
  GetBoundaryCondition(int level) const noexcept;

  [[nodiscard]] AnyBoundaryCondition& GetBoundaryCondition(int level) noexcept;

  [[nodiscard]] const AnyInitialData<GriddingAlgorithm>& GetInitialCondition() const noexcept;

  [[nodiscard]] const Tagging& GetTagging() const noexcept;

  void FillMultiFabFromLevel(::amrex::MultiFab& mf, int level_number);

private:
  [[nodiscard]] ::amrex::DistributionMapping
  LoadBalance(int level, const ::amrex::BoxArray& box_array,
              const ::amrex::DistributionMapping& distribution_mapping) const;

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

  void PostProcessBaseGrids(::amrex::BoxArray& box_array) const override;

  void ClearLevel(int level) override;

  PatchHierarchy hierarchy_;
  AnyInitialData<GriddingAlgorithm> initial_condition_;
  Tagging tagging_;
  std::vector<AnyBoundaryCondition> boundary_condition_;
};

} // namespace fub::amrex::cutcell

#endif
