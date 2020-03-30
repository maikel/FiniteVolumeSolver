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

#ifndef FUB_AMREX_GRIDDING_ALGORITHM_HPP
#define FUB_AMREX_GRIDDING_ALGORITHM_HPP

#include "fub/AMReX/BoundaryCondition.hpp"
#include "fub/AMReX/InitialData.hpp"
#include "fub/AMReX/PatchHierarchy.hpp"
#include "fub/AMReX/Tagging.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>

#include <memory>

namespace fub {
namespace amrex {

/// \defgroup GriddingAlgorithm
/// This modules summarizes all GriddingAlgorithms.

/// \ingroup GriddingAlgorithm
/// \brief This class modifies and initializes a PatchLevel in a PatchHierarchy.
class GriddingAlgorithm : private ::amrex::AmrCore {
public:
  //static constexpr int Rank = AMREX_SPACEDIM;

  /// @{
  /// \name Constructors

  /// \brief The copy constructor makes a deep copy of the all data for each MPI rank.
  GriddingAlgorithm(const GriddingAlgorithm& other);

  /// \brief The copy assignment makes a deep copy of the all data for each MPI rank.
  GriddingAlgorithm& operator=(const GriddingAlgorithm& other);

  /// \brief The move constructor moves a gridding algorithm without allocating any memory.
  GriddingAlgorithm(GriddingAlgorithm&& other) noexcept;

  /// \brief The move assignment moves a gridding algorithm without allocating any memory.
  GriddingAlgorithm& operator=(GriddingAlgorithm&& other) noexcept;

  GriddingAlgorithm(PatchHierarchy hier, InitialData initial_data,
                    Tagging tagging);

  GriddingAlgorithm(PatchHierarchy hier, InitialData initial_data,
                    Tagging tagging, AnyBoundaryCondition boundary);

  /// @}

  /// @{
  /// \name Accessors

  PatchHierarchy& GetPatchHierarchy() noexcept { return hierarchy_; }
  const PatchHierarchy& GetPatchHierarchy() const noexcept {
    return hierarchy_;
  }

  [[nodiscard]] const AnyBoundaryCondition& GetBoundaryCondition(int level) const
      noexcept;
  [[nodiscard]] AnyBoundaryCondition& GetBoundaryCondition(int level) noexcept;

  [[nodiscard]] const InitialData& GetInitialCondition() const noexcept;

  [[nodiscard]] const Tagging& GetTagging() const noexcept;
  /// @}


  /// @{
  /// \name Observers
  [[nodiscard]] std::ptrdiff_t GetCycles() const noexcept {
    return hierarchy_.GetCycles();
  }

  [[nodiscard]] Duration GetTimePoint() const noexcept {
    return hierarchy_.GetTimePoint();
  }
  /// @}

  /// @{
  /// \name Modifiers
  /// \brief Attempt to regrid all finer level than the specified `which_level`.
  ///
  /// \return Returns the coarsest level which was regrid. If no level changed
  /// this function returns the maximum number of levels.
  int RegridAllFinerlevels(int which_level);

  void InitializeHierarchy(double level_time);

  void SetBoundaryCondition(int level, const AnyBoundaryCondition& condition);
  void SetBoundaryCondition(int level, AnyBoundaryCondition&& condition);
  /// @}

  /// @{
  /// \name Actions
  void FillMultiFabFromLevel(::amrex::MultiFab& mf, int level_number);
  /// @}

private:
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

  void ClearLevel([[maybe_unused]] int level) override;

  PatchHierarchy hierarchy_;
  InitialData initial_data_;
  Tagging tagging_;
  std::vector<AnyBoundaryCondition> boundary_condition_;
};

} // namespace amrex
} // namespace fub

#endif
