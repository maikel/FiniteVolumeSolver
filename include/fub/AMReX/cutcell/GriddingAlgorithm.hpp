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
#include "fub/AMReX/cutcell/PatchHierarchy.hpp"
#include "fub/AnyBoundaryCondition.hpp"
#include "fub/AnyInitialData.hpp"
#include "fub/AnyTaggingMethod.hpp"

#include <AMReX_AmrCore.H>
// #include <AMReX_MultiFabUtil.H>

#include <memory>

namespace fub::amrex::cutcell {
class GriddingAlgorithm;
}

namespace fub {
template <> struct GridTraits<amrex::cutcell::GriddingAlgorithm> {
  using PatchLevel = ::fub::amrex::cutcell::PatchLevel;
  using TagDataHandle = ::amrex::TagBoxArray&;
  using DataReference = ::amrex::MultiFab&;
};
} // namespace fub

namespace fub::amrex::cutcell {
using AnyInitialData = ::fub::AnyInitialData<GriddingAlgorithm>;
using AnyTaggingMethod = ::fub::AnyTaggingMethod<GriddingAlgorithm>;
using AnyBoundaryCondition = ::fub::AnyBoundaryCondition<GriddingAlgorithm>;

/// \ingroup GriddingAlgorithm
///
/// \brief This class modifies and initializes a cutcell::PatchLevel in a
/// cutcell::PatchHierarchy.
class GriddingAlgorithm : private ::amrex::AmrCore {
public:
  /// @{
  /// \name Constructors, Assignments and Destructors

  /// \brief Constructs an empty and invalid GriddingAlgorithm
  GriddingAlgorithm();

  /// \brief Constructs a gridding algorithm and defines all customization
  /// points
  ///
  /// \param hier The base PatchHierarchy that will be modified. It also
  /// contains options that influence the box generation step.
  ///
  /// \param initial_data An InitialCondition policy object
  ///
  /// \param tagging A tagging routine masks cells that need further refinement.
  ///
  /// \param boundary The boundary condition fills the ghost layer for physical
  /// domain boundary.
  ///
  /// \throw Throws `std::bad_alloc` if a memory allocation fails.
  GriddingAlgorithm(PatchHierarchy hier, AnyInitialData initial_data,
                    AnyTaggingMethod tagging, AnyBoundaryCondition boundary);

  /// \brief The copy constructor makes a deep copy of the all data for each MPI
  /// rank.
  GriddingAlgorithm(const GriddingAlgorithm&);

  /// \brief The copy constructor makes a deep copy of the all data for each MPI
  /// rank.
  GriddingAlgorithm& operator=(const GriddingAlgorithm&);

  /// \brief The move constructor moves a gridding algorithm without allocating
  /// any memory.
  GriddingAlgorithm(GriddingAlgorithm&&) noexcept = default;

  /// \brief The move constructor moves a gridding algorithm without allocating
  /// any memory.
  GriddingAlgorithm& operator=(GriddingAlgorithm&&) noexcept = default;

  ~GriddingAlgorithm() noexcept override = default;
  /// @}

  /// @{
  /// \name Accessors

  [[nodiscard]] const PatchHierarchy& GetPatchHierarchy() const noexcept;
  PatchHierarchy& GetPatchHierarchy() noexcept;

  [[nodiscard]] const AnyBoundaryCondition&
  GetBoundaryCondition() const noexcept;

  [[nodiscard]] AnyBoundaryCondition& GetBoundaryCondition() noexcept;

  [[nodiscard]] const AnyInitialData& GetInitialCondition() const noexcept;

  [[nodiscard]] const AnyTaggingMethod& GetTagging() const noexcept;
  /// @}

  /// @{
  /// \name Observers

  /// \brief Returns the number of time steps taken on the coarsest refinement
  /// level.
  [[nodiscard]] std::ptrdiff_t GetCycles() const noexcept {
    return hierarchy_.GetCycles();
  }

  /// \brief Returns the current time point on the coarsest refinement level.
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

  /// \brief Initializes the underlying patch hierarchy using the stored initial
  /// data method.
  ///
  /// \note This function might recreate a refinement level if nesting
  /// conditions are violated upon creating a new refinement level. This implies
  /// that the initial condition might be called multiple times for one
  /// refinement level.
  void InitializeHierarchy(double level_time);
  /// @}

  /// @{
  /// \name Actions

  /// \brief Fill the ghost layer boundary specified of the specifed MultiFab
  /// `mf`.
  void FillMultiFabFromLevel(::amrex::MultiFab& mf, int level_number);
  /// @}

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
  AnyInitialData initial_condition_;
  AnyTaggingMethod tagging_;
  AnyBoundaryCondition boundary_condition_;
};

} // namespace fub::amrex::cutcell

#endif
