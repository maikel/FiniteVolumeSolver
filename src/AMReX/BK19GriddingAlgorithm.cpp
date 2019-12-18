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

#include "fub/AMReX/BK19GriddingAlgorithm.hpp"

#include <AMReX_FillPatchUtil.H>

namespace fub::amrex {

void BK19GriddingAlgorithm::AllocatePvs(int level, Direction dir) {
  const int dir_v = static_cast<int>(dir);
  const std::size_t lvl = static_cast<std::size_t>(level);
  FUB_ASSERT(level <= Pvs_[dir_v].size());
  if (Pvs_[dir_v].size() == level) {
    Pvs_[dir_v].emplace_back(box_array, distribution_mapping, 1,
                            pv_ghost_cell_width_);
  } else {
    Pvs_[dir_v][lvl] = ::amrex::MultiFab(box_array, distribution_mapping, 1,
                                         pv_ghost_cell_width_);
  }
}

void BK19GriddingAlgorithm::AllocatePvs(int level) {
  for (int dir_v = 0; dir_v < space_dimension; ++dir_v) {
    AllocatePvs(level, Direction(dir_v));
  }
}

void BK19GriddingAlgorithm::AllocatePvs() {
  std::size_t nlevel =
      static_cast<std::size_t>(GetPatchHierarchy().GetNumberOfLevels());
  const std::size_t sRank = static_cast<std::size_t>(Rank);
  for (std::size_t dir_v = 0; d < sRank; ++dir_v) {
    Pvs_[dir_v].resize(nlevel);
    for (std::size_t level = 0; level < nlevel; ++level) {
      Pvs_[dir_v][level] = ::amrex::MultiFab(box_array, distribution_mapping, 1,
                            pv_ghost_cell_width_);
    }
  }
}

BK19GriddingAlgorithm::BK19GriddingAlgorithm(const BK19GriddingAlgorithm& other)
    : GriddingAlgorithm(other), pv_ghost_cell_width_{
                                    other.pv_ghost_cell_width_} {
  AllocatePvs();
}

BK19GriddingAlgorithm& BK19GriddingAlgorithm::
operator=(const BK19GriddingAlgorithm& other) {
  BK19GriddingAlgorithm tmp{other};
  return *this = std::move(tmp);
}

BK19GriddingAlgorithm::BK19GriddingAlgorithm(
    BK19GriddingAlgorithm&& other) noexcept
    : GriddingAlgorithm(std::move(other)), pv_ghost_cell_width_{
                                               other.pv_ghost_cell_width_} {
  AllocatePvs();
}

BK19GriddingAlgorithm& BK19GriddingAlgorithm::
operator=(BK19GriddingAlgorithm&& other) noexcept {
  *this = static_cast<GriddingAlgorithm&&>(other);
  pv_ghost_cell_width_ = other.pv_ghost_cell_width_;
  Pvs_ = std::move(other.Pvs_);
  return *this;
}

BK19GriddingAlgorithm::BK19GriddingAlgorithm(PatchHierarchy hier,
                                             InitialData initial_data,
                                             Tagging tagging, int pv_gcw)
    : GriddingAlgorithm(std::move(hier), std::move(initial_data),
                        std::move(tagging)),
      pv_ghost_cell_width_{pv_gcw} {
  AllocatePvs();
}

BK19GriddingAlgorithm::BK19GriddingAlgorithm(PatchHierarchy hier,
                                             InitialData initial_data,
                                             Tagging tagging,
                                             BoundaryCondition bc, int pv_gcw)
    : GriddingAlgorithm(std::move(hier), std::move(initial_data),
                        std::move(tagging), std::move(bc)),
      pv_ghost_cell_width_{pv_gcw} {
  AllocatePvs();
}

void BK19GriddingAlgorithm::MakeNewLevelFromScratch(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  // Allocate level data.
  GriddingAlgorithm::MakeNewLevelFromScratch(level, time_point, box_array,
                                             distribution_mapping);
  AllocatePvs(level);
}

void BK19GriddingAlgorithm::MakeNewLevelFromCoarse(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  GriddingAlgorithm::MakeNewLevelFromCoarse(level, time_point, box_array,
                                            distribution_mapping);
  AllocatePvs(level);
}

void BK19GriddingAlgorithm::RemakeLevel(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  GriddingAlgorithm::RemakeLevel(level, time_point, box_array,
                                 distribution_mapping);
  AllocatePvs(level);
}

void BK19GriddingAlgorithm::ClearLevel(int level) {
  GriddingAlgorithm::ClearLevel(level, time_point, box_array,
                                 distribution_mapping);
  Pvs_.pop_back();
}

} // namespace fub::amrex
