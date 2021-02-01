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

#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm2.hpp"
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/zip.hpp>

namespace fub::amrex {

MultiBlockGriddingAlgorithm2::MultiBlockGriddingAlgorithm2(
    const MultiBlockGriddingAlgorithm2& other)
    : tubes_(other.tubes_.size()), plena_(other.plena_.size()),
      connectivity_(other.connectivity_), boundaries_(other.boundaries_),
      tube_bcs_(other.tube_bcs_), plenum_bcs_(other.plenum_bcs_) {
  for (std::size_t i = 0; i < tubes_.size(); ++i) {
    tubes_[i] = std::make_shared<GriddingAlgorithm>(*other.tubes_[i]);
  }
  for (std::size_t i = 0; i < plena_.size(); ++i) {
    plena_[i] = std::make_shared<cutcell::GriddingAlgorithm>(*other.plena_[i]);
  }
}

MultiBlockGriddingAlgorithm2& MultiBlockGriddingAlgorithm2::
operator=(const MultiBlockGriddingAlgorithm2& other) {
  MultiBlockGriddingAlgorithm2 tmp(other);
  *this = std::move(tmp);
  return *this;
}

void MultiBlockGriddingAlgorithm2::PreAdvanceHierarchy() {
  // Pre-compute ghost cells for each multi block boundary
  for (std::vector<AnyMultiBlockBoundary>& boundaries : boundaries_) {
    for (AnyMultiBlockBoundary& boundary : boundaries) {
      boundary.PreAdvanceHierarchy(*this);
    }
  }
  // Add multi block boundary to the local boundaries of tubes
  for (auto&& [tube_id, original_tube_bc] :
       ranges::view::enumerate(tube_bcs_)) {
    std::vector<AnyBoundaryCondition> all_tube_conditions{original_tube_bc};
    for (auto& boundaries : boundaries_) {
      for (auto&& [conn, bc] : ranges::view::zip(connectivity_, boundaries)) {
        if (conn.tube.id == tube_id) {
          all_tube_conditions.emplace_back(bc);
        }
      }
    }
    amrex::BoundarySet new_tube_bc{};
    new_tube_bc.conditions = all_tube_conditions;
    tubes_[tube_id]->GetBoundaryCondition() = new_tube_bc;
  }
  // Add multi block boundary to the local boundaries of plena
  for (auto&& [plenum_id, original_plenum_bc] :
       ranges::view::enumerate(plenum_bcs_)) {
    std::vector<cutcell::AnyBoundaryCondition> all_plenum_conditions{
        original_plenum_bc};
    for (auto& boundaries : boundaries_) {
      for (auto&& [conn, bc] : ranges::view::zip(connectivity_, boundaries)) {
        if (conn.plenum.id == plenum_id) {
          all_plenum_conditions.emplace_back(bc);
        }
      }
    }
    amrex::cutcell::BoundarySet new_plenum_bc{};
    new_plenum_bc.conditions = all_plenum_conditions;
    plena_[plenum_id]->GetBoundaryCondition() = new_plenum_bc;
  }
}

span<const std::shared_ptr<GriddingAlgorithm>>
MultiBlockGriddingAlgorithm2::GetTubes() const noexcept {
  return tubes_;
}

span<const std::shared_ptr<cutcell::GriddingAlgorithm>>
MultiBlockGriddingAlgorithm2::GetPlena() const noexcept {
  return plena_;
}

span<const BlockConnection>
MultiBlockGriddingAlgorithm2::GetConnectivity() const noexcept {
  return connectivity_;
}

void MultiBlockGriddingAlgorithm2::RegridAllFinerLevels(int which_level) {
  for (const std::shared_ptr<GriddingAlgorithm>& tube : tubes_) {
    tube->RegridAllFinerlevels(which_level);
  }
  for (const std::shared_ptr<cutcell::GriddingAlgorithm>& plenum : plena_) {
    plenum->RegridAllFinerlevels(which_level);
  }
}

} // namespace fub::amrex
