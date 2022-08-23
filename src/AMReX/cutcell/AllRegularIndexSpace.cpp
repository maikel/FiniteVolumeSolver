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

#include "fub/AMReX/cutcell/AllRegularIndexSpace.hpp"

#include "fub/core/assert.hpp"

namespace fub::amrex::cutcell {

AllRegularIndexSpace::AllRegularIndexSpace(const ::amrex::Geometry& geom,
                                           const ::amrex::IntVect& refine_ratio,
                                           int required_coarsening_level,
                                           int max_coarsening_level,
                                           int ngrow) {
  max_coarsening_level =
      std::max(required_coarsening_level, max_coarsening_level);
  max_coarsening_level = std::min(30, max_coarsening_level);

  int ngrow_finest = std::max(ngrow, 0);
  for (int i = 1; i <= required_coarsening_level; ++i) {
    ngrow_finest *= 2;
  }
  auto gshop = ::amrex::EB2::makeShop(::amrex::EB2::AllRegularIF());
  geometries_.push_back(geom);
  domains_.push_back(geom.Domain());
  n_grows_.push_back(ngrow_finest);
  levels_.reserve(static_cast<std::size_t>(max_coarsening_level + 1));
  levels_.emplace_back(this, gshop, geom, ::amrex::EB2::max_grid_size,
                       ngrow_finest, false, 0);

  for (int ilev = 1; ilev <= max_coarsening_level; ++ilev) {
    bool coarsenable = geometries_.back().Domain().coarsenable(refine_ratio);
    if (!coarsenable) {
      if (ilev <= required_coarsening_level) {
        ::amrex::Abort("IndexSpaceImp: domain is not coarsenable at level " +
                       std::to_string(ilev));
      } else {
        break;
      }
    }

    int ng = (ilev > required_coarsening_level) ? 0 : n_grows_.back() / 2;

    ::amrex::Box cdomain =
        ::amrex::coarsen(geometries_.back().Domain(), refine_ratio);
    ::amrex::Geometry cgeom =
        ::amrex::coarsen(geometries_.back(), refine_ratio);
    levels_.emplace_back(this, ilev, ::amrex::EB2::max_grid_size, ng, cgeom,
                         levels_[static_cast<std::size_t>(ilev - 1)]);
    if (!levels_.back().isOK()) {
      levels_.pop_back();
      if (ilev <= required_coarsening_level) {
        ::amrex::Abort("Failed to build required coarse EB level " +
                       std::to_string(ilev));
      } else {
        break;
      }
    }
    geometries_.push_back(cgeom);
    domains_.push_back(cdomain);
    n_grows_.push_back(ng);
  }
}

void AllRegularIndexSpace::addFineLevels(int) {}

const ::amrex::EB2::Level&
AllRegularIndexSpace::getLevel(const ::amrex::Geometry& geom) const {
  auto it = std::find(domains_.begin(), domains_.end(), geom.Domain());
  const auto i = static_cast<std::size_t>(std::distance(domains_.begin(), it));
  return levels_[i];
}

const ::amrex::Geometry&
AllRegularIndexSpace::getGeometry(const ::amrex::Box& dom) const {
  auto it = std::find(std::begin(domains_), std::end(domains_), dom);
  const auto i = static_cast<std::size_t>(std::distance(domains_.begin(), it));
  return geometries_[i];
}

const ::amrex::Box& AllRegularIndexSpace::coarsestDomain() const {
  return domains_.back();
}

std::vector<const ::amrex::EB2::IndexSpace*>
BuildRegularSpace(const ::amrex::Geometry& coarse_geom,
                  const ::amrex::IntVect& refine_ratio, int n_required_levels) {
  FUB_ASSERT(n_required_levels > 0);
  std::vector<const ::amrex::EB2::IndexSpace*> index_spaces(
      static_cast<std::size_t>(n_required_levels));
  ::amrex::Geometry geom = coarse_geom;
  for (int level = 0; level < n_required_levels; ++level) {
    std::unique_ptr index_space = std::make_unique<AllRegularIndexSpace>(
        geom, refine_ratio, level, level + 1);
    ::amrex::EB2::IndexSpace::push(index_space.release());
    index_spaces[static_cast<std::size_t>(level)] =
        &::amrex::EB2::IndexSpace::top();
    geom.refine(refine_ratio);
  }
  return index_spaces;
}

} // namespace fub::amrex::cutcell