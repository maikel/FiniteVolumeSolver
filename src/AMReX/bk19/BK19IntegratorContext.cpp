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

#include "fub/AMReX/bk19/BK19IntegratorContext.hpp"

namespace fub::amrex {

BK19IntegratorContext::BK19IntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, HyperbolicMethod method)
    : IntegratorContext(std::move(gridding), std::move(method)) {
  ResetHierarchyConfiguration();
}

BK19IntegratorContext::BK19IntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, HyperbolicMethod method,
    int cell_gcw, int face_gcw)
    : IntegratorContext(std::move(gridding), std::move(method), cell_gcw,
                        face_gcw) {
  ResetHierarchyConfiguration();
}

/// \brief Deeply copies a context and all its distributed data for all MPI
/// ranks.
BK19IntegratorContext::BK19IntegratorContext(const BK19IntegratorContext& other)
    : IntegratorContext(other) {
  ResetHierarchyConfiguration();
  const int nlevel = GetPatchHierarchy().GetNumberOfLevels();
  const int rank = Rank();
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::Geometry& geom = GetGeometry(level);
    Pv_[level].on_cells.ParallelCopy(other.Pv_[level].on_cells, geom.periodicity());
    for (int d = 0; d < rank; ++d) {
      Pv_[level].on_faces[d].ParallelCopy(other.Pv_[level].on_faces[d], geom.periodicity());
    }
  }
}

BK19AdvectiveFluxes& BK19IntegratorContext::GetAdvectiveFluxes(int level) {
  return Pv_[level];
}

const BK19AdvectiveFluxes&
BK19IntegratorContext::GetAdvectiveFluxes(int level) const {
  return Pv_[level];
}

void BK19IntegratorContext::ResetHierarchyConfiguration(
    std::shared_ptr<GriddingAlgorithm> gridding) {
  Pv_.resize(gridding->GetPatchHierarchy().GetMaxNumberOfLevels());
  IntegratorContext::ResetHierarchyConfiguration(std::move(gridding));
}

void BK19IntegratorContext::ResetHierarchyConfiguration(int coarsest_level) {
  if (Pv_.size() == 0) {
    Pv_.resize(GetPatchHierarchy().GetMaxNumberOfLevels());
  }
  IntegratorContext::ResetHierarchyConfiguration(coarsest_level);
  PatchHierarchy& hier = GetPatchHierarchy();
  const int nlevel = GetPatchHierarchy().GetNumberOfLevels();
  for (int level = coarsest_level; level < nlevel; ++level) {
    const int ngrow_numeric_flux = GetFluxes(level, Direction::X).nGrow();
    const int ngrow_Pv_on_faces = ngrow_numeric_flux + 1;
    const int ngrow_Pv_on_cells = ngrow_Pv_on_faces + 1;
    {
      const ::amrex::BoxArray& ba = hier.GetPatchLevel(level).box_array;
      const ::amrex::DistributionMapping& dm =
          hier.GetPatchLevel(level).distribution_mapping;
      Pv_[level].on_cells = ::amrex::MultiFab(ba, dm, 1, ngrow_Pv_on_cells);
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const ::amrex::BoxArray& ba = GetFluxes(level, Direction(dir)).boxArray();
      const ::amrex::DistributionMapping& dm =
          GetFluxes(level, Direction(dir)).DistributionMap();
      ::amrex::IntVect ngrow_Pv(ngrow_Pv_on_cells);
      ngrow_Pv[dir] = ngrow_Pv_on_faces;
      Pv_[level].on_faces[dir] =
          ::amrex::MultiFab(ba, dm, 1, ngrow_Pv);
    }
  }
}

} // namespace fub::amrex
