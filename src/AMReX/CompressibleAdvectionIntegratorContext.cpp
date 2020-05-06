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

#include "fub/AMReX/CompressibleAdvectionIntegratorContext.hpp"

namespace fub::amrex {

CompressibleAdvectionIntegratorContext::CompressibleAdvectionIntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, HyperbolicMethod method)
    : IntegratorContext(std::move(gridding), std::move(method)) {
  ResetHierarchyConfiguration();
}

CompressibleAdvectionIntegratorContext::CompressibleAdvectionIntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, HyperbolicMethod method,
    int cell_gcw, int face_gcw)
    : IntegratorContext(std::move(gridding), std::move(method), cell_gcw,
                        face_gcw) {
  ResetHierarchyConfiguration();
}

/// \brief Deeply copies a context and all its distributed data for all MPI
/// ranks.
CompressibleAdvectionIntegratorContext::CompressibleAdvectionIntegratorContext(
    const CompressibleAdvectionIntegratorContext& other)
    : IntegratorContext(other) {
  ResetHierarchyConfiguration();
  const int nlevel = GetPatchHierarchy().GetNumberOfLevels();
  const std::size_t rank = static_cast<std::size_t>(Rank());
  for (int level = 0; level < nlevel; ++level) {
    const std::size_t levelli = static_cast<std::size_t>(level);
    const ::amrex::Geometry& geom = GetGeometry(level);
    Pv_[levelli].on_cells.ParallelCopy(other.Pv_[levelli].on_cells,
                                     geom.periodicity());
    for (std::size_t d = 0; d < rank; ++d) {
      Pv_[levelli].on_faces[d].ParallelCopy(other.Pv_[levelli].on_faces[d],
                                          geom.periodicity());
    }
  }
}

CompressibleAdvectionAdvectiveFluxes&
CompressibleAdvectionIntegratorContext::GetAdvectiveFluxes(int level) {
  return Pv_[static_cast<std::size_t>(level)];
}

const CompressibleAdvectionAdvectiveFluxes&
CompressibleAdvectionIntegratorContext::GetAdvectiveFluxes(int level) const {
  return Pv_[static_cast<std::size_t>(level)];
}

void CompressibleAdvectionIntegratorContext::ResetHierarchyConfiguration(
    std::shared_ptr<GriddingAlgorithm> gridding) {
  Pv_.resize(static_cast<std::size_t>(gridding->GetPatchHierarchy().GetMaxNumberOfLevels()));
  IntegratorContext::ResetHierarchyConfiguration(std::move(gridding));
  PatchHierarchy& hier = GetPatchHierarchy();
  const int nlevel = GetPatchHierarchy().GetNumberOfLevels();
  for (int level = 0; level < nlevel; ++level) {
    const std::size_t levelli = static_cast<std::size_t>(level);
    const int ngrow_Pv_on_cells = GetScratch(level).nGrow() + 1;
    {
      const ::amrex::BoxArray& ba = hier.GetPatchLevel(level).box_array;
      const ::amrex::DistributionMapping& dm =
          hier.GetPatchLevel(level).distribution_mapping;
      Pv_[levelli].on_cells = ::amrex::MultiFab(ba, dm, 2, ngrow_Pv_on_cells);
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const std::size_t dirli = static_cast<std::size_t>(dir);
      ::amrex::IntVect ngrow_Pv_on_faces =
          GetFluxes(level, Direction(dir)).nGrowVect();
      ngrow_Pv_on_faces[dir] += 1;
      const ::amrex::BoxArray& ba = GetFluxes(level, Direction(dir)).boxArray();
      const ::amrex::DistributionMapping& dm =
          GetFluxes(level, Direction(dir)).DistributionMap();
      Pv_[levelli].on_faces[dirli] =
          ::amrex::MultiFab(ba, dm, 1, ngrow_Pv_on_faces);
    }
  }
}

void CompressibleAdvectionIntegratorContext::ResetHierarchyConfiguration(
    int coarsest_level) {
  if (Pv_.size() == 0) {
    Pv_.resize(static_cast<std::size_t>(GetPatchHierarchy().GetMaxNumberOfLevels()));
  }
  IntegratorContext::ResetHierarchyConfiguration(coarsest_level);
  PatchHierarchy& hier = GetPatchHierarchy();
  const int nlevel = GetPatchHierarchy().GetNumberOfLevels();
  for (int level = coarsest_level; level < nlevel; ++level) {
    const std::size_t levelli = static_cast<std::size_t>(level);
    const int ngrow_Pv_on_cells = GetScratch(level).nGrow() + 1;
    {
      const ::amrex::BoxArray& ba = hier.GetPatchLevel(level).box_array;
      const ::amrex::DistributionMapping& dm =
          hier.GetPatchLevel(level).distribution_mapping;
      Pv_[levelli].on_cells = ::amrex::MultiFab(ba, dm, 2, ngrow_Pv_on_cells);
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      ::amrex::IntVect ngrow_Pv_on_faces =
          GetFluxes(level, Direction(dir)).nGrowVect();
      ngrow_Pv_on_faces[dir] += 1;
      const ::amrex::BoxArray& ba = GetFluxes(level, Direction(dir)).boxArray();
      const ::amrex::DistributionMapping& dm =
          GetFluxes(level, Direction(dir)).DistributionMap();
      Pv_[levelli].on_faces[static_cast<std::size_t>(dir)] =
          ::amrex::MultiFab(ba, dm, 1, ngrow_Pv_on_faces);
    }
  }
}

} // namespace fub::amrex
