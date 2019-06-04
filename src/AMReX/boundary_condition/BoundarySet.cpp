#include "fub/AMReX/boundary_condition/BoundarySet.hpp"

namespace fub::amrex {

void BoundarySet::FillBoundary(::amrex::MultiFab& mf,
                               const ::amrex::Geometry& geom,
                               Duration timepoint,
                               const GriddingAlgorithm& gridding) {
  for (BoundaryCondition& condition : conditions) {
    condition.FillBoundary(mf, geom, timepoint, gridding);
  }
}

} // namespace fub::amrex
