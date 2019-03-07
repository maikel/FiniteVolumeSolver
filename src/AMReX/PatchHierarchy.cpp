#include "fub/grid/AMReX/PatchHierarchy.hpp"

namespace fub {
namespace amrex {
PatchLevel::PatchLevel(int level, Duration tp,
                       const ::amrex::BoxArray& box_array,
                       const ::amrex::DistributionMapping& distribution_mapping,
                       int n_components)
    : level_number{level}, time_point{tp}, data{box_array, distribution_mapping,
                                                n_components, 0} {}

PatchHierarchy::PatchHierarchy(DataDescription desc,
                               const CartesianGridGeometry& geometry,
                               const PatchHierarchyOptions& options)
    : description_{std::move(desc)}, grid_geometry_{geometry},
      options_{options}, patch_level_{}, patch_level_geometry_{} {
  patch_level_.resize(static_cast<std::size_t>(options.max_number_of_levels));
  patch_level_geometry_.resize(static_cast<std::size_t>(options.max_number_of_levels));
  ::amrex::IntVect lower{};
  ::amrex::IntVect upper{};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    upper[i] = geometry.cell_dimensions[static_cast<std::size_t>(i)] - 1;
  }
  ::amrex::Box level_box(lower, upper);
  for (::amrex::Geometry& geom : patch_level_geometry_) {
    geom = ::amrex::Geometry(level_box, &grid_geometry_.coordinates, -1,
                             grid_geometry_.periodicity.data());
    level_box.refine(2);
  }
}

} // namespace amrex
} // namespace fub
