#include "fub/grid/AMReX/PatchHierarchy.hpp"

namespace fub {
namespace amrex {
PatchLevel::PatchLevel(int level, Duration tp,
                       const ::amrex::BoxArray& box_array,
                       const ::amrex::DistributionMapping& distribution_mapping,
                       int n_components)
    : level_number{level}, time_point{tp}, data{box_array, distribution_mapping,
                                                n_components, 0} {}

PatchLevel::PatchLevel(int level, Duration tp,
                       const ::amrex::BoxArray& box_array,
                       const ::amrex::DistributionMapping& distribution_mapping,
                       int n_components,
                       const ::amrex::FabFactory<::amrex::FArrayBox>& factory)
    : level_number{level},
      time_point{tp}, data{box_array, distribution_mapping, n_components,
                           0,         ::amrex::MFInfo(),    factory} {}

PatchHierarchy::PatchHierarchy(DataDescription desc,
                               const CartesianGridGeometry& geometry,
                               const PatchHierarchyOptions& options)
    : description_{std::move(desc)}, grid_geometry_{geometry},
      options_{options}, patch_level_{}, patch_level_geometry_{} {
  patch_level_.resize(static_cast<std::size_t>(options.max_number_of_levels));
  patch_level_geometry_.resize(
      static_cast<std::size_t>(options.max_number_of_levels));
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

void WriteCheckpointFile(const std::string checkpointname,
                         const fub::amrex::PatchHierarchy& hier) {
  const int nlevels = hier.GetNumberOfLevels();
  ::amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);
  // write Header file
  if (::amrex::ParallelDescriptor::IOProcessor()) {
    const std::string header(checkpointname + "/Header");
    ::amrex::VisMF::IO_Buffer io_buffer(::amrex::VisMF::IO_Buffer_Size);
    std::ofstream hout;
    hout.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    hout.open(header.c_str(), std::ofstream::out | std::ofstream::trunc |
                                  std::ofstream::binary);
    if (!hout.good()) {
      ::amrex::FileOpenFailed(header);
    }

    hout.precision(17);

    // write out title line
    hout << "Checkpoint file for AmrCoreAdv\n";

    // write out finest_level
    hout << nlevels - 1 << '\n';

    // write out array of istep
    for (int level = 0; level < nlevels; ++level) {
      hout << hier.GetCycles(level) << ' ';
    }
    hout << '\n';

    // write out array of t_new
    for (int level = 0; level < nlevels; ++level) {
      hout << hier.GetTimePoint(level).count() << ' ';
    }
    hout << '\n';

    // write the BoxArray at each level
    for (int level = 0; level < nlevels; ++level) {
      hier.GetPatchLevel(level).data.boxArray().writeOn(hout);
      hout << '\n';
    }
  }

  // write the MultiFab data to, e.g., chk00010/Level_0/
  for (int level = 0; level < nlevels; ++level) {
    ::amrex::VisMF::Write(hier.GetPatchLevel(level).data,
                          ::amrex::MultiFabFileFullPrefix(level, checkpointname,
                                                          "Level_", "data"));
  }
}

} // namespace amrex
} // namespace fub
