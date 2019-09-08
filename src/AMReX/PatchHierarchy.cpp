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

#include "fub/AMReX/PatchHierarchy.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

#ifdef AMREX_USE_EB
#include "fub/AMReX/cutcell/IndexSpace.hpp"
#include <AMReX_EB2.H>
#endif

namespace fub {
namespace amrex {

PatchLevel::PatchLevel(const PatchLevel& other)
    : level_number(other.level_number), time_point(other.time_point),
      cycles(other.cycles), box_array(other.box_array.boxList()),
      distribution_mapping(other.distribution_mapping.ProcessorMap()),
      data(box_array, distribution_mapping, other.data.nComp(),
           other.data.nGrowVect(), ::amrex::MFInfo(), other.data.Factory()) {
  data.copy(other.data);
}

PatchLevel& PatchLevel::operator=(const PatchLevel& other) {
  PatchLevel tmp(other);
  return *this = std::move(tmp);
}

PatchLevel::PatchLevel(int level, Duration tp, const ::amrex::BoxArray& ba,
                       const ::amrex::DistributionMapping& dm, int n_components)
    : level_number{level}, time_point{tp}, box_array{ba},
      distribution_mapping{dm}, data{box_array, distribution_mapping,
                                     n_components, 0} {}

PatchLevel::PatchLevel(int level, Duration tp, const ::amrex::BoxArray& ba,
                       const ::amrex::DistributionMapping& dm, int n_components,
                       const ::amrex::FabFactory<::amrex::FArrayBox>& factory)
    : level_number{level}, time_point{tp}, box_array{ba},
      distribution_mapping{dm}, data{box_array,         distribution_mapping,
                                     n_components,      0,
                                     ::amrex::MFInfo(), factory} {}

PatchHierarchy::PatchHierarchy(DataDescription desc,
                               const CartesianGridGeometry& geometry,
                               const PatchHierarchyOptions& options)
    : description_{std::move(desc)}, grid_geometry_{geometry},
      options_{options}, patch_level_{}, patch_level_geometry_{} {
  patch_level_.reserve(static_cast<std::size_t>(options.max_number_of_levels));
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
    level_box.refine(options_.refine_ratio);
  }
#ifdef AMREX_USE_EB
  auto shop = ::amrex::EB2::makeShop(::amrex::EB2::AllRegularIF());
  index_spaces_ = cutcell::MakeIndexSpaces(shop, patch_level_geometry_[0], options.max_number_of_levels, options.refine_ratio);
#endif
}

int PatchHierarchy::GetRatioToCoarserLevel(int level, Direction dir) const
    noexcept {
  if (level == 0) {
    return 1;
  }
  return options_.refine_ratio[static_cast<int>(dir)];
}

::amrex::IntVect PatchHierarchy::GetRatioToCoarserLevel(int level) const
    noexcept {
  if (level == 0) {
    return ::amrex::IntVect::TheUnitVector();
  }
  return options_.refine_ratio;
}

const ::amrex::Geometry& PatchHierarchy::GetGeometry(int level) const {
  return patch_level_geometry_[static_cast<std::size_t>(level)];
}

const PatchHierarchyOptions& PatchHierarchy::GetOptions() const noexcept {
  return options_;
}

const CartesianGridGeometry& PatchHierarchy::GetGridGeometry() const noexcept {
  return grid_geometry_;
}

std::ptrdiff_t PatchHierarchy::GetCycles(int level) const {
  return GetPatchLevel(level).cycles;
}

Duration PatchHierarchy::GetTimePoint(int level) const {
  return GetPatchLevel(level).time_point;
}

int PatchHierarchy::GetNumberOfLevels() const noexcept {
  return static_cast<int>(patch_level_.size());
}

int PatchHierarchy::GetMaxNumberOfLevels() const noexcept {
  return GetOptions().max_number_of_levels;
}

PatchLevel& PatchHierarchy::GetPatchLevel(int level) {
  return patch_level_.at(static_cast<std::size_t>(level));
}

span<const ::amrex::EB2::IndexSpace*> PatchHierarchy::GetIndexSpaces() noexcept {
  return index_spaces_;
}

const PatchLevel& PatchHierarchy::GetPatchLevel(int level) const {
  return patch_level_[static_cast<std::size_t>(level)];
}

void PatchHierarchy::PushBack(const PatchLevel& level) {
  patch_level_.push_back(level);
}

void PatchHierarchy::PushBack(PatchLevel&& level) {
  patch_level_.push_back(std::move(level));
}

void PatchHierarchy::PopBack() { patch_level_.pop_back(); }

const DataDescription& PatchHierarchy::GetDataDescription() const noexcept {
  return description_;
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

PatchHierarchy ReadCheckpointFile(const std::string& checkpointname,
                                  DataDescription desc,
                                  const CartesianGridGeometry& geometry,
                                  const PatchHierarchyOptions& options) {
  std::string File(checkpointname + "/Header");
  ::amrex::VisMF::IO_Buffer io_buffer(
      static_cast<std::size_t>(::amrex::VisMF::GetIOBufferSize()));
  ::amrex::Vector<char> fileCharPtr;
  ::amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
  std::string fileCharPtrString(fileCharPtr.dataPtr());
  std::istringstream is(fileCharPtrString, std::istringstream::in);

  PatchHierarchy hierarchy(desc, geometry, options);

  // read in title line
  std::string line;
  std::getline(is, line);

  // read in finest_level
  std::getline(is, line);
  const int finest_level = std::stoi(line);
  std::vector<std::ptrdiff_t> cycles(
      static_cast<std::size_t>(finest_level + 1));
  std::vector<double> time_points(static_cast<std::size_t>(finest_level + 1));

  // read in array of istep
  std::getline(is, line);
  {
    std::istringstream lis(line);
    for (std::size_t i = 0; i <= static_cast<std::size_t>(finest_level); ++i) {
      lis >> cycles[i];
    }
  }

  // read in array of t_new
  std::getline(is, line);
  {
    std::istringstream lis(line);
    for (std::size_t i = 0; i <= static_cast<std::size_t>(finest_level); ++i) {
      lis >> time_points[i];
    }
  }

  for (int lev = 0; lev <= finest_level; ++lev) {
    // read in level 'lev' BoxArray from Header
    ::amrex::BoxArray ba;
    ba.readFrom(is);

    // create a distribution mapping
    ::amrex::DistributionMapping dm{ba, ::amrex::ParallelDescriptor::NProcs()};

    hierarchy.PushBack(
        PatchLevel(lev, Duration(time_points[static_cast<std::size_t>(lev)]),
                   ba, dm, desc.n_state_components));
    hierarchy.GetPatchLevel(lev).cycles = cycles[static_cast<std::size_t>(lev)];
  }

  // read in the MultiFab data
  for (int lev = 0; lev <= finest_level; ++lev) {
    ::amrex::VisMF::Read(
        hierarchy.GetPatchLevel(lev).data,
        ::amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "data"));
  }

  return hierarchy;
}


  void WriteMatlabData(std::ostream& out, const ::amrex::FArrayBox& fab,
                     const fub::IdealGasMix<1>& eq,
                       const ::amrex::Geometry& geom) {
  auto view = MakeView<const Complete<IdealGasMix<1>>>(fab, eq);
  out << fmt::format("{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}",
      "X", "Density", "VelocityX", "SpeedOfSound", "Temperature",
      "Pressure", "Gamma", "HeatCapacityAtConstantPressure");
  const int nspecies = eq.GetReactor().GetNSpecies();
  span<const std::string> names = eq.GetReactor().GetSpeciesNames();
  for (int s = 0; s < nspecies - 1; ++s) {
    out << names[s] << ' ';
  }
  out << names[nspecies - 1] << '\n';
  ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i) {
    double x[3] = {0.0, 0.0, 0.0};
    ::amrex::IntVect iv{AMREX_D_DECL(int(i), 0, 0)};
    geom.CellCenter(iv, x);
    const double density = view.density(i);
    const double velocity_x =
        density > 0.0 ? view.momentum(i, 0) / density : 0.0;
    const double temperature = view.temperature(i);
    const double speed_of_sound = view.speed_of_sound(i);
    const double gamma = view.gamma(i);
    const double c_p = view.c_p(i);
    const double pressure = view.pressure(i);
    out << fmt::format("{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< "
                       "24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}",
                       x[0], density, velocity_x, speed_of_sound, temperature,
                       pressure, gamma, c_p);
    for (int s = 0; s < nspecies; ++s) {
      out << fmt::format("{:< 24.15e}", view.species(i, s));
    }
    out << '\n';
  });
}

#if AMREX_SPACEDIM == 3
void WriteMatlabData(std::ostream& out, const ::amrex::FArrayBox& fab,
                     const fub::IdealGasMix<3>& eq,
                       const ::amrex::Geometry& geom) {
  auto view = MakeView<const Complete<IdealGasMix<3>>>(fab, eq);
  out << fmt::format("{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}",
      "X", "Y", "Density", "VelocityX", "VelocityY", "VelocityZ", "SpeedOfSound", "Temperature",
      "Pressure", "Gamma", "HeatCapacityAtConstantPressure");
  const int nspecies = eq.GetReactor().GetNSpecies();
  span<const std::string> names = eq.GetReactor().GetSpeciesNames();
  for (int s = 0; s < nspecies - 1; ++s) {
    out << fmt::format("{:24s}", names[s]);
  }
  out << fmt::format("{:24s}\n", names[nspecies - 1]);
  ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i, std::ptrdiff_t j,
                                 std::ptrdiff_t k) {
    double x[3] = {0.0, 0.0, 0.0};
    ::amrex::IntVect iv{int(i), int(j), int(k)};
    geom.CellCenter(iv, x);
    const double density = view.density(i, j, k);
    const double velocity_x =
        density > 0.0 ? view.momentum(i, j, k, 0) / density : 0.0;
    const double velocity_y =
        density > 0.0 ? view.momentum(i, j, k, 1) / density : 0.0;
    const double velocity_z =
        density > 0.0 ? view.momentum(i, j, k, 2) / density : 0.0;
    const double temperature = view.temperature(i, j, k);
    const double speed_of_sound = view.speed_of_sound(i, j, k);
    const double gamma = view.gamma(i, j, k);
    const double c_p = view.c_p(i, j, k);
    const double pressure = view.pressure(i, j, k);
    out << fmt::format("{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}"
                       "{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}"
                       "{:< 24.15e}{:< 24.15e}{:< 24.15e}",
                       x[0], x[1], density, velocity_x, velocity_y, velocity_z,
                       speed_of_sound, temperature, pressure, gamma, c_p);
    for (int s = 0; s < nspecies; ++s) {
      out << fmt::format("{:< 24.15e}", view.species(i, j, k, s));
    }
    out << '\n';
  });
}
#endif

std::vector<double>
GatherStates(const PatchHierarchy& hierarchy,
             basic_mdspan<const double, extents<AMREX_SPACEDIM, dynamic_extent>> xs,
             MPI_Comm comm) {
  const int nlevel = hierarchy.GetNumberOfLevels();
  const int finest_level = nlevel - 1;
  const int ncomp = hierarchy.GetDataDescription().n_state_components;
  std::vector<double> buffer(xs.extent(1) * ncomp * nlevel);
  mdspan<double, 3> states(buffer.data(), xs.extent(1), ncomp, nlevel);
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(level).data;
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(level);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      ForEachIndex(mfi.tilebox(), [&](auto... is) {
        double lo[AMREX_SPACEDIM]{};
        double hi[AMREX_SPACEDIM]{};
        const ::amrex::IntVect iv{int(is)...};
        level_geom.LoNode(iv, lo);
        level_geom.HiNode(iv, hi);
        for (int k = 0; k < xs.extent(1); ++k) {
          bool is_in_range = true;
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            is_in_range = is_in_range && lo[d] <= xs(d, k) && xs(d, k) < hi[d];
          }
          if (is_in_range) {
            for (int comp = 0; comp < level_data.nComp(); ++comp) {
              states(k, comp, level) = level_data[mfi](iv, comp);
            }
          }
        }
      });
    });
  }
  std::vector<double> global_buffer(buffer.size());
  ::MPI_Allreduce(buffer.data(), global_buffer.data(), global_buffer.size(),
                  MPI_DOUBLE, MPI_SUM, comm);
  states = mdspan<double, 3>(global_buffer.data(), xs.extent(1), ncomp, nlevel);
  for (int level = 1; level < nlevel; ++level) {
    for (int comp = 0; comp < ncomp; ++comp) {
      for (int i = 0; i < xs.extent(1); ++i) {
        if (states(i, comp, level) == 0.0) {
          states(i, comp, level) = states(i, comp, level - 1);
        }
      }
    }
  }
  std::vector<double> result(&states(0, 0, finest_level),
                             &states(0, 0, finest_level) +
                                 xs.extent(1) * ncomp);
  return result;
}

  
void WriteTubeData(std::ostream& out, const PatchHierarchy& hierarchy,
                   const IdealGasMix<1>& eq, fub::Duration time_point,
                   std::ptrdiff_t cycle_number, MPI_Comm comm) {
  const std::size_t n_level =
      static_cast<std::size_t>(hierarchy.GetNumberOfLevels());
  std::vector<::amrex::FArrayBox> fabs{};
  fabs.reserve(n_level);
  for (std::size_t level = 0; level < n_level; ++level) {
    const int ilvl = static_cast<int>(level);
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(ilvl);
    ::amrex::Box domain = level_geom.Domain();
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(ilvl).data;
    ::amrex::FArrayBox local_fab(domain, level_data.nComp());
    local_fab.setVal(0.0);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& patch_data = level_data[mfi];
      local_fab.copy(patch_data);
    });
    int rank = -1;
    ::MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      ::amrex::FArrayBox& fab = fabs.emplace_back(domain, level_data.nComp());
      fab.setVal(0.0);
      ::MPI_Reduce(local_fab.dataPtr(), fab.dataPtr(), local_fab.size(),
                   MPI_DOUBLE, MPI_SUM, 0, comm);
      if (level > 0) {
        for (int comp = 1; comp < level_data.nComp(); ++comp) {
          for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
            ::amrex::IntVect fine_i{AMREX_D_DECL(i, domain.smallEnd(1), domain.smallEnd(2))};
            ::amrex::IntVect coarse_i = fine_i;
            coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
            if (fab(fine_i, 0) == 0.0) {
              fab(fine_i, comp) = fabs[level - 1](coarse_i, comp);
            }
          }
        }
        for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
          ::amrex::IntVect fine_i{AMREX_D_DECL(i, domain.smallEnd(1), domain.smallEnd(2))};
          ::amrex::IntVect coarse_i = fine_i;
          coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
          if (fab(fine_i, 0) == 0.0) {
            fab(fine_i, 0) = fabs[level - 1](coarse_i, 0);
          }
        }
      }
      if (level == n_level - 1) {
        out << fmt::format("nx = {}\n", domain.length(0));
        out << fmt::format("t = {}\n", time_point.count());
        out << fmt::format("cycle = {}\n", cycle_number);
        WriteMatlabData(out, fab, eq, level_geom);
        out.flush();
      }
    } else {
      ::MPI_Reduce(local_fab.dataPtr(), nullptr, local_fab.size(), MPI_DOUBLE,
                   MPI_SUM, 0, comm);
    }
  }
}

} // namespace amrex
} // namespace fub
