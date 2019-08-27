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

} // namespace amrex
} // namespace fub
