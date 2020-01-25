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

#include "fub/AMReX/cutcell/output/WriteHdf5.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

#include <array>

namespace fub::amrex::cutcell {

WriteHdf5::WriteHdf5(const std::map<std::string, pybind11::object>& vm)
    : OutputAtFrequencyOrInterval<GriddingAlgorithm>(vm) {
  path_to_file_ = GetOptionOr(vm, "path", std::string("grid.h5"));
  auto it = vm.find("box");
  if (it != vm.end()) {
    auto box = ToMap(it->second);
    std::array<int, AMREX_SPACEDIM> lo =
        GetOptionOr(box, "lower", std::array<int, AMREX_SPACEDIM>{});
    std::array<int, AMREX_SPACEDIM> hi =
        GetOptionOr(box, "upper", std::array<int, AMREX_SPACEDIM>{});
    output_box_.emplace(::amrex::IntVect(lo), ::amrex::IntVect(hi));
  }
  if (boost::filesystem::is_regular_file(path_to_file_)) {
    boost::log::sources::severity_logger<boost::log::trivial::severity_level>
        log(boost::log::keywords::severity = boost::log::trivial::info);
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "HDF5");
    for (int i = 1; i < std::numeric_limits<int>::max(); ++i) {
      std::string new_name = fmt::format("{}.{}", path_to_file_, i);
      if (!boost::filesystem::exists(new_name)) {
        BOOST_LOG(log) << fmt::format(
            "Old output file '{}' detected. Rename old file to '{}'",
            path_to_file_, new_name);
        boost::filesystem::rename(path_to_file_, new_name);
        break;
      }
    }
  } else if (boost::filesystem::exists(path_to_file_)) {
    boost::log::sources::severity_logger<boost::log::trivial::severity_level>
        log(boost::log::keywords::severity = boost::log::trivial::warning);
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "HDF5");
    for (int i = 1; i < std::numeric_limits<int>::max(); ++i) {
      std::string new_name = fmt::format("{}.{}", path_to_file_, i);
      if (!boost::filesystem::exists(new_name)) {
        BOOST_LOG(log) << fmt::format(
            "Path'{}' points to some non-file. Output will be directory to "
            "'{}' instead.",
            path_to_file_, new_name);
        path_to_file_ = new_name;
        break;
      }
    }
  }
}

namespace {
void WriteHdf5UnRestricted(const std::string& name,
                           const PatchHierarchy& hierarchy) {
  MPI_Comm comm = ::amrex::ParallelContext::CommunicatorAll();
  int rank = -1;
  ::MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    boost::filesystem::path path(name);
    boost::filesystem::path dir = boost::filesystem::absolute(
        path.parent_path(), boost::filesystem::current_path());
    if (!boost::filesystem::exists(dir)) {
      boost::filesystem::create_directories(dir);
    }
  }
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
    if (rank == 0) {
      ::amrex::FArrayBox& fab = fabs.emplace_back(domain, level_data.nComp());
      fab.setVal(0.0);
      ::MPI_Reduce(local_fab.dataPtr(), fab.dataPtr(), local_fab.size(),
                   MPI_DOUBLE, MPI_SUM, 0, comm);
      if (level > 0) {
        for (int comp = 1; comp < level_data.nComp(); ++comp) {
          ForEachIndex(domain, [&](auto... is) {
            ::amrex::IntVect fine_i{int(is)...};
            ::amrex::IntVect coarse_i = fine_i;
            coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
            if (fab(fine_i, 0) == 0.0) {
              fab(fine_i, comp) = fabs[level - 1](coarse_i, comp);
            }
          });
        }
        ForEachIndex(domain, [&](auto... is) {
          ::amrex::IntVect fine_i{int(is)...};
          ::amrex::IntVect coarse_i = fine_i;
          coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
          if (fab(fine_i, 0) == 0.0) {
            fab(fine_i, 0) = fabs[level - 1](coarse_i, 0);
          }
        });
      }
      if (level == n_level - 1) {
        const fub::Duration time_point = hierarchy.GetTimePoint();
        const std::ptrdiff_t cycle_number = hierarchy.GetCycles();
        WriteToHDF5(name, fab, level_geom, time_point, cycle_number);
      }
    } else {
      ::MPI_Reduce(local_fab.dataPtr(), nullptr, local_fab.size(), MPI_DOUBLE,
                   MPI_SUM, 0, comm);
    }
  }
}

void WriteHdf5RestrictedToBox(const std::string& name,
                              const PatchHierarchy& hierarchy,
                              const ::amrex::Box& finest_box) {
  MPI_Comm comm = ::amrex::ParallelContext::CommunicatorAll();
  int rank = -1;
  ::MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    boost::filesystem::path path(name);
    boost::filesystem::path dir = boost::filesystem::absolute(
        path.parent_path(), boost::filesystem::current_path());
    if (!boost::filesystem::exists(dir)) {
      boost::filesystem::create_directories(dir);
    }
  }
  const std::size_t n_level =
      static_cast<std::size_t>(hierarchy.GetNumberOfLevels());
  std::vector<::amrex::Geometry> geoms{};
  std::vector<::amrex::MultiFab> data{};
  std::vector<::amrex::FArrayBox> fabs{};
  std::vector<::amrex::Box> boxes{};
  boxes.reserve(n_level);
  boxes.push_back(finest_box);
  ::amrex::Box box = finest_box;
  for (std::size_t level = n_level - 1; level > 0; --level) {
    ::amrex::IntVect refine_ratio = hierarchy.GetRatioToCoarserLevel(level);
    box.coarsen(refine_ratio);
    boxes.push_back(box);
  }
  std::reverse(boxes.begin(), boxes.end());
  geoms.reserve(n_level);
  data.reserve(n_level);
  fabs.reserve(n_level);
  for (std::size_t level = 0; level < n_level; ++level) {
    const int ilvl = static_cast<int>(level);
    ::amrex::Box domain = boxes[level];
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(ilvl).data;
    ::amrex::FArrayBox local_fab(domain, level_data.nComp());
    local_fab.setVal(0.0);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& patch_data = level_data[mfi];
      const ::amrex::Box box = mfi.tilebox() & domain;
      local_fab.copy(patch_data, box);
    });
    if (rank == 0) {
      ::amrex::FArrayBox& fab = fabs.emplace_back(domain, level_data.nComp());
      fab.setVal(0.0);
      ::MPI_Reduce(local_fab.dataPtr(), fab.dataPtr(), local_fab.size(),
                   MPI_DOUBLE, MPI_SUM, 0, comm);
      if (level > 0) {
        for (int comp = 1; comp < level_data.nComp(); ++comp) {
          ForEachIndex(domain, [&](auto... is) {
            ::amrex::IntVect fine_i{int(is)...};
            ::amrex::IntVect coarse_i = fine_i;
            coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
            if (fab(fine_i, 0) == 0.0) {
              fab(fine_i, comp) = fabs[level - 1](coarse_i, comp);
            }
          });
        }
        ForEachIndex(domain, [&](auto... is) {
          ::amrex::IntVect fine_i{int(is)...};
          ::amrex::IntVect coarse_i = fine_i;
          coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
          if (fab(fine_i, 0) == 0.0) {
            fab(fine_i, 0) = fabs[level - 1](coarse_i, 0);
          }
        });
      }
      if (level == n_level - 1) {
        const fub::Duration time_point = hierarchy.GetTimePoint();
        const std::ptrdiff_t cycle_number = hierarchy.GetCycles();
        const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(ilvl);
        WriteToHDF5(name, fab, level_geom, time_point, cycle_number);
      }
    } else {
      ::MPI_Reduce(local_fab.dataPtr(), nullptr, local_fab.size(), MPI_DOUBLE,
                   MPI_SUM, 0, comm);
    }
  }
}
} // namespace

void WriteHdf5::operator()(const GriddingAlgorithm& grid) {
  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "HDF5");
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", grid.GetTimePoint().count());
  BOOST_LOG(log) << fmt::format("Write Hdf5 output to '{}'.", path_to_file_);
  if (output_box_) {
    WriteHdf5RestrictedToBox(path_to_file_, grid.GetPatchHierarchy(),
                             *output_box_);
  } else {
    WriteHdf5UnRestricted(path_to_file_, grid.GetPatchHierarchy());
  }
}

} // namespace fub::amrex::cutcell
