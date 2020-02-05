// Copyright (c) 2020 Maikel Nadolski
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

#include "fub/AMReX/cutcell/output/DebugOutput.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

#include <algorithm>
#include <numeric>

namespace fub::amrex::cutcell {

DebugOutput::DebugOutput(const ProgramOptions& opts)
    : OutputAtFrequencyOrInterval(opts) {
  directory_ = GetOptionOr(opts, "directory", directory_);
}

void DebugOutput::operator()(const GriddingAlgorithm& grid) {
  DebugStorage& storage = *grid.GetPatchHierarchy().GetDebugStorage();

  std::vector<
      std::pair<std::vector<::amrex::MultiFab>, std::vector<std::string>>>
      hierarchies = storage.GatherFields(::amrex::IndexType::TheCellType());

  int partition_counter = 0;
  for (auto&& [hierarchy, names] : hierarchies) {
    const std::string plotfilename =
        fmt::format("{}/partition_{}_plt{:09}", directory_, partition_counter,
                    grid.GetPatchHierarchy().GetCycles());
    const std::size_t size = hierarchy.size();
    const double time_point = grid.GetTimePoint().count();
    ::amrex::Vector<const ::amrex::MultiFab*> mf(size);
    ::amrex::Vector<::amrex::Geometry> geoms(size);
    ::amrex::Vector<int> level_steps(size);
    ::amrex::Vector<::amrex::IntVect> ref_ratio(size);
    for (std::size_t i = 0; i < size; ++i) {
      mf[i] = &hierarchy[i];
      const int ii = static_cast<int>(i);
      geoms[i] = grid.GetPatchHierarchy().GetGeometry(ii);
      level_steps[i] = static_cast<int>(grid.GetPatchHierarchy().GetCycles(ii));
      ref_ratio[i] = grid.GetPatchHierarchy().GetRatioToCoarserLevel(ii);
    }
    ::amrex::Vector<std::string> vnames(names.begin(), names.end());
    ::amrex::WriteMultiLevelPlotfile(
        plotfilename, size, mf, vnames, geoms, time_point, level_steps,
        ref_ratio, "HyperCLaw-V1.1", "Level_", "Cell", {"raw_fields"});
    partition_counter += 1;
  }

  const std::vector<std::vector<::amrex::MultiFab>>& all_hierarchies =
      storage.GetHierarchies();
  const std::vector<std::vector<std::string>>& all_names = storage.GetNames();

  auto hier_fields = all_names.begin();
  const std::string level_prefix = "Level_";
  partition_counter = 0;
  ::amrex::VisMF::SetHeaderVersion(::amrex::VisMF::Header::Version_v1);

  for (const std::vector<::amrex::MultiFab>& hierarchy : all_hierarchies) {
    if (hierarchy[0].ixType() != ::amrex::IndexType::TheCellType()) {
      const std::string plotfilename = fmt::format(
          "{}/{}_{}_plt{:09}", directory_, all_names[partition_counter][0],
          partition_counter, grid.GetPatchHierarchy().GetCycles());
      const std::string raw_data = plotfilename + "/raw_fields";
      const int nlevels = hierarchy.size();
      ::amrex::PreBuildDirectorHierarchy(plotfilename, "Level_", nlevels, true);
      ::amrex::PreBuildDirectorHierarchy(raw_data, "Level_", nlevels, true);
      // write Header file
      if (::amrex::ParallelDescriptor::IOProcessor()) {
        const std::string header(plotfilename + "/Header");
        ::amrex::Vector<::amrex::BoxArray> box_arrays(hierarchy.size());
        ::amrex::Vector<::amrex::Geometry> geometries(hierarchy.size());
        ::amrex::Vector<int> level_steps(hierarchy.size());
        ::amrex::Vector<::amrex::IntVect> ref_ratio(hierarchy.size());
        ::amrex::Vector<std::string> varnames(
            all_names[partition_counter].begin(),
            all_names[partition_counter].end());
        double time_point = grid.GetPatchHierarchy().GetTimePoint(0).count();

        for (int level = 0; level < nlevels; ++level) {
          box_arrays[level] = hierarchy[level].boxArray();
          geometries[level] = grid.GetPatchHierarchy().GetGeometry(level);
          level_steps[level] = grid.GetPatchHierarchy().GetCycles(level);
          ref_ratio[level] =
              grid.GetPatchHierarchy().GetRatioToCoarserLevel(level);
        }

        ::amrex::VisMF::IO_Buffer io_buffer(::amrex::VisMF::IO_Buffer_Size);
        std::ofstream header_out;
        header_out.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        header_out.open(header.c_str(), std::ofstream::out |
                                            std::ofstream::trunc |
                                            std::ofstream::binary);
        ::amrex::WriteGenericPlotfileHeader(header_out, nlevels, box_arrays,
                                            varnames, geometries, time_point,
                                            level_steps, ref_ratio);
      }
      // Force processors to wait until directory has been built.
      ::amrex::ParallelDescriptor::Barrier();
      for (int level = 0; level < nlevels; ++level) {
        std::string prefix =
            ::amrex::MultiFabFileFullPrefix(level, plotfilename);
        std::string prefix_raw = ::amrex::MultiFabFileFullPrefix(
            level, raw_data, level_prefix, hier_fields->front() + "_raw");
        // if (::amrex::ParallelDescriptor::IOProcessor()) {
        //   if (!::amrex::UtilCreateDirectory(prefix, 0755)) {
        //     ::amrex::CreateDirectoryFailed(prefix);
        //   }
        // }
        //
        // Force other processors to wait until directory is built.
        //
        // ::amrex::ParallelDescriptor::Barrier();
        ::amrex::BoxArray ba = hierarchy[level].boxArray();
        ba.enclosedCells();
        ::amrex::DistributionMapping dm = hierarchy[level].DistributionMap();
        ::amrex::MultiFab on_cell(ba, dm, hierarchy[level].nComp(), 0);
        if (hierarchy[level].ixType() == ::amrex::IndexType::TheNodeType()) {
          ::amrex::average_node_to_cellcenter(on_cell, 0, hierarchy[level], 0,
                                              hierarchy[level].nComp(), 0);
        } else {
          ForEachFab(fub::execution::openmp, on_cell, [&](const ::amrex::MFIter& mfi) {
            ::amrex::FArrayBox& cells = on_cell[mfi];
            const ::amrex::FArrayBox& faces = hierarchy[level][mfi];
            for (int icomp = 0; icomp < hierarchy[level].nComp(); ++icomp) {
            ForEachIndex(mfi.tilebox(), [&](auto... is) {
              ::amrex::IntVect iv{int(is)...};
              cells(iv, icomp) = faces(iv, icomp);
            });
            }
          });
        }
        ::amrex::VisMF::Write(on_cell, prefix);
        ::amrex::VisMF::Write(hierarchy[level], prefix_raw);
      }
    }
    ++partition_counter;
    ++hier_fields;
  }

  storage.ClearAll();
}

} // namespace fub::amrex
