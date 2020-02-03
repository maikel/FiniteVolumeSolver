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

#include "fub/AMReX/output/DebugOutput.hpp"

#include <algorithm>
#include <numeric>

namespace fub::amrex {
/// \brief Returns all the hierarchies which are stored via SaveData
const std::vector<DebugStorage::Hierarchy>& DebugStorage::GetHierarchies() const
    noexcept {
  return saved_hierarchies_;
}

/// \brief Returns all the component names which are stored via SaveData
const std::vector<DebugStorage::ComponentNames>& DebugStorage::GetNames() const
    noexcept {
  return names_per_hierarchy_;
}

void DebugStorage::SaveData(const ::amrex::MultiFab& mf,
                            const std::string& name,
                            ::amrex::SrcComp component) {
  Hierarchy single_level_hierarchy{};
  ::amrex::BoxArray ba = mf.boxArray();
  ::amrex::DistributionMapping dm = mf.DistributionMap();
  ::amrex::MultiFab& dest = single_level_hierarchy.emplace_back(ba, dm, 1, 0);
  ::amrex::MultiFab::Copy(dest, mf, component.i, 0, 1, 0);
  saved_hierarchies_.push_back(std::move(single_level_hierarchy));
  names_per_hierarchy_.push_back(std::vector<std::string>{name});
}

namespace {

struct equal_to {
  int i;
  const std::vector<std::vector<::amrex::BoxArray>>* bas_;
  bool operator()(int j) { return (*bas_)[i] == (*bas_)[j]; }
};

std::vector<std::vector<int>>
Partitions_(const std::vector<std::vector<::amrex::BoxArray>>& boxes) {
  std::vector<int> selection(boxes.size());
  std::vector<std::vector<int>> partitions{};
  auto first = selection.begin();
  auto last = selection.end();
  if (first == last) {
    return partitions;
  }
  std::iota(first, last, 0);
  auto next = std::stable_partition(first, last, equal_to{*first, &boxes});
  partitions.emplace_back(first, next);
  while (next != last) {
    first = next;
    next = std::stable_partition(first, last, equal_to{*first, &boxes});
    partitions.emplace_back(first, next);
  }
  return partitions;
}

template <typename Proj>
auto GatherBy(const std::vector<std::vector<::amrex::MultiFab>>& hierarchies,
              ::amrex::IndexType location, Proj projection) {
  using T = std::decay_t<decltype(projection(hierarchies[0][0]))>;
  std::vector<std::vector<T>> projected{};
  auto do_projection = [p = std::move(projection)](
                           const std::vector<::amrex::MultiFab>& hierarchy) {
    std::vector<T> projected(hierarchy.size());
    std::transform(hierarchy.begin(), hierarchy.end(), projected.begin(),
                   [proj = std::move(p)](const ::amrex::MultiFab& mf) {
                     return proj(mf);
                   });
    return projected;
  };
  for (const std::vector<::amrex::MultiFab>& hierarchy : hierarchies) {
    if (hierarchy[0].boxArray().ixType() == location) {
      projected.push_back(do_projection(hierarchy));
    }
  }
  return projected;
}

std::vector<std::vector<::amrex::BoxArray>>
GatherBoxArrays_(const std::vector<std::vector<::amrex::MultiFab>>& hierarchies,
                 ::amrex::IndexType location) {
  return GatherBy(hierarchies, location,
                  [](const ::amrex::MultiFab& mf) { return mf.boxArray(); });
}

std::vector<std::vector<::amrex::DistributionMapping>> GatherDistributionMaps_(
    const std::vector<std::vector<::amrex::MultiFab>>& hierarchies,
    ::amrex::IndexType location) {
  return GatherBy(hierarchies, location, [](const ::amrex::MultiFab& mf) {
    return mf.DistributionMap();
  });
}

std::vector<const std::vector<::amrex::MultiFab>*> GatherHierarchies_(
    const std::vector<std::vector<::amrex::MultiFab>>& hierarchies,
    ::amrex::IndexType location) {
  std::vector<const std::vector<::amrex::MultiFab>*> gathered{};
  for (const std::vector<::amrex::MultiFab>& hierarchy : hierarchies) {
    if (hierarchy[0].boxArray().ixType() == location) {
      gathered.push_back(&hierarchy);
    }
  }
  return gathered;
}

std::vector<std::vector<std::string>>
GatherNames_(const std::vector<std::vector<std::string>>& names,
             const std::vector<std::vector<::amrex::MultiFab>>& hierarchies,
             ::amrex::IndexType location) {
  std::vector<std::vector<std::string>> gathered{};
  auto components = names.begin();
  for (const std::vector<::amrex::MultiFab>& hierarchy : hierarchies) {
    if (hierarchy[0].boxArray().ixType() == location) {
      gathered.push_back(*components);
    }
    ++components;
  }
  return gathered;
}

std::vector<std::string>
Select_(const std::vector<std::vector<std::string>>& names,
        const std::vector<int>& partition) {
  using namespace std::literals;
  std::vector<std::string> selected;
  for (int i : partition) {
    selected.insert(selected.end(), names[i].begin(), names[i].end());
  }
  auto last = selected.end();
  for (std::size_t k = 0; k < selected.size(); ++k) {
    const std::string candidate = selected[k];
    auto first = selected.begin() + k + 1;
    auto pos = std::find(first, last, candidate);
    if (last != pos) {
      selected[k] += "_0"s;
      int counter = 1;
      do {
        *pos = fmt::format("{}_{}", *pos, counter);
        counter += 1;
        pos = std::find(first, last, candidate);
      } while (last != pos);
    }
  }
  return selected;
}
} // namespace

std::vector<std::pair<std::vector<::amrex::MultiFab>, std::vector<std::string>>>
DebugStorage::GatherFields(::amrex::IndexType location) const {
  using namespace ::amrex;
  std::vector<std::vector<BoxArray>> box_arrays =
      GatherBoxArrays_(saved_hierarchies_, location);
  std::vector<std::vector<DistributionMapping>> distribution_maps =
      GatherDistributionMaps_(saved_hierarchies_, location);
  std::vector<const std::vector<::amrex::MultiFab>*> filtered_hierarchies =
      GatherHierarchies_(saved_hierarchies_, location);
  std::vector<std::vector<int>> partitions = Partitions_(box_arrays);
  std::vector<
      std::pair<std::vector<::amrex::MultiFab>, std::vector<std::string>>>
      hierarchies{};
  for (const std::vector<int>& partition : partitions) {
    std::vector<std::vector<std::string>> all_components =
        GatherNames_(names_per_hierarchy_, saved_hierarchies_, location);
    std::vector<std::string> components = Select_(all_components, partition);
    int n_components = static_cast<int>(components.size());
    std::vector<::amrex::MultiFab> hierarchy{};
    std::vector<BoxArray> box_array = box_arrays[partition[0]];
    std::vector<DistributionMapping> distribution =
        distribution_maps[partition[0]];
    for (std::size_t level = 0; level < box_array.size(); ++level) {
      MultiFab& dest = hierarchy.emplace_back(
          box_array[level], distribution[level], n_components, 0);
      int comp = 0;
      for (int i : partition) {
        const MultiFab& src = (*filtered_hierarchies[i])[level];
        const int n = src.nComp();
        MultiFab::Copy(dest, src, 0, comp, n, 0);
        comp += n;
      }
      FUB_ASSERT(comp == n_components);
    }
    hierarchies.emplace_back(std::move(hierarchy), std::move(components));
  }
  return hierarchies;
}

void DebugStorage::ClearAll() {
  names_per_hierarchy_.clear();
  saved_hierarchies_.clear();
}

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
        ::amrex::Vector<std::string> varnames(all_names[partition_counter]);
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
          on_cell.setVal(0.0);
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
