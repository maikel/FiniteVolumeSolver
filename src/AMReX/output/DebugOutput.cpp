// Copyright (c) 2020 Maikel Nadolski
// Copyright (c) 2020 Stefan Vater
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
#include "fub/AMReX/ForEachFab.hpp"

#include <algorithm>
#include <numeric>

namespace fub::amrex {

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
auto GatherBy(const std::vector<DebugSnapshot::Hierarchy>& hierarchies,
              ::amrex::IndexType location, Proj projection) {
  using T = std::decay_t<decltype(projection(hierarchies[0][0]))>;
  std::vector<std::vector<T>> projected{};
  auto do_projection =
      [p = std::move(projection)](const DebugSnapshot::Hierarchy& hierarchy) {
        std::vector<T> projected(hierarchy.size());
        std::transform(hierarchy.begin(), hierarchy.end(), projected.begin(),
                       [proj = std::move(p)](const ::amrex::MultiFab& mf) {
                         return proj(mf);
                       });
        return projected;
      };
  for (const DebugSnapshot::Hierarchy& hierarchy : hierarchies) {
    if (hierarchy[0].boxArray().ixType() == location) {
      projected.push_back(do_projection(hierarchy));
    }
  }
  return projected;
}

std::vector<std::vector<::amrex::BoxArray>>
GatherBoxArrays_(const std::vector<DebugSnapshot::Hierarchy>& hierarchies,
                 ::amrex::IndexType location) {
  return GatherBy(hierarchies, location,
                  [](const ::amrex::MultiFab& mf) { return mf.boxArray(); });
}

std::vector<std::vector<::amrex::DistributionMapping>>
GatherDistributionMaps_(const std::vector<DebugSnapshot::Hierarchy>& hierarchies,
                        ::amrex::IndexType location) {
  return GatherBy(hierarchies, location, [](const ::amrex::MultiFab& mf) {
    return mf.DistributionMap();
  });
}

std::vector<const DebugSnapshot::Hierarchy*>
GatherHierarchies_(const std::vector<DebugSnapshot::Hierarchy>& hierarchies,
                   ::amrex::IndexType location) {
  std::vector<const DebugSnapshot::Hierarchy*> gathered{};
  for (const DebugSnapshot::Hierarchy& hierarchy : hierarchies) {
    if (hierarchy[0].boxArray().ixType() == location) {
      gathered.push_back(&hierarchy);
    }
  }
  return gathered;
}

std::vector<DebugSnapshot::ComponentNames>
GatherNames_(const std::vector<DebugSnapshot::ComponentNames>& names,
             const std::vector<DebugSnapshot::Hierarchy>& hierarchies,
             ::amrex::IndexType location) {
  std::vector<DebugSnapshot::ComponentNames> gathered{};
  auto components = names.begin();
  for (const DebugSnapshot::Hierarchy& hierarchy : hierarchies) {
    if (hierarchy[0].boxArray().ixType() == location) {
      gathered.push_back(*components);
    }
    ++components;
  }
  return gathered;
}

DebugSnapshot::ComponentNames
Select_(const std::vector<DebugSnapshot::ComponentNames>& names,
        const std::vector<int>& partition) {
  using namespace std::literals;
  DebugSnapshot::ComponentNames selected;
  for (int i : partition) {
    selected.insert(selected.end(), names[i].begin(), names[i].end());
  }
  return selected;
}

// Apply face to cell average to a specified face_component on mf_faces and
// write its result into cell_component of mf_cells.
void
AverageFaceToCell(::amrex::MultiFab& mf_cells, int cell_component,
                  const ::amrex::MultiFab& mf_faces, int face_component) {
  ::amrex::IndexType facetype = mf_faces.ixType();
  Direction dir;
  if (facetype == ::amrex::IndexType(::amrex::IntVect(AMREX_D_DECL(1, 0, 0)))) {
    dir = Direction::X;
  } else if (facetype == ::amrex::IndexType(::amrex::IntVect(AMREX_D_DECL(0, 1, 0)))) {
    dir = Direction::Y;
  } else if (facetype == ::amrex::IndexType(::amrex::IntVect(AMREX_D_DECL(0, 0, 1)))) {
    dir = Direction::Z;
  } else {
    throw std::runtime_error("Invalid Index type of face variable!");
  }
  const int dir_v = static_cast<int>(dir);
  fub::amrex::ForEachFab(mf_cells, [&](const ::amrex::MFIter& mfi) {
    ::amrex::Box cell_box = mfi.tilebox();
    ::amrex::Box face_box = ::amrex::surroundingNodes(cell_box, dir_v);
    auto cells = MakePatchDataView(mf_cells[mfi], cell_component, cell_box);
    auto faces = MakePatchDataView(mf_faces[mfi], face_component, face_box);
    fub::ForEachIndex(cells.Box(), [&](auto... is) {
      std::array<std::ptrdiff_t, AMREX_SPACEDIM> cell{is...};
      std::array<std::ptrdiff_t, AMREX_SPACEDIM> face_left{is...};
      std::array<std::ptrdiff_t, AMREX_SPACEDIM> face_right = Shift(face_left, dir, 1);
      cells(cell) = 0.5 * faces(face_left) + 0.5 * faces(face_right);
    });
  });

}


::amrex::IntVect GetRefRatio_(::amrex::Geometry geom0, ::amrex::Geometry geom1) {
  ::amrex::Box dom0 = geom0.Domain();
  ::amrex::Box dom1 = geom1.Domain();

  ::amrex::IntVect size0 = dom0.size();
  ::amrex::IntVect size1 = dom1.size();

  if (size1 < size0) {
    throw std::runtime_error("geom0 must be coarser than geom1!");
  }

  return size1 / size0;
}

} // namespace

/// \brief Returns all the hierarchies which are stored via SaveData
const std::vector<DebugSnapshot::Hierarchy>& DebugSnapshot::GetHierarchies() const
    noexcept {
  return saved_hierarchies_;
}

/// \brief Returns all the component names which are stored via SaveData
const std::vector<DebugSnapshot::ComponentNames>& DebugSnapshot::GetNames() const
    noexcept {
  return names_per_hierarchy_;
}

void DebugSnapshot::SaveData(const ::amrex::MultiFab& mf,
                            const std::string& name,
                            ::amrex::SrcComp component) {
    Hierarchy single_level_hierarchy{};
    ::amrex::BoxArray ba = mf.boxArray();
    ::amrex::DistributionMapping dm = mf.DistributionMap();
    ::amrex::MultiFab& dest = single_level_hierarchy.emplace_back(ba, dm, 1, 0);
    ::amrex::MultiFab::Copy(dest, mf, component.i, 0, 1, 0);
    saved_hierarchies_.push_back(std::move(single_level_hierarchy));
    names_per_hierarchy_.push_back(ComponentNames{name});
}

void DebugSnapshot::SaveData(const ::amrex::MultiFab& mf,
                            const ComponentNames& names,
                            ::amrex::SrcComp first_component) {
    Hierarchy single_level_hierarchy{};
    ::amrex::BoxArray ba = mf.boxArray();
    ::amrex::DistributionMapping dm = mf.DistributionMap();
    int n_components = static_cast<int>(names.size());
    ::amrex::MultiFab& dest = single_level_hierarchy.emplace_back(ba, dm, n_components, 0);
    ::amrex::MultiFab::Copy(dest, mf, first_component.i, 0, n_components, 0);
    saved_hierarchies_.push_back(std::move(single_level_hierarchy));
    names_per_hierarchy_.push_back(ComponentNames{names});
}

void DebugSnapshot::SaveData(const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
              const std::string& name,
              ::amrex::SrcComp component) {
  SaveData(hierarchy, ComponentNames{name}, component);
}

void DebugSnapshot::SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
              const std::string& name,
              ::amrex::SrcComp component) {
  SaveData(::amrex::GetVecOfConstPtrs(hierarchy), ComponentNames{name}, component);
}

void DebugSnapshot::SaveData(
    const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
    const ComponentNames& names, ::amrex::SrcComp first_component) {
    std::size_t nlevels = hierarchy.size();
    Hierarchy& multi_level_hierarchy = saved_hierarchies_.emplace_back();
    multi_level_hierarchy.reserve(nlevels);
    for (const ::amrex::MultiFab* mf_pointer : hierarchy) {
      FUB_ASSERT(mf_pointer != nullptr);
      const ::amrex::MultiFab& mf = *mf_pointer;
      ::amrex::BoxArray ba = mf.boxArray();
      ::amrex::DistributionMapping dm = mf.DistributionMap();
      ::amrex::MultiFab& dest =
          multi_level_hierarchy.emplace_back(ba, dm, names.size(), 0);
      ::amrex::MultiFab::Copy(dest, mf, first_component.i, 0, names.size(), 0);
    }
    names_per_hierarchy_.push_back(names);
}

void DebugSnapshot::SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                            const ComponentNames& names,
                            ::amrex::SrcComp first_component) {
    SaveData(::amrex::GetVecOfConstPtrs(hierarchy), names, first_component);
}

void DebugSnapshot::MakeUniqueComponentNames() {
  using namespace std::literals;

  std::size_t h = 0;
  for (DebugSnapshot::ComponentNames& names : names_per_hierarchy_) {

    auto last = names.end();
    for (int k = 0; k < names.size(); ++k) {
      const std::string candidate = names[k];
      int counter = 0;
      // search for candidate in current hierarchy
      auto first = names.begin() + k + 1;
      auto pos = std::find(first, last, candidate);
      if (last != pos) {
        names[k] += "_0"s;
        counter += 1;
        do {
          *pos = fmt::format("{}_{}", *pos, counter);
          counter += 1;
          pos = std::find(first, last, candidate);
        } while (last != pos);
      }
      // search for candidate in other hierarchies
      auto nexthier = names_per_hierarchy_.begin() + h + 1;
      auto lasthier = names_per_hierarchy_.end();
      for (auto p = nexthier; p != lasthier; ++p) {
        DebugSnapshot::ComponentNames& hnames = *p;
        auto first = hnames.begin();
        auto last = hnames.end();
        auto pos = std::find(first, last, candidate);
        if (last != pos) {
          if (counter == 0) {
            names[k] += "_0"s;
            counter += 1;
          }
          do {
            *pos = fmt::format("{}_{}", *pos, counter);
            counter += 1;
            pos = std::find(first, last, candidate);
          } while (last != pos);
        }
      }
    }
    ++h;
  }
}

std::vector<std::pair<DebugSnapshot::Hierarchy, DebugSnapshot::ComponentNames>>
DebugSnapshot::GatherFields(::amrex::IndexType location) const {
  using namespace ::amrex;
  std::vector<std::vector<BoxArray>> box_arrays =
      GatherBoxArrays_(saved_hierarchies_, location);
  std::vector<std::vector<DistributionMapping>> distribution_maps =
      GatherDistributionMaps_(saved_hierarchies_, location);
  std::vector<const DebugSnapshot::Hierarchy*> filtered_hierarchies =
      GatherHierarchies_(saved_hierarchies_, location);
  std::vector<std::vector<int>> partitions = Partitions_(box_arrays);
  std::vector<std::pair<DebugSnapshot::Hierarchy, DebugSnapshot::ComponentNames>>
      hiers_with_names{};
  for (const std::vector<int>& partition : partitions) {
    std::vector<ComponentNames> all_components =
        GatherNames_(names_per_hierarchy_, saved_hierarchies_, location);
    ComponentNames components = Select_(all_components, partition);
    int n_components = static_cast<int>(components.size());
    DebugSnapshot::Hierarchy hierarchy{};
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
    hiers_with_names.emplace_back(std::move(hierarchy), std::move(components));
  }
  return hiers_with_names;
}

void DebugSnapshot::ClearAll() {
  names_per_hierarchy_.clear();
  saved_hierarchies_.clear();
}


void DebugSnapshotProxy::SaveData(const ::amrex::MultiFab& mf,
                            const std::string& name,
                            ::amrex::SrcComp component) {
  if (m_snapshot) {
    m_snapshot->SaveData(mf, name, component);
  }
}

void DebugSnapshotProxy::SaveData(const ::amrex::MultiFab& mf,
                            const DebugSnapshot::ComponentNames& names,
                            ::amrex::SrcComp first_component) {
  if (m_snapshot) {
    m_snapshot->SaveData(mf, names, first_component);
  }
}

void DebugSnapshotProxy::SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                            const std::string& name,
                            ::amrex::SrcComp first_component) {
  if (m_snapshot) {
    m_snapshot->SaveData(hierarchy, name, first_component);
  }
}

void DebugSnapshotProxy::SaveData(
    const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
    const std::string& name, ::amrex::SrcComp first_component) {
  if (m_snapshot) {
    m_snapshot->SaveData(hierarchy, name, first_component);
  }
}

void DebugSnapshotProxy::SaveData(
    const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
    const DebugSnapshot::ComponentNames& names, ::amrex::SrcComp first_component) {
  if (m_snapshot) {
    m_snapshot->SaveData(hierarchy, names, first_component);
  }
}

void DebugSnapshotProxy::SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                            const DebugSnapshot::ComponentNames& names,
                            ::amrex::SrcComp first_component) {
  if (m_snapshot) {
    m_snapshot->SaveData(hierarchy, names, first_component);
  }
}


/// \brief Returns all the snapshots which are stored via FlushData
std::list<DebugSnapshot>& DebugStorage::GetSnapshots()
    noexcept {
  return saved_snapshots_;
}

DebugSnapshotProxy DebugStorage::AddSnapshot(const std::string& snapshot_directory) {
  if (is_enabled_) {
    DebugSnapshot& snapshot = saved_snapshots_.emplace_back();
    snapshot.SetSnapshotDirectory(snapshot_directory);
    DebugSnapshotProxy snapshotproxy(snapshot);
    return snapshotproxy;
  }
  return DebugSnapshotProxy{};
}

DebugOutput::DebugOutput(const ProgramOptions& opts,
                         const std::shared_ptr<DebugStorage>& storage)
    : OutputAtFrequencyOrInterval(opts) {
  OutputAtFrequencyOrInterval::frequencies_ = std::vector<std::ptrdiff_t>{1LL};
  directory_ = GetOptionOr(opts, "directory", directory_);
  storage->Enable();
}

void DebugOutput::operator()(const GriddingAlgorithm& grid) {
  using namespace std::literals;
  DebugStorage& storage = *grid.GetPatchHierarchy().GetDebugStorage();
  std::list<DebugSnapshot>& saved_snapshots = storage.GetSnapshots();

  for (DebugSnapshot& snapshot : saved_snapshots) {
    snapshot.MakeUniqueComponentNames();

    const std::vector<DebugSnapshot::Hierarchy>& all_hierarchies =
        snapshot.GetHierarchies();
    const std::vector<DebugSnapshot::ComponentNames>& all_names =
        snapshot.GetNames();

    // average node/face data to cell data
    std::vector<DebugSnapshot::Hierarchy> avg_hierarchies{};
    std::vector<DebugSnapshot::ComponentNames> names_avg_hierarchies{};

    auto hier_fields = all_names.begin();
    for (const DebugSnapshot::Hierarchy& hierarchy : all_hierarchies) {
      if (hierarchy[0].ixType() != ::amrex::IndexType::TheCellType()) {
        const int nlevels = hierarchy.size();
        const int nComp = hierarchy[0].nComp();
        DebugSnapshot::Hierarchy cell_average_hierarchy(nlevels);
        DebugSnapshot::ComponentNames names_avg(*hier_fields);

        if (hierarchy[0].ixType() == ::amrex::IndexType::TheNodeType()) {
          for (int level = 0; level < nlevels; ++level) {
            ::amrex::BoxArray ba = hierarchy[level].boxArray();
            ba.enclosedCells();
            ::amrex::DistributionMapping dm = hierarchy[level].DistributionMap();
            cell_average_hierarchy[level].define(ba, dm, nComp, 0);

            ::amrex::average_node_to_cellcenter(cell_average_hierarchy[level], 0,
                                                hierarchy[level], 0, nComp, 0);
          }
          for (std::string& vname : names_avg) {
            vname += "_nd2cellavg"s;
          }
        } else {
          for (int level = 0; level < nlevels; ++level) {
            ::amrex::BoxArray ba = hierarchy[level].boxArray();
            ba.enclosedCells();
            ::amrex::DistributionMapping dm = hierarchy[level].DistributionMap();
            cell_average_hierarchy[level].define(ba, dm, nComp, 0);

            for (int Comp = 0; Comp < nComp; ++Comp) {
              AverageFaceToCell(cell_average_hierarchy[level], Comp,
                      hierarchy[level], Comp);
            }
          }
          for (std::string& vname : names_avg) {
            vname += "_fc2cellavg"s;
          }
        }

        avg_hierarchies.emplace_back(std::move(cell_average_hierarchy));
        names_avg_hierarchies.emplace_back(std::move(names_avg));
      }
      ++hier_fields;
    }

    auto avg_fields = names_avg_hierarchies.begin();
    for (const DebugSnapshot::Hierarchy& hierarchy : avg_hierarchies) {
      DebugSnapshot::ComponentNames names_avg(*avg_fields);
      snapshot.SaveData(hierarchy, names_avg);
      ++avg_fields;
    }

    avg_hierarchies.clear();
    names_avg_hierarchies.clear();

    // save data from snapshot by partition (same patch hierarchy)
    std::vector<std::pair<DebugSnapshot::Hierarchy, DebugSnapshot::ComponentNames>>
        hiers_with_names =
            snapshot.GatherFields(::amrex::IndexType::TheCellType());
    const std::string snapshot_directory = snapshot.GetSnapshotDirectory();

    int partition_counter = 0;
    for (auto&& [hierarchy, names] : hiers_with_names) {
      const std::string plotfilename =
          fmt::format("{}/{}/partition_{}_plt{:09}", directory_, snapshot_directory, partition_counter,
                      grid.GetPatchHierarchy().GetCycles());
      const std::size_t size = hierarchy.size();
      const double time_point = grid.GetTimePoint().count();
      ::amrex::Vector<const ::amrex::MultiFab*> mf(size);
      ::amrex::Vector<::amrex::Geometry> geoms(size);
      ::amrex::Vector<int> level_steps(size);
      ::amrex::Vector<::amrex::IntVect> ref_ratio(size-1);
      std::vector<::amrex::BoxArray> cell_box_array(size);

      for (std::size_t i = 0; i < size; ++i) {
        mf[i] = &hierarchy[i];
        const int ii = static_cast<int>(i);
        geoms[i] = grid.GetPatchHierarchy().GetGeometry(ii);
        level_steps[i] = static_cast<int>(grid.GetPatchHierarchy().GetCycles(ii));
        cell_box_array[i] = hierarchy[i].boxArray();
      }

      for (std::size_t i = 0; i < size-1; ++i) {
        ref_ratio[i] = GetRefRatio_(geoms[i], geoms[i+1]);
      }

      ::amrex::Vector<std::string> vnames(names.begin(), names.end());
      ::amrex::WriteMultiLevelPlotfile(
          plotfilename, size, mf, vnames, geoms, time_point, level_steps,
          ref_ratio, "HyperCLaw-V1.1", "Level_", "Cell", {"raw_fields"});

      ::amrex::VisMF::SetHeaderVersion(::amrex::VisMF::Header::Version_v1);

      // write nodal/face hierarchies as raw data
      auto hier_fields = all_names.begin();

      for (const DebugSnapshot::Hierarchy& hierarchy : all_hierarchies) {
        if (hierarchy[0].ixType() != ::amrex::IndexType::TheCellType()) {
          const std::size_t size = hierarchy.size();
          const int nComp = hierarchy[0].nComp();
          std::vector<std::string> names(*hier_fields);

          std::vector<::amrex::BoxArray> noncell_box_array(size);
          for (std::size_t i = 0; i < size; ++i) {
            noncell_box_array[i] = hierarchy[i].boxArray();
            noncell_box_array[i].enclosedCells();
          }
          // check if we are working on the correct partition
          if (cell_box_array == noncell_box_array) {
            for (std::size_t i = 0; i < size; ++i) {
              const int ii = static_cast<int>(i);
                ::amrex::BoxArray ba = hierarchy[i].boxArray();
                ::amrex::DistributionMapping dm = hierarchy[i].DistributionMap();
              const std::string raw_pltname = fmt::format("{}/raw_fields", plotfilename);
              const std::string level_prefix = "Level_";
              for (int Comp = 0; Comp < nComp; ++Comp) {
                const std::string full_prefix =
                ::amrex::MultiFabFileFullPrefix(ii, raw_pltname, level_prefix, names[Comp]);
                ::amrex::MultiFab data(ba, dm, 1, 0);
                ::amrex::MultiFab::Copy(data, hierarchy[i], Comp, 0, 1, 0);
                ::amrex::VisMF::Write(data, full_prefix);
              }
            }
          }
        }
        ++hier_fields;
      }
      partition_counter += 1;
    }
  }
  storage.ClearAll();
}

} // namespace fub::amrex
