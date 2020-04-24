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
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/ext/Log.hpp"

#include <algorithm>
#include <limits>
#include <numeric>

namespace fub::amrex {

namespace {

template <class T>
::amrex::Vector<const T*> GetVecOfConstPtrs(const ::amrex::Vector<T>& a) {
  ::amrex::Vector<const T*> r;
  r.reserve(static_cast<std::size_t>(a.size()));
  for (const auto& x : a) {
    r.push_back(&x);
  }
  return r;
}

struct equal_to {
  int i;
  const std::vector<std::vector<::amrex::BoxArray>>* bas_;
  bool operator()(int j) {
    return (*bas_)[static_cast<std::size_t>(i)] ==
           (*bas_)[static_cast<std::size_t>(j)];
  }
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
        std::vector<T> projected(static_cast<std::size_t>(hierarchy.size()));
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

std::vector<std::vector<::amrex::DistributionMapping>> GatherDistributionMaps_(
    const std::vector<DebugSnapshot::Hierarchy>& hierarchies,
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

std::vector<const DebugSnapshot::GeomHierarchy*>
GatherGeometries_(const std::vector<DebugSnapshot::GeomHierarchy>& geoms,
                  const std::vector<DebugSnapshot::Hierarchy>& hierarchies,
                  ::amrex::IndexType location) {
  std::vector<const DebugSnapshot::GeomHierarchy*> gathered{};
  auto geom = geoms.begin();
  for (const DebugSnapshot::Hierarchy& hierarchy : hierarchies) {
    if (hierarchy[0].boxArray().ixType() == location) {
      gathered.emplace_back(&(*geom));
    }
    ++geom;
  }
  return gathered;
}

DebugSnapshot::ComponentNames
Select_(const std::vector<DebugSnapshot::ComponentNames>& names,
        const std::vector<int>& partition) {
  using namespace std::literals;
  DebugSnapshot::ComponentNames selected;
  for (int i : partition) {
    const auto s = static_cast<std::size_t>(i);
    selected.insert(selected.end(), names[s].begin(), names[s].end());
  }
  return selected;
}

// Apply face to cell average to a specified face_component on mf_faces and
// write its result into cell_component of mf_cells.
void AverageFaceToCell(::amrex::MultiFab& mf_cells, int cell_component,
                       const ::amrex::MultiFab& mf_faces, int face_component) {
  ::amrex::IndexType facetype = mf_faces.ixType();
  Direction dir;
  if (facetype == ::amrex::IndexType(::amrex::IntVect(AMREX_D_DECL(1, 0, 0)))) {
    dir = Direction::X;
  } else if (facetype ==
             ::amrex::IndexType(::amrex::IntVect(AMREX_D_DECL(0, 1, 0)))) {
    dir = Direction::Y;
  } else if (facetype ==
             ::amrex::IndexType(::amrex::IntVect(AMREX_D_DECL(0, 0, 1)))) {
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
      std::array<std::ptrdiff_t, AMREX_SPACEDIM> face_right =
          Shift(face_left, dir, 1);
      cells(cell) = 0.5 * faces(face_left) + 0.5 * faces(face_right);
    });
  });
}

::amrex::IntVect GetRefRatio_(::amrex::Geometry geom0,
                              ::amrex::Geometry geom1) {
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
const std::vector<DebugSnapshot::Hierarchy>&
DebugSnapshot::GetHierarchies() const noexcept {
  return saved_hierarchies_;
}

/// \brief Returns all the component names which are stored via SaveData
const std::vector<DebugSnapshot::ComponentNames>&
DebugSnapshot::GetNames() const noexcept {
  return names_per_hierarchy_;
}

/// \brief Returns all the associated geometries which are stored via SaveData
const std::vector<DebugSnapshot::GeomHierarchy>&
DebugSnapshot::GetGeometries() const noexcept {
  return saved_geometries_;
}

void DebugSnapshot::SaveData(const ::amrex::MultiFab& mf,
                             const std::string& name,
                             const ::amrex::Geometry& geom,
                             ::amrex::SrcComp component) {
  Hierarchy single_level_hierarchy{};
  ::amrex::BoxArray ba = mf.boxArray();
  ::amrex::DistributionMapping dm = mf.DistributionMap();
  ::amrex::MultiFab& dest = single_level_hierarchy.emplace_back(ba, dm, 1, 0);
  ::amrex::MultiFab::Copy(dest, mf, component.i, 0, 1, 0);
  saved_hierarchies_.push_back(std::move(single_level_hierarchy));
  names_per_hierarchy_.push_back(ComponentNames{name});
  saved_geometries_.push_back(GeomHierarchy{geom});
}

void DebugSnapshot::SaveData(const ::amrex::MultiFab& mf,
                             const ComponentNames& names,
                             const ::amrex::Geometry& geom,
                             ::amrex::SrcComp first_component) {
  Hierarchy single_level_hierarchy{};
  ::amrex::BoxArray ba = mf.boxArray();
  ::amrex::DistributionMapping dm = mf.DistributionMap();
  int n_components = static_cast<int>(names.size());
  ::amrex::MultiFab& dest =
      single_level_hierarchy.emplace_back(ba, dm, n_components, 0);
  ::amrex::MultiFab::Copy(dest, mf, first_component.i, 0, n_components, 0);
  saved_hierarchies_.push_back(std::move(single_level_hierarchy));
  names_per_hierarchy_.push_back(ComponentNames{names});
  saved_geometries_.push_back(GeomHierarchy{geom});
}

void DebugSnapshot::SaveData(
    const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
    const std::string& name,
    const ::amrex::Vector<const ::amrex::Geometry*>& geomhier,
    ::amrex::SrcComp component) {
  SaveData(hierarchy, ComponentNames{name}, geomhier, component);
}

void DebugSnapshot::SaveData(
    const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
    const std::string& name, const ::amrex::Vector<::amrex::Geometry>& geomhier,
    ::amrex::SrcComp component) {
  SaveData(::amrex::GetVecOfConstPtrs(hierarchy), ComponentNames{name},
           GetVecOfConstPtrs(geomhier), component);
}

void DebugSnapshot::SaveData(
    const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
    const ComponentNames& names,
    const ::amrex::Vector<const ::amrex::Geometry*>& geomhier,
    ::amrex::SrcComp first_component) {
  const auto nlevels = static_cast<std::size_t>(hierarchy.size());
  Hierarchy& multi_level_hierarchy = saved_hierarchies_.emplace_back();
  multi_level_hierarchy.reserve(nlevels);
  for (const ::amrex::MultiFab* mf_pointer : hierarchy) {
    FUB_ASSERT(mf_pointer != nullptr);
    const ::amrex::MultiFab& mf = *mf_pointer;
    ::amrex::BoxArray ba = mf.boxArray();
    ::amrex::DistributionMapping dm = mf.DistributionMap();
    ::amrex::MultiFab& dest =
        multi_level_hierarchy.emplace_back(ba, dm, names.size(), 0);
    const auto n_names = static_cast<int>(names.size());
    ::amrex::MultiFab::Copy(dest, mf, first_component.i, 0, n_names, 0);
  }
  names_per_hierarchy_.push_back(names);
  GeomHierarchy& multi_level_geometry = saved_geometries_.emplace_back();
  multi_level_geometry.reserve(nlevels);
  for (const ::amrex::Geometry* geom_pointer : geomhier) {
    multi_level_geometry.emplace_back(*geom_pointer);
  }
}

void DebugSnapshot::SaveData(
    const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
    const ComponentNames& names,
    const ::amrex::Vector<::amrex::Geometry>& geomhier,
    ::amrex::SrcComp first_component) {
  SaveData(::amrex::GetVecOfConstPtrs(hierarchy), names,
           GetVecOfConstPtrs(geomhier), first_component);
}

void DebugSnapshot::MakeUniqueComponentNames() {
  using namespace std::literals;

  std::size_t h = 0;
  for (DebugSnapshot::ComponentNames& names : names_per_hierarchy_) {
    auto last = names.end();
    const auto n_names = static_cast<std::size_t>(names.size());
    for (std::size_t k = 0; k < n_names; ++k) {
      const std::string candidate = names[k];
      int counter = 0;
      // search for candidate in current hierarchy
      auto first = names.begin() + static_cast<std::ptrdiff_t>(k + 1);
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
      auto nexthier =
          names_per_hierarchy_.begin() + static_cast<std::ptrdiff_t>(h + 1);
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

std::vector<std::tuple<DebugSnapshot::Hierarchy, DebugSnapshot::ComponentNames,
                       DebugSnapshot::GeomHierarchy>>
DebugSnapshot::GatherFields(::amrex::IndexType location) const {
  using namespace ::amrex;
  std::vector<std::vector<BoxArray>> box_arrays =
      GatherBoxArrays_(saved_hierarchies_, location);
  std::vector<std::vector<DistributionMapping>> distribution_maps =
      GatherDistributionMaps_(saved_hierarchies_, location);
  std::vector<const DebugSnapshot::Hierarchy*> filtered_hierarchies =
      GatherHierarchies_(saved_hierarchies_, location);
  std::vector<const DebugSnapshot::GeomHierarchy*> filtered_geoms =
      GatherGeometries_(saved_geometries_, saved_hierarchies_, location);
  std::vector<std::vector<int>> partitions = Partitions_(box_arrays);
  std::vector<
      std::tuple<DebugSnapshot::Hierarchy, DebugSnapshot::ComponentNames,
                 DebugSnapshot::GeomHierarchy>>
      hiers_tuple{};
  for (const std::vector<int>& partition : partitions) {
    std::vector<ComponentNames> all_components =
        GatherNames_(names_per_hierarchy_, saved_hierarchies_, location);
    ComponentNames components = Select_(all_components, partition);
    const auto p = static_cast<std::size_t>(partition[0]);
    GeomHierarchy geoms = *filtered_geoms[p];
    int n_components = static_cast<int>(components.size());
    DebugSnapshot::Hierarchy hierarchy{};
    const std::vector<BoxArray>& box_array = box_arrays[p];
    const std::vector<DistributionMapping>& distribution = distribution_maps[p];
    for (std::size_t level = 0; level < box_array.size(); ++level) {
      MultiFab& dest = hierarchy.emplace_back(
          box_array[level], distribution[level], n_components, 0);
      int comp = 0;
      for (int i : partition) {
        const auto j = static_cast<std::size_t>(i);
        const MultiFab& src = (*filtered_hierarchies[j])[level];
        const int n = src.nComp();
        MultiFab::Copy(dest, src, 0, comp, n, 0);
        comp += n;
      }
      FUB_ASSERT(comp == n_components);
    }
    hiers_tuple.emplace_back(std::move(hierarchy), std::move(components),
                             std::move(geoms));
  }
  return hiers_tuple;
}

void DebugSnapshot::ClearAll() {
  saved_hierarchies_.clear();
  names_per_hierarchy_.clear();
  saved_geometries_.clear();
}

void DebugSnapshotProxy::SaveData(const ::amrex::MultiFab& mf,
                                  const std::string& name,
                                  const ::amrex::Geometry& geom,
                                  ::amrex::SrcComp component) {
  if (snapshot_) {
    snapshot_->SaveData(mf, name, geom, component);
  }
}

void DebugSnapshotProxy::SaveData(const ::amrex::MultiFab& mf,
                                  const DebugSnapshot::ComponentNames& names,
                                  const ::amrex::Geometry& geom,
                                  ::amrex::SrcComp first_component) {
  if (snapshot_) {
    snapshot_->SaveData(mf, names, geom, first_component);
  }
}

void DebugSnapshotProxy::SaveData(
    const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
    const std::string& name,
    const ::amrex::Vector<const ::amrex::Geometry*>& geomhier,
    ::amrex::SrcComp first_component) {
  if (snapshot_) {
    snapshot_->SaveData(hierarchy, name, geomhier, first_component);
  }
}

void DebugSnapshotProxy::SaveData(
    const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
    const std::string& name, const ::amrex::Vector<::amrex::Geometry>& geomhier,
    ::amrex::SrcComp first_component) {
  if (snapshot_) {
    snapshot_->SaveData(hierarchy, name, geomhier, first_component);
  }
}

void DebugSnapshotProxy::SaveData(
    const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
    const DebugSnapshot::ComponentNames& names,
    const ::amrex::Vector<const ::amrex::Geometry*>& geomhier,
    ::amrex::SrcComp first_component) {
  if (snapshot_) {
    snapshot_->SaveData(hierarchy, names, geomhier, first_component);
  }
}

void DebugSnapshotProxy::SaveData(
    const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
    const DebugSnapshot::ComponentNames& names,
    const ::amrex::Vector<::amrex::Geometry>& geomhier,
    ::amrex::SrcComp first_component) {
  if (snapshot_) {
    snapshot_->SaveData(hierarchy, names, geomhier, first_component);
  }
}

DebugSnapshotProxy
DebugStorage::AddSnapshot(const std::string& snapshot_directory) {
  if (is_enabled_) {
    DebugSnapshot& snapshot = saved_snapshots_.emplace_back();
    snapshot.SetSnapshotDirectory(snapshot_directory);
    DebugSnapshotProxy snapshotproxy(snapshot);
    return snapshotproxy;
  }
  return DebugSnapshotProxy{};
}

void DebugStorage::FlushData(const std::string& directory, int cycle,
                             Duration time_point) {
  using namespace std::literals;

  if (cycle < 0) {
    ++cycle_;
  } else {
    cycle_ = cycle;
  }

  for (DebugSnapshot& snapshot : saved_snapshots_) {
    snapshot.MakeUniqueComponentNames();

    const std::vector<DebugSnapshot::Hierarchy>& all_hierarchies =
        snapshot.GetHierarchies();
    const std::vector<DebugSnapshot::ComponentNames>& all_names =
        snapshot.GetNames();
    const std::vector<DebugSnapshot::GeomHierarchy>& all_geometries =
        snapshot.GetGeometries();

    // average node/face data to cell data
    std::vector<DebugSnapshot::Hierarchy> avg_hierarchies{};
    std::vector<DebugSnapshot::ComponentNames> names_avg_hierarchies{};
    std::vector<DebugSnapshot::GeomHierarchy> avg_geoms{};

    auto hier_fields = all_names.begin();
    auto hier_geoms = all_geometries.begin();
    for (const DebugSnapshot::Hierarchy& hierarchy : all_hierarchies) {
      if (hierarchy[0].ixType() != ::amrex::IndexType::TheCellType()) {
        const auto nlevels = static_cast<std::size_t>(hierarchy.size());
        const int nComp = hierarchy[0].nComp();
        DebugSnapshot::Hierarchy cell_average_hierarchy(nlevels);
        DebugSnapshot::ComponentNames names_avg(*hier_fields);
        DebugSnapshot::GeomHierarchy geoms(*hier_geoms);

        if (hierarchy[0].ixType() == ::amrex::IndexType::TheNodeType()) {
          for (std::size_t level = 0; level < nlevels; ++level) {
            ::amrex::BoxArray ba = hierarchy[level].boxArray();
            ba.enclosedCells();
            ::amrex::DistributionMapping dm =
                hierarchy[level].DistributionMap();
            cell_average_hierarchy[level].define(ba, dm, nComp, 0);

            ::amrex::average_node_to_cellcenter(cell_average_hierarchy[level],
                                                0, hierarchy[level], 0, nComp,
                                                0);
          }
          for (std::string& vname : names_avg) {
            vname += "_nd2cellavg"s;
          }
        } else {
          for (std::size_t level = 0; level < nlevels; ++level) {
            ::amrex::BoxArray ba = hierarchy[level].boxArray();
            ba.enclosedCells();
            ::amrex::DistributionMapping dm =
                hierarchy[level].DistributionMap();
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
        avg_geoms.emplace_back(std::move(geoms));
      }
      ++hier_fields;
      ++hier_geoms;
    }

    auto avg_fieldit = names_avg_hierarchies.begin();
    auto avg_geomit = avg_geoms.begin();
    for (const DebugSnapshot::Hierarchy& hierarchy : avg_hierarchies) {
      DebugSnapshot::ComponentNames names_avg(*avg_fieldit);
      snapshot.SaveData(hierarchy, names_avg, *avg_geomit);
      ++avg_fieldit;
      ++avg_geomit;
    }

    avg_hierarchies.clear();
    names_avg_hierarchies.clear();
    avg_geoms.clear();

    // save data from snapshot by partition (same patch hierarchy)
    std::vector<
        std::tuple<DebugSnapshot::Hierarchy, DebugSnapshot::ComponentNames,
                   DebugSnapshot::GeomHierarchy>>
        hiers_tuple = snapshot.GatherFields(::amrex::IndexType::TheCellType());
    const std::string snapshot_directory = snapshot.GetSnapshotDirectory();

    int partition_counter = 0;
    for (auto&& [hierarchy, names, geoms] : hiers_tuple) {
      const std::string plotfilename =
          fmt::format("{}/{}/partition_{}_plt{:09}", directory,
                      snapshot_directory, partition_counter, cycle_);
      const auto size = static_cast<std::size_t>(hierarchy.size());
      ::amrex::Vector<const ::amrex::MultiFab*> mf(size);
      ::amrex::Vector<int> level_steps(size);
      ::amrex::Vector<::amrex::IntVect> ref_ratio(size - 1);
      std::vector<::amrex::BoxArray> cell_box_array(size);

      for (std::size_t i = 0; i < size; ++i) {
        mf[i] = &hierarchy[i];
        level_steps[i] = 0;
        //         const int ii = static_cast<int>(i);
        //         level_steps[i] =
        //             static_cast<int>(grid.GetPatchHierarchy().GetCycles(ii));
        cell_box_array[i] = hierarchy[i].boxArray();
      }

      for (std::size_t i = 0; i < size - 1; ++i) {
        ref_ratio[i] = GetRefRatio_(geoms[i], geoms[i + 1]);
      }

      ::amrex::Vector<std::string> vnames(names.begin(), names.end());
      ::amrex::WriteMultiLevelPlotfile(plotfilename, static_cast<int>(size), mf,
                                       vnames, geoms, time_point.count(),
                                       level_steps, ref_ratio, "HyperCLaw-V1.1",
                                       "Level_", "Cell", {"raw_fields"});

      ::amrex::VisMF::SetHeaderVersion(::amrex::VisMF::Header::Version_v1);

      // write nodal/face hierarchies as raw data
      auto hier_fields = all_names.begin();

      for (const DebugSnapshot::Hierarchy& hierarchy : all_hierarchies) {
        if (hierarchy[0].ixType() != ::amrex::IndexType::TheCellType()) {
          const std::size_t size = static_cast<std::size_t>(hierarchy.size());
          std::vector<std::string> names(*hier_fields);
          std::vector<::amrex::BoxArray> noncell_box_array(size);
          for (std::size_t i = 0; i < size; ++i) {
            noncell_box_array[i] = hierarchy[i].boxArray();
            noncell_box_array[i].enclosedCells();
          }
          // Check if we are working on the correct partition,
          // Note that a correct partition must exist, since we have created
          // a hierarchy with cell averaged data before.
          if (cell_box_array == noncell_box_array) {
            for (std::size_t i = 0; i < size; ++i) {
              const int ii = static_cast<int>(i);
              ::amrex::BoxArray ba = hierarchy[i].boxArray();
              ::amrex::DistributionMapping dm = hierarchy[i].DistributionMap();
              const std::string raw_pltname =
                  fmt::format("{}/raw_fields", plotfilename);
              const std::string level_prefix = "Level_";
              const auto nComp = static_cast<std::size_t>(hierarchy[0].nComp());
              for (std::size_t c = 0; c < nComp; ++c) {
                const std::string full_prefix = ::amrex::MultiFabFileFullPrefix(
                    ii, raw_pltname, level_prefix, names[c]);
                ::amrex::MultiFab data(ba, dm, 1, 0);
                ::amrex::MultiFab::Copy(data, hierarchy[i], static_cast<int>(c),
                                        0, 1, 0);
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
  ClearAll();
}

DebugOutput::DebugOutput(const ProgramOptions& opts,
                         const std::shared_ptr<DebugStorage>& storage)
    : OutputAtFrequencyOrInterval(opts) {
  OutputAtFrequencyOrInterval::frequencies_ = std::vector<std::ptrdiff_t>{1LL};
  directory_ = GetOptionOr(opts, "directory", directory_);
  storage->Enable();
  SeverityLogger log = GetInfoLogger();
  BOOST_LOG(log) << "DebugOutput configured:";
  BOOST_LOG(log) << fmt::format(" - directory = '{}'", directory_);
  OutputAtFrequencyOrInterval::Print(log);
}

void DebugOutput::operator()(const GriddingAlgorithm& grid) {
  DebugStorage& storage = *grid.GetPatchHierarchy().GetDebugStorage();
  const std::ptrdiff_t cycles = grid.GetPatchHierarchy().GetCycles();
  const auto int_max =
      static_cast<std::ptrdiff_t>(std::numeric_limits<int>::max());
  FUB_ASSERT(cycles < int_max);
  storage.FlushData(directory_, static_cast<int>(cycles), grid.GetTimePoint());
}

} // namespace fub::amrex
