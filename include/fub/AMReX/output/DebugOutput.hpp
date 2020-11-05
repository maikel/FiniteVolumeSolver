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

#ifndef FUB_AMREX_DEBUG_OUTPUT_HPP
#define FUB_AMREX_DEBUG_OUTPUT_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

namespace fub::amrex {

/// \brief This class stores debug data for a debug output for a single
/// hierarchy state.
///
/// This class will be created by \a DebugStorage and modified through \a
/// DebugSnapshotProxy.
///
/// \see DebugStorage
/// \see DebugSnapshotProxy
class DebugSnapshot {
public:
  using Hierarchy = ::amrex::Vector<::amrex::MultiFab>;
  using GeomHierarchy = ::amrex::Vector<::amrex::Geometry>;
  using ComponentNames = ::amrex::Vector<std::string>;

  /// \brief Initializes an empty snapshot
  DebugSnapshot() = default;

  /// @{
  /// \brief Saves a current hierarchy state with given component names
  ///
  /// The actual output will be handled by the DebugOutput class and will
  /// usually happen at a later time point.
  void SaveData(const ::amrex::MultiFab& mf, const std::string& name,
                const ::amrex::Geometry& geom,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::MultiFab& mf, const ComponentNames& names,
                const ::amrex::Geometry& geom,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
                const std::string& name,
                const ::amrex::Vector<const ::amrex::Geometry*>& geomhier,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                const std::string& name,
                const ::amrex::Vector<::amrex::Geometry>& geomhier,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
                const ComponentNames& names,
                const ::amrex::Vector<const ::amrex::Geometry*>& geomhier,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                const ComponentNames& names,
                const ::amrex::Vector<::amrex::Geometry>& geomhier,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));
  /// @}

  /// \brief Deletes all currently stored data
  void ClearAll();

  /// \brief Returns all the hierarchies which are stored via SaveData
  const std::vector<Hierarchy>& GetHierarchies() const noexcept;

  /// \brief Returns all the component names which are stored via SaveData
  const std::vector<ComponentNames>& GetNames() const noexcept;

  /// \brief Returns all the associated geometries which are stored via SaveData
  const std::vector<GeomHierarchy>& GetGeometries() const noexcept;

  /// \brief Changes component names that appear more than once by appending
  /// assending numbers to each. This also makes names unique, which appear
  /// more than once in one hierarchy.
  void MakeUniqueComponentNames();

  /// \brief Collects all hierachies and associated names which are stored on
  /// the specified location.
  ///
  /// This function will join all compatible hierachies to a common big one if
  /// possible.
  std::vector<std::tuple<Hierarchy, ComponentNames, GeomHierarchy>>
  GatherFields(::amrex::IndexType location) const;

  /// \brief Sets the directory in which the output will be saved to.
  void SetSnapshotDirectory(const std::string& snapshot_directory) {
    snapshot_directory_ = snapshot_directory;
  }

  /// \brief Returns the directory in which the output will be saved to.
  const std::string GetSnapshotDirectory() const { return snapshot_directory_; }

private:
  /// \brief Each SaveData will append a hierarchy
  std::vector<Hierarchy> saved_hierarchies_{};

  /// \brief Each SaveData will append a list of names on the hierachy
  std::vector<ComponentNames> names_per_hierarchy_{};

  /// \brief Each SaveData will append a hierarchy of corresponding geometries
  std::vector<GeomHierarchy> saved_geometries_{};

  /// \brief output directory of storage;
  std::string snapshot_directory_{};
};

/// \brief This class is a possibly empty handle to a existing
/// DebugSnapshotProxy.
///
/// DebugSnapshotProxy maybe in a valid or invalid
class DebugSnapshotProxy {
public:
  /// \brief Initializes an empty proxy snapshot
  DebugSnapshotProxy() = default;

  /// \brief Initializes a proxy snapshot for a real snapshot
  DebugSnapshotProxy(DebugSnapshot& snapshot) : snapshot_(&snapshot){};

  /// @{
  /// \brief Saves a current hierarchy state with given component names
  ///
  /// The actual output will be handled by the DebugOutput class and will
  /// usually happen at a later time point.
  void SaveData(const ::amrex::MultiFab& mf, const std::string& name,
                const ::amrex::Geometry& geom,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::MultiFab& mf,
                const DebugSnapshot::ComponentNames& names,
                const ::amrex::Geometry& geom,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
                const std::string& name,
                const ::amrex::Vector<const ::amrex::Geometry*>& geomhier,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                const std::string& name,
                const ::amrex::Vector<::amrex::Geometry>& geomhier,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
                const DebugSnapshot::ComponentNames& names,
                const ::amrex::Vector<const ::amrex::Geometry*>& geomhier,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                const DebugSnapshot::ComponentNames& names,
                const ::amrex::Vector<::amrex::Geometry>& geomhier,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));
  /// @}

  /// \brief Returns true if this object is valid.
  explicit operator bool() const noexcept { return snapshot_ == nullptr; }

  /// \brief Returns true if both proxy objects point to the same snapshot.
  friend inline bool operator==(const DebugSnapshotProxy& p1,
                                const DebugSnapshotProxy& p2) noexcept {
    return p1.snapshot_ == p2.snapshot_;
  }

private:
  /// \brief Pointer to the real snapshot
  DebugSnapshot* snapshot_{nullptr};
};

/// \brief Returns true if both proxy objects point to different snapshots.
inline bool operator!=(const DebugSnapshotProxy& p1,
                       const DebugSnapshotProxy& p2) noexcept {
  return !(p1 == p2);
}

/// \brief This class stores a list of snapshots and returns proxy objects to
/// those.
///
/// As long as this class is in a disabled mode the returned proxy objects are
/// invalid and subsequent calls to SaveData will have no effect. Only if this
/// the storage is enabled debug data will be stored.
class DebugStorage {
public:
  /// \brief Initializes the storage
  DebugStorage() = default;

  /// \brief Checks if debug storage is enabled
  bool IsEnabled() const noexcept { return is_enabled_; }

  void Enable() { is_enabled_ = true; }
  void Disable() { is_enabled_ = false; }

  /// \brief Adds snapshot with given directory name to storage.
  DebugSnapshotProxy AddSnapshot(const std::string& snapshot_directory);

  /// \brief Writes data to disk and clears all data present in the storage.
  void FlushData(const std::string& directory, int cycle = -1,
                 Duration time_point = Duration(0.0));

  /// \brief Deletes all currently stored snapshots in the debug storage
  void ClearAll() { saved_snapshots_.clear(); }

private:
  /// \brief List of snapshots
  std::list<DebugSnapshot> saved_snapshots_;

  /// \brief If this is false no data will be saved.
  bool is_enabled_{false};

  int cycle_{-1};
};

/// \brief This output method enables a debug storage and manages its output in
/// every time step.
class DebugOutput : public OutputAtFrequencyOrInterval<GriddingAlgorithm> {
public:
  /// \brief Read program options from opts and enable the storage.
  explicit DebugOutput(const ProgramOptions& opts,
                       const std::shared_ptr<DebugStorage>& storage);

  /// \brief Write out the debug storage on the grid.
  void operator()(const GriddingAlgorithm& grid) override;

private:
  /// \brief This is the base directory where the snapshots will be output to.
  std::string directory_;
};

namespace cutcell {
class DebugOutput : public OutputAtFrequencyOrInterval<GriddingAlgorithm> {
public:
  /// \brief Read program options from opts and enable the storage.
  explicit DebugOutput(const ProgramOptions& opts,
                       const std::shared_ptr<DebugStorage>& storage);

  /// \brief Write out the debug storage on the grid.
  void operator()(const GriddingAlgorithm& grid) override;

private:
  /// \brief This is the base directory where the snapshots will be output to.
  std::string directory_;
};
}

} // namespace fub::amrex

#endif
