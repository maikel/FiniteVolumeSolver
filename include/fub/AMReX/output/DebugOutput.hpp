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
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

namespace fub::amrex {

class DebugSnapshot {
public:
  using Hierarchy = ::amrex::Vector<::amrex::MultiFab>;
  using ComponentNames = ::amrex::Vector<std::string>;

  /// \brief Initializes an empty snapshot
  DebugSnapshot() = default;

  /// @{
  /// \brief Saves a current hierarchy state with given component names
  ///
  /// The actual output will be handled by the DebugOutput class and will
  /// usually happen at a later time point.
  void SaveData(const ::amrex::MultiFab& mf, const std::string& name,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::MultiFab& mf, const ComponentNames& names,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
                const std::string& name,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                const std::string& name,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
                const ComponentNames& names,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                const ComponentNames& names,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));
  /// @}

  /// \brief Deletes all currently stored data
  void ClearAll();

  /// \brief Returns all the hierarchies which are stored via SaveData
  const std::vector<Hierarchy>& GetHierarchies() const noexcept;

  /// \brief Returns all the component names which are stored via SaveData
  const std::vector<ComponentNames>& GetNames() const noexcept;

  /// \brief Changes component names that appear more than once by appending
  /// assending numbers to each. This also makes names unique, which appear
  /// more than once in one hierarchy.
  void MakeUniqueComponentNames();

  /// \brief Collects all hierachies and associated names which are stored on
  /// the specified location.
  ///
  /// This function will join all compatible hierachies to a common big one if
  /// possible.
  std::vector<std::pair<Hierarchy, ComponentNames>>
  GatherFields(::amrex::IndexType location) const;

  /// \brief
  void SetSnapshotDirectory(const std::string& snapshot_directory) {
    snapshot_directory_ = snapshot_directory;
  }

  /// \brief
  const std::string GetSnapshotDirectory() const {
    return snapshot_directory_;
  }

private:
  /// \brief Each SaveData will append a hierarchy
  std::vector<Hierarchy> saved_hierarchies_{};

  /// \brief Each SaveData will append a list of names on the hierachy
  std::vector<ComponentNames> names_per_hierarchy_{};

  /// \brief output directory of storage;
  std::string snapshot_directory_{};

};

class DebugSnapshotProxy {
public:
  /// \brief Initializes an empty proxy snapshot
  DebugSnapshotProxy() = default;

  /// \brief Initializes a proxy snapshot for a real snapshot
  DebugSnapshotProxy(DebugSnapshot& snapshot) : m_snapshot(&snapshot) {};

  /// @{
  /// \brief Saves a current hierarchy state with given component names
  ///
  /// The actual output will be handled by the DebugOutput class and will
  /// usually happen at a later time point.
  void SaveData(const ::amrex::MultiFab& mf, const std::string& name,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::MultiFab& mf, const DebugSnapshot::ComponentNames& names,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
                const std::string& name,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                const std::string& name,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<const ::amrex::MultiFab*>& hierarchy,
                const DebugSnapshot::ComponentNames& names,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::Vector<::amrex::MultiFab>& hierarchy,
                const DebugSnapshot::ComponentNames& names,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));
  /// @}

private:
  /// \brief Pointer to the real snapshot
  DebugSnapshot* m_snapshot{nullptr};

};

class DebugStorage {
public:

  /// \brief Initializes the storage
  DebugStorage() = default;

  /// \brief Checks if debug storage is enabled
  bool IsEnabled() const noexcept { return is_enabled_; }

  void Enable() { is_enabled_ = true; }
  void Disable() { is_enabled_ = false; }

  /// \brief Returns all the snapshots which are stored in the debug storage
  std::list<DebugSnapshot>& GetSnapshots() noexcept;

  /// \brief Adds snapshot with given directory name to storage.
  DebugSnapshotProxy AddSnapshot(const std::string& snapshot_directory);

  /// \brief Deletes all currently stored snapshots in the debug storage
  void ClearAll() { saved_snapshots_.clear(); }

private:
  /// \brief List of snapshots
  std::list<DebugSnapshot> saved_snapshots_;

  /// \brief If this is false no data will be saved.
  bool is_enabled_{false};
};

class DebugOutput : public OutputAtFrequencyOrInterval<GriddingAlgorithm> {
public:
  explicit DebugOutput(const ProgramOptions& opts,
                       const std::shared_ptr<DebugStorage>& storage);

  void operator()(const GriddingAlgorithm& grid) override;

private:
  std::string directory_;
};

} // namespace fub::amrex

#endif
