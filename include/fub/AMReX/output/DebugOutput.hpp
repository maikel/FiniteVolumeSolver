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

#ifndef FUB_AMREX_DEBUG_OUTPUT_HPP
#define FUB_AMREX_DEBUG_OUTPUT_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

namespace fub::amrex {

class DebugStorage {
public:
  using Hierarchy = std::vector<::amrex::MultiFab>;
  using ComponentNames = std::vector<std::string>;

  /// \brief Initializes an empty storage.
  DebugStorage() = default;

  /// @{
  /// \brief Saves a current hierarchy state with given component names.
  ///
  /// The actual output will be handled by the DebugOutput class and will
  /// usually happen at a later time point.
  void SaveData(const ::amrex::MultiFab& mf, const std::string& name,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::MultiFab& mf, const ComponentNames& names,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void SaveData(const std::vector<const ::amrex::MultiFab*>& hierarchy,
                const std::string& name,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const std::vector<const ::amrex::MultiFab*>& hierarchy,
                const ComponentNames& names,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));
  /// @}

  /// \brief Deletes all currently stored data.
  void ClearAll();

  /// \brief Collects all hierachies and associated names which are stored on the
  /// specified location.
  ///
  /// This function will join all compatible hierachies to a common big one if
  /// possible. If two component names collide they will be renamed by numbering.
  std::vector<std::pair<Hierarchy, ComponentNames>>
  GatherFields(::amrex::IndexType location) const;

private:
  /// \brief Each SaveData will append a hierarchy
  std::vector<Hierarchy> saved_hierarchies_{};

  /// \brief Each SaveData will append a list of names on the hierachy
  std::vector<ComponentNames> names_per_hierarchy_{};
};

class DebugOutput : public OutputAtFrequencyOrInterval<GriddingAlgorithm> {
public:
  DebugOutput(const ProgramOptions& opts,
              std::shared_ptr<DebugStorage> storage);

  const std::shared_ptr<DebugStorage>& GetStorage() noexcept {
    return storage_;
  }

  void operator()(const GriddingAlgorithm& grid) override;

private:
  std::shared_ptr<DebugStorage> storage_{};
  std::string directory_;
};

} // namespace fub::amrex

#endif