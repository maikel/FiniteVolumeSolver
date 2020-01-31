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

namespace fub::amrex {
namespace {
void SaveCells_(std::vector<ComponentNames>& names,
                std::vector<Hierarchy>& cells, const ::amrex::MultiFab& src,
                const std::string& name, ::amrex::SrcComp component) {
  Hierarchy single_level_hierarchy{};
  ::amrex::BoxArray ba = cells.boxArray();
  ::amrex::DistributionMapping dm = cells.DistributionMap();
  ::amrex::MultiFab& dest = single_level_hierarchy.emplace_back(ba, dm, 1, 0);
  dest.copy(src);
  cells.push_back(std::move(single_level_hierarchy));
  names.push_back(ComponentNames{name});
}
} // namespace

void DebugOutput::SaveData(const ::amrex::MultiFab& mf, const std::string& name,
                           ::amrex::SrcComp component = ::amrex::SrcComp(0)) {
  if (mf.boxArray.boxType() == ::amrex::IntVect{}) {
    SaveCells_(names_, cells_, mf, name, component);
  } else {
    throw std::runtime_error{
        "Debugging face or node cenetered data is not supported yet."};
  }
}

namespace {
struct SrcComp_ {
  explicit SrcComp_(int ii) : i{ii} {}
  int i;
  operator int() const noexcept { return i; }
};

struct DestComp_ {
  explicit DestComp_(int ii) : i{ii} {}
  int i;
  operator int() const noexcept { return i; }
};

struct NumComps_ {
  explicit NumComps_(int ii) : i{ii} {}
  int i;
  operator int() const noexcept { return i; }
};

struct equal_to {
  int i;
  const std::vector<std::vector<::amrex::BoxArray>>* bas_;
  bool operator()(int j) { return bas_[i] == bas_[j]; }
};

std::vector<std::vector<int>>
Partitions_(const std::vector<std::vector<::amrex::BoxArray>>& boxes) {
  std::vector<int> selection(boxes.size());
  std::vector<std::vector<int>> partitions{};
  auto first = selection.begin();
  auto last = selection.end();
  std::iota(first, last);
  auto next = std::partition(first, last, equal_to(*first, &boxes));
  partitions.emplace_back(first, next);
  while (next != last) {
    first = next;
    next = std::partition(first, last, equal_to(*first, &boxes));
    partitions.emplace_back(first, next);
  }
  return partitions;
}

} // namespace

std::vector<std::pair<Hierarchy, ComponentNames>>
DebugOutput::GatherCellCentered() const {
  using namespace ::amrex;
  std::vector<std::vector<BoxArray>> box_arrays = GatherBoxArrays_(cells_);
  std::vector<std::vector<DistributionMapping>> distribution_maps =
      GatherDistributionMappings_(cells_);
  std::vector<std::vector<int>> partitions = Partitions_(box_arrays);
  std::vector<std::pair<Hierarchy, ComponentNames>> hierarchies{};
  for (const std::vector<int>& partition : partitions) {
    ComponentNames all_components = Select(names_, partition);
    int n_components = static_cast<int>(all_components.size());
    Hierarchy hierarchy{};
    std::vector<BoxArray> box_array = box_arrays[partition[0]];
    std::vector<DistributionMapping> distribution =
        distribution_maps[partition[0]];
    for (std::size_t level = 0; level < box_array.size(); ++level) {
      MultiFab& dest = hierarchy.emplace_back(
          box_array[level], distribution[level], n_components, no_ghosts);
      int comp = 0;
      for (int i : partition) {
        const MultiFab& src = cells_[i][level];
        const int n = src.nComp();
        MultiFab::Copy(dest, src, SrcComp_(0), DestComp_(comp), NumComps_(n));
        comp += n;
      }
      FUB_ASSERT(comp == n_components);
    }
    hierarchies.emplace_back(std::move(hierarchy), std::move(all_components));
  }
  return hierarchies;
}

void DebugOutput::ClearAll() {
  names_.clear();
  cells_.clear();
}

} // namespace fub::amrex