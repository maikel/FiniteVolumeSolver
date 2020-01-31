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
void DebugStorage::SaveData(const ::amrex::MultiFab& mf,
                            const std::string& name,
                            ::amrex::SrcComp component) {
  Hierarchy single_level_hierarchy{};
  ::amrex::BoxArray ba = mf.boxArray();
  ::amrex::DistributionMapping dm = mf.DistributionMap();
  ::amrex::MultiFab& dest =
      single_level_hierarchy.emplace_back(ba, dm, component.i, 0);
  dest.copy(mf);
  saved_hierarchies_.push_back(std::move(single_level_hierarchy));
  names_per_hierarchy_.push_back(std::vector<std::string>{name});
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
    }
    int counter = 1;
    do {
      *pos = fmt::format("{}_{}", *pos, counter);
      counter += 1;
      pos = std::find(first, last, candidate);
    } while (last != pos);
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
  std::vector<std::vector<int>> partitions = Partitions_(box_arrays);
  std::vector<
      std::pair<std::vector<::amrex::MultiFab>, std::vector<std::string>>>
      hierarchies{};
  for (const std::vector<int>& partition : partitions) {
    std::vector<std::string> all_components =
        Select_(names_per_hierarchy_, partition);
    int n_components = static_cast<int>(all_components.size());
    std::vector<::amrex::MultiFab> hierarchy{};
    std::vector<BoxArray> box_array = box_arrays[partition[0]];
    std::vector<DistributionMapping> distribution =
        distribution_maps[partition[0]];
    for (std::size_t level = 0; level < box_array.size(); ++level) {
      MultiFab& dest = hierarchy.emplace_back(
          box_array[level], distribution[level], n_components, 0);
      int comp = 0;
      for (int i : partition) {
        const MultiFab& src = saved_hierarchies_[i][level];
        const int n = src.nComp();
        MultiFab::Copy(dest, src, SrcComp_(0), DestComp_(comp), NumComps_(n),
                       0);
        comp += n;
      }
      FUB_ASSERT(comp == n_components);
    }
    hierarchies.emplace_back(std::move(hierarchy), std::move(all_components));
  }
  return hierarchies;
}

void DebugStorage::ClearAll() {
  names_per_hierarchy_.clear();
  saved_hierarchies_.clear();
}

} // namespace fub::amrex