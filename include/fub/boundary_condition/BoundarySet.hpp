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

#ifndef FUB_BOUNDARY_CONDITION_BOUNDARIES_HPP
#define FUB_BOUNDARY_CONDITION_BOUNDARIES_HPP

#include "fub/AMReX/BoundaryCondition.hpp"
#include "fub/Direction.hpp"
#include "fub/PatchDataView.hpp"

#include <functional>
#include <vector>

namespace fub {

template <typename Hierarchy> class BoundarySet {
public:
  using PatchHandle = typename Hierarchy::PatchHandle;
  using BoundaryCondition = std::function<void(
      const PatchDataView<double, AMREX_SPACEDIM + 1>&, const Hierarchy&,
      typename Hierarchy::PatchHandle, Location, int, Duration)>;

  BoundarySet() : conditions_(2 * AMREX_SPACEDIM) {}

  void operator()(const PatchDataView<double, AMREX_SPACEDIM + 1>& data,
                  const Hierarchy& hierarchy, PatchHandle patch, Location loc,
                  int fillwidth, Duration timepoint) const {
    const std::size_t i = static_cast<std::size_t>(
        loc.direction * 2 + static_cast<std::size_t>(loc.side));
    FUB_ASSERT(i < conditions_.size());
    if (conditions_[i]) {
      conditions_[i](data, hierarchy, patch, loc, fillwidth, timepoint);
    }
  }

  void SetBoundaryCondition(Location loc, BoundaryCondition condition) {
    const std::size_t i = static_cast<std::size_t>(
        loc.direction * 2 + static_cast<std::size_t>(loc.side));
    FUB_ASSERT(i < conditions_.size());
    conditions_[i] = std::move(condition);
  }

private:
  std::vector<BoundaryCondition> conditions_;
};

} // namespace fub

#endif
