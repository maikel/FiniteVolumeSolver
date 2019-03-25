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

#ifndef FUB_GRID_AMREX_CUTCELL_TAGGING_HPP
#define FUB_GRID_AMREX_CUTCELL_TAGGING_HPP

#include "fub/CartesianCoordinates.hpp"
#include "fub/core/mdspan.hpp"
#include "fub/grid/AMReX/PatchHandle.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"
#include "fub/grid/AMReX/cutcell/PatchHierarchy.hpp"

#include <AMReX.H>

namespace fub {
namespace amrex {
namespace cutcell {

template <typename Equation, typename... Tagging> struct AdaptTagging {
  AdaptTagging(Equation equation,
               std::shared_ptr<const PatchHierarchy> hierarchy,
               Tagging... tagging)
      : tagging_{std::move(tagging)...}, equation_{std::move(equation)},
        hierarchy_{std::move(hierarchy)} {}

  void TagCellsForRefinement(
      const PatchDataView<char, AMREX_SPACEDIM>& tags,
      const PatchDataView<const double, AMREX_SPACEDIM + 1>& states,
      const PatchHandle& patch) {
    View<const Complete<Equation>> state_view =
        MakeView<View<Complete<Equation>>>(states, equation_);
    const ::amrex::Geometry& geom = hierarchy_->GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->growntilebox();
    const ::amrex::EBFArrayBoxFactory& factory =
        hierarchy_->GetEmbeddedBoundary(patch.level);
    const auto& flags = factory.getMultiEBCellFlagFab()[*patch.iterator];
    const ::amrex::FabType type = flags.getType(box);
    if (type == ::amrex::FabType::regular) {
      boost::hana::for_each(tagging_, [&](auto&& tagging) {
        tagging.TagCellsForRefinement(tags, state_view,
                                      GetCartesianCoordinates(geom, box));
      });
    } else {
      FUB_ASSERT(type == ::amrex::FabType::singlevalued);
      auto cutcell_data = hierarchy_->GetCutCellData(patch, Direction::X);
      boost::hana::for_each(tagging_, [&](auto&& tagging) {
        tagging.TagCellsForRefinement(tags, state_view, cutcell_data,
                                      GetCartesianCoordinates(geom, box));
      });
    }
  }

  boost::hana::tuple<Tagging...> tagging_;
  Equation equation_;
  std::shared_ptr<const PatchHierarchy> hierarchy_;
};

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
