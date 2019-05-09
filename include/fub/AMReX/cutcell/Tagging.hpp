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

#include "fub/AMReX/PatchHandle.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/cutcell/PatchHierarchy.hpp"
#include "fub/CartesianCoordinates.hpp"
#include "fub/core/mdspan.hpp"

#include <AMReX.H>

namespace fub {
namespace amrex {
namespace cutcell {

struct TaggingStrategy {
  virtual ~TaggingStrategy() = default;
  virtual std::unique_ptr<TaggingStrategy> Clone() const = 0;
  virtual void TagCellsForRefinement(
      const PatchDataView<char, AMREX_SPACEDIM>& tags,
      const PatchDataView<const double, AMREX_SPACEDIM + 1>& states,
      const PatchHierarchy& hierarchy, const PatchHandle& patch) = 0;
};

template <typename T> struct TaggingWrapper : public TaggingStrategy {
  TaggingWrapper(const T& tag) : tag_{tag} {}
  TaggingWrapper(T&& tag) : tag_{std::move(tag)} {}

  void TagCellsForRefinement(
      const PatchDataView<char, AMREX_SPACEDIM>& tags,
      const PatchDataView<const double, AMREX_SPACEDIM + 1>& states,
      const PatchHierarchy& hierarchy, const PatchHandle& patch) override {
    tag_.TagCellsForRefinement(tags, states, hierarchy, patch);
  }

  std::unique_ptr<TaggingStrategy> Clone() const override {
    return std::make_unique<TaggingWrapper<T>>(tag_);
  };

  T tag_;
};

struct Tagging {
  Tagging() = default;

  Tagging(const Tagging& other) : tag_(other.tag_->Clone()) {}
  Tagging(Tagging&&) noexcept = default;

  Tagging& operator=(const Tagging& other) {
    Tagging tmp(other);
    return *this = std::move(tmp);
  }
  Tagging& operator=(Tagging&&) noexcept = default;

  template <typename T>
  Tagging(const T& tag) : tag_{std::make_unique<TaggingWrapper<T>>(tag)} {}
  template <typename T>
  Tagging(T&& tag)
      : tag_{std::make_unique<TaggingWrapper<remove_cvref_t<T>>>(
            std::move(tag))} {}

  void TagCellsForRefinement(
      const PatchDataView<char, AMREX_SPACEDIM>& tags,
      const PatchDataView<const double, AMREX_SPACEDIM + 1>& states,
      const PatchHierarchy& hierarchy, const PatchHandle& patch) {
    if (tag_) {
      return tag_->TagCellsForRefinement(tags, states, hierarchy, patch);
    }
  }

  std::unique_ptr<TaggingStrategy> tag_;
};

template <typename Equation, typename... Tagging> struct AdaptTagging {
  AdaptTagging(Equation equation, Tagging... tagging)
      : tagging_{std::move(tagging)...}, equation_{std::move(equation)} {}

  void TagCellsForRefinement(
      const PatchDataView<char, AMREX_SPACEDIM>& tags,
      const PatchDataView<const double, AMREX_SPACEDIM + 1>& states,
      const PatchHierarchy& hierarchy, const PatchHandle& patch) {
    BasicView<const Complete<Equation>> state_view =
        MakeView<BasicView<const Complete<Equation>>>(states, equation_);
    const ::amrex::Geometry& geom = hierarchy.GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->growntilebox();
    const std::shared_ptr<::amrex::EBFArrayBoxFactory>& factory =
        hierarchy.GetEmbeddedBoundary(patch.level);
    const auto& flags = factory->getMultiEBCellFlagFab()[*patch.iterator];
    const ::amrex::FabType type = flags.getType(box);
    if (type == ::amrex::FabType::regular) {
      boost::mp11::tuple_for_each(tagging_, [&](auto&& tagging) {
        tagging.TagCellsForRefinement(tags, state_view,
                                      GetCartesianCoordinates(geom, box));
      });
    } else if (type != ::amrex::FabType::covered) {
      FUB_ASSERT(type == ::amrex::FabType::singlevalued);
      auto cutcell_data = hierarchy.GetCutCellData(patch, Direction::X);
      boost::mp11::tuple_for_each(tagging_, [&](auto&& tagging) {
        tagging.TagCellsForRefinement(tags, state_view, cutcell_data,
                                      GetCartesianCoordinates(geom, box));
      });
    }
  }

  std::tuple<Tagging...> tagging_;
  Equation equation_;
};

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
