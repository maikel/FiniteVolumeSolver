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

#ifndef FUB_AMREX_TAGGING_HPP
#define FUB_AMREX_TAGGING_HPP

#include "fub/CartesianCoordinates.hpp"
#include "fub/core/mdspan.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"

#include <AMReX.H>

namespace fub {
namespace amrex {
struct TaggingStrategy {
  virtual ~TaggingStrategy() = default;
  virtual void
  TagCellsForRefinement(mdspan<char, AMREX_SPACEDIM, layout_stride> tags,
                        mdspan<const double, AMREX_SPACEDIM + 1> states,
                        const CartesianCoordinates& coords) = 0;
};

template <typename T> struct TaggingWrapper : public TaggingStrategy {
  TaggingWrapper(const T& tag) : tag_{tag} {}
  TaggingWrapper(T&& tag) : tag_{std::move(tag)} {}

  void
  TagCellsForRefinement(mdspan<char, AMREX_SPACEDIM, layout_stride> tags,
                        mdspan<const double, AMREX_SPACEDIM + 1> states,
                        const CartesianCoordinates& coords) override {
    tag_.TagCellsForRefinement(tags, states, coords);
  }

  T tag_;
};

struct Tagging {
  template <typename T>
  Tagging(const T& tag) : tag_{std::make_unique<TaggingWrapper<T>>(tag)} {}
  template <typename T>
  Tagging(T&& tag)
      : tag_{std::make_unique<TaggingWrapper<remove_cvref_t<T>>>(
            std::move(tag))} {}

  void
  TagCellsForRefinement(mdspan<char, AMREX_SPACEDIM, layout_stride> tags,
                        mdspan<const double, AMREX_SPACEDIM + 1> states,
                        const CartesianCoordinates& coords) {
    if (tag_) {
      return tag_->TagCellsForRefinement(tags, states, coords);
    }
  }

  std::unique_ptr<TaggingStrategy> tag_;
};

template <typename Tagging, typename Equation> struct AdaptTagging {
  AdaptTagging(Tagging tagging, Equation equation)
      : tagging_{std::move(tagging)}, equation_{std::move(equation)} {}

  void
  TagCellsForRefinement(mdspan<char, AMREX_SPACEDIM, layout_stride> tags,
                        mdspan<const double, AMREX_SPACEDIM + 1> states,
                        const CartesianCoordinates& coords) {
    auto state_view = MakeView(boost::hana::type_c<View<Complete<Equation>>>,
                           states, equation_);
    tagging_.TagCellsForRefinement(tags, state_view, coords);
  }

  Tagging tagging_;
  Equation equation_;
};

} // namespace amrex
} // namespace fub

#endif