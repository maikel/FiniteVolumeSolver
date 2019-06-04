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

#ifndef FUB_SAMRAI_TAGGING_HPP
#define FUB_SAMRAI_TAGGING_HPP

#include "fub/CartesianCoordinates.hpp"
#include "fub/SAMRAI/ViewPatch.hpp"
#include "fub/core/mdspan.hpp"

namespace fub {
namespace samrai {
struct TaggingStrategy {
  virtual ~TaggingStrategy() = default;
  virtual void
  TagCellsForRefinement(SAMRAI::pdat::CellData<int>& tags,
                        span<SAMRAI::pdat::CellData<double>*> states,
                        const CartesianCoordinates& coords) = 0;

  virtual std::unique_ptr<TaggingStrategy> Clone() const = 0;
};

template <typename T> struct TaggingWrapper : public TaggingStrategy {
  TaggingWrapper(const T& tag) : tag_{tag} {}
  TaggingWrapper(T&& tag) : tag_{std::move(tag)} {}

  void TagCellsForRefinement(SAMRAI::pdat::CellData<int>& tags,
                             span<SAMRAI::pdat::CellData<double>*> states,
                             const CartesianCoordinates& coords) override {
    tag_.TagCellsForRefinement(tags, states, coords);
  }

  std::unique_ptr<TaggingStrategy> Clone() const override {
    return std::make_unique<TaggingWrapper<T>>(tag_);
  }

  T tag_;
};

struct Tagging {
  Tagging(const Tagging& other) : tag_{other.tag_->Clone()} {}

  Tagging& operator=(const Tagging& other) {
    tag_ = other.tag_->Clone();
    return *this;
  }

  Tagging(Tagging&& other) = default;
  Tagging& operator=(Tagging&& other) = default;

  template <typename T>
  Tagging(const T& tag) : tag_{std::make_unique<TaggingWrapper<T>>(tag)} {}
  template <typename T>
  Tagging(T&& tag)
      : tag_{std::make_unique<TaggingWrapper<remove_cvref_t<T>>>(
            std::move(tag))} {}

  void TagCellsForRefinement(SAMRAI::pdat::CellData<int>& tags,
                             span<SAMRAI::pdat::CellData<double>*> states,
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

  void TagCellsForRefinement(SAMRAI::pdat::CellData<int>& tags,
                             span<SAMRAI::pdat::CellData<double>*> states,
                             const CartesianCoordinates& coords) {
    auto state_view = MakeView<Complete<Equation>>(states, equation_);
    auto tags_view = MakeMdSpan<Equation::Rank()>(tags.getArrayData());
    tagging_.TagCellsForRefinement(tags_view, state_view, coords);
  }

  Tagging tagging_;
  Equation equation_;
};

} // namespace samrai
} // namespace fub

#endif
