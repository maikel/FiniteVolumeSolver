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

#include "GriddingAlgorithm.hpp"
#include "fub/CartesianCoordinates.hpp"
#include "fub/SAMRAI/ViewPatch.hpp"
#include "fub/core/mdspan.hpp"

namespace fub {
namespace samrai {

struct GriddingAlgorithm;

struct TaggingStrategy {
  virtual ~TaggingStrategy() = default;
  virtual void TagCellsForRefinement(GriddingAlgorithm& gridding, int level, int tag_id, Duration time_point) = 0;

  virtual std::unique_ptr<TaggingStrategy> Clone() const = 0;
};

template <typename T> struct TaggingWrapper : public TaggingStrategy {
  TaggingWrapper(const T& tag) : tag_{tag} {}
  TaggingWrapper(T&& tag) : tag_{std::move(tag)} {}

  void TagCellsForRefinement(GriddingAlgorithm& gridding, int level, int tag_id, Duration time_point) override {
    tag_.TagCellsForRefinement(gridding, level, tag_id, time_point);
  }

  std::unique_ptr<TaggingStrategy> Clone() const override {
    return std::make_unique<TaggingWrapper<T>>(tag_);
  }

  T tag_;
};

struct Tagging {
  Tagging() = default;

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

  void TagCellsForRefinement(GriddingAlgorithm& gridding, int level, int tag_id, Duration time_point) {
    if (tag_) {
      return tag_->TagCellsForRefinement(gridding, level, tag_id, time_point);
    }
  }

  std::unique_ptr<TaggingStrategy> tag_;
};

} // namespace samrai
} // namespace fub

#endif
