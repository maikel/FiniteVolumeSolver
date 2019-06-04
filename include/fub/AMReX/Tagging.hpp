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

#include <AMReX.H>
#include <AMReX_TagBox.H>

namespace fub::amrex {

class GriddingAlgorithm;

struct TaggingStrategy {
  virtual ~TaggingStrategy() = default;
  virtual std::unique_ptr<TaggingStrategy> Clone() const = 0;
  virtual void TagCellsForRefinement(::amrex::TagBoxArray& tags,
                                     Duration time_point, int level,
                                     GriddingAlgorithm& gridding) = 0;
};

template <typename T> struct TaggingWrapper : public TaggingStrategy {
  TaggingWrapper(const T& tag) : tag_{tag} {}
  TaggingWrapper(T&& tag) : tag_{std::move(tag)} {}

  void TagCellsForRefinement(::amrex::TagBoxArray& tags, Duration time_point,
                             int level, GriddingAlgorithm& gridding) override {
    tag_.TagCellsForRefinement(tags, time_point, level, gridding);
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
  Tagging(T&& tag)
      : tag_{std::make_unique<TaggingWrapper<remove_cvref_t<T>>>(
            std::forward<T>(tag))} {}

  void TagCellsForRefinement(::amrex::TagBoxArray& tags, Duration time_point,
                             int level, GriddingAlgorithm& gridding) {
    return tag_->TagCellsForRefinement(tags, time_point, level, gridding);
  }

  std::unique_ptr<TaggingStrategy> tag_;
};

template <typename... Tagging> struct TagAllOf {
  TagAllOf(Tagging... tagging) : tagging_{std::move(tagging)...} {}

  void TagCellsForRefinement(::amrex::TagBoxArray& tags, Duration time_point,
                             int level, GriddingAlgorithm& gridding) {
    std::apply(
        [this, &tags, time_point, level, &gridding](Tagging&... tagging) {
          (tagging.TagCellsForRefinement(tags, time_point, level, gridding),
           ...);
        },
        tagging_);
  }

  std::tuple<Tagging...> tagging_;
};

} // namespace fub::amrex

#endif
