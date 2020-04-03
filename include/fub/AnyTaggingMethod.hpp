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

#ifndef FUB_TAGGING_METHOD_HPP
#define FUB_TAGGING_METHOD_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Meta.hpp"

#include <memory>

namespace fub {

/// \defgroup TaggingMethod Tagging Methods
/// \ingroup GriddingAlgorithm
/// This modules collects all components that tag cells for refinement.

template <typename GriddingAlgorithm> struct TaggingMethodStrategy_ {
  using TagDataHandle = typename GridTraits<GriddingAlgorithm>::TagDataHandle;
  virtual ~TaggingMethodStrategy_() = default;
  virtual std::unique_ptr<TaggingMethodStrategy_> Clone() const = 0;
  virtual void TagCellsForRefinement(TagDataHandle tags,
                                     GriddingAlgorithm& gridding, int level,
                                     Duration time_point) = 0;
};

/// \ingroup TaggingMethod PolymorphicValueType
/// \brief This class is a polymorphic value type that stores objects which
/// satisfies the TaggingMethod<GriddingAlgorithm> concept.
template <typename GriddingAlgorithm> class AnyTaggingMethod {
public:
  using TagDataHandle = typename GridTraits<GriddingAlgorithm>::TagDataHandle;

  /// @{
  /// \name Constructors

  /// \brief This constructs a method that does nothing on invocation.
  AnyTaggingMethod() = default;

  /// \brief Stores any object which satisfies the tagging method concept.
  template <typename T,
            typename = std::enable_if_t<!decays_to<T, AnyTaggingMethod>()>>
  AnyTaggingMethod(T&& tag);

  /// \brief Copies the implementation.
  AnyTaggingMethod(const AnyTaggingMethod& other);

  /// \brief Copies the implementation.
  AnyTaggingMethod& operator=(const AnyTaggingMethod& other);

  /// \brief Moves the `other` object without allocating and leaves an empty
  /// method.
  AnyTaggingMethod(AnyTaggingMethod&& other) noexcept = default;

  /// \brief Moves the `other` object without allocating and leaves an empty
  /// method.
  AnyTaggingMethod& operator=(AnyTaggingMethod&& other) noexcept = default;
  /// @}

  /// @{
  /// \name Member functions

  /// \brief Mask cells that need further refinement in a regridding procedure.
  void TagCellsForRefinement(TagDataHandle tags, GriddingAlgorithm& gridding,
                             int level, Duration t);
  /// @}

private:
  std::unique_ptr<TaggingMethodStrategy_<GriddingAlgorithm>> tag_;
};

///////////////////////////////////////////////////////////////////////////////
//                                                              Implementation

template <typename GriddingAlgorithm, typename T>
struct TaggingMethodWrapper_
    : public TaggingMethodStrategy_<GriddingAlgorithm> {
  using TagDataHandle = typename GridTraits<GriddingAlgorithm>::TagDataHandle;

  TaggingMethodWrapper_(const T& tag) : tag_{tag} {}
  TaggingMethodWrapper_(T&& tag) : tag_{std::move(tag)} {}

  void TagCellsForRefinement(TagDataHandle tags, GriddingAlgorithm& gridding,
                             int level, Duration time_point) override {
    tag_.TagCellsForRefinement(tags, gridding, level, time_point);
  }

  std::unique_ptr<TaggingMethodStrategy_<GriddingAlgorithm>>
  Clone() const override {
    return std::make_unique<TaggingMethodWrapper_<GriddingAlgorithm, T>>(tag_);
  };

  T tag_;
};

template <typename GriddingAlgorithm>
AnyTaggingMethod<GriddingAlgorithm>::AnyTaggingMethod(
    const AnyTaggingMethod& other)
    : tag_(other.tag_ ? other.tag_->Clone() : nullptr) {}

template <typename GriddingAlgorithm>
AnyTaggingMethod<GriddingAlgorithm>&
AnyTaggingMethod<GriddingAlgorithm>::operator=(const AnyTaggingMethod& other) {
  AnyTaggingMethod tmp(other);
  return *this = std::move(tmp);
}

template <typename GriddingAlgorithm>
template <typename T, typename>
AnyTaggingMethod<GriddingAlgorithm>::AnyTaggingMethod(T&& tag)
    : tag_{std::make_unique<
          TaggingMethodWrapper_<GriddingAlgorithm, remove_cvref_t<T>>>(
          std::forward<T>(tag))} {}

template <typename GriddingAlgorithm>
void AnyTaggingMethod<GriddingAlgorithm>::TagCellsForRefinement(
    TagDataHandle tags, GriddingAlgorithm& gridding, int level,
    Duration time_point) {
  if (tag_) {
    return tag_->TagCellsForRefinement(tags, gridding, level, time_point);
  }
}

} // namespace fub

#endif
