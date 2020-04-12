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

#ifndef FUB_AMREX_INITIAL_DATA_HPP
#define FUB_AMREX_INITIAL_DATA_HPP

#include "fub/Meta.hpp"
#include "fub/core/type_traits.hpp"
#include "fub/Duration.hpp"

#include <memory>

namespace fub {

/// \defgroup InitialData Initial Data Methods
/// \brief This group contains all components which initialize simulation data.
/// \ingroup GriddingAlgorithm

namespace detail {
template <typename GriddingAlgorithm> struct InitialDataStrategy {
  using PatchLevel = typename GridTraits<GriddingAlgorithm>::PatchLevel;

  virtual ~InitialDataStrategy() = default;

  virtual void InitializeData(PatchLevel& patch_level,
                              const GriddingAlgorithm& grid, int level,
                              Duration time) = 0;

  virtual std::unique_ptr<InitialDataStrategy> Clone() const = 0;
};
} // namespace detail

/// \ingroup InitialData PolymorphicValueType
/// \brief This class is a polymoprhic value type which stores components to
/// initialize a gridding algorithm during its initialization procedure.
template <typename GriddingAlgorithm> class AnyInitialData {
public:
  using PatchLevel = typename GridTraits<GriddingAlgorithm>::PatchLevel;

  /// @{
  /// \name Constructors

  /// \brief Constructs an empty object that does no initialization.
  AnyInitialData() = default;

  /// \brief Stores and wraps the `initial_data` object
  template <typename T, typename = std::enable_if_t<
                            !decays_to<T, AnyInitialData<GriddingAlgorithm>>()>>
  AnyInitialData(T&& initial_data);

  /// \brief Copies the `other` implementation and invokes a memory allocation.
  AnyInitialData(const AnyInitialData& other);

  /// \brief Copies the `other` implementation and invokes a memory allocation.
  AnyInitialData& operator=(const AnyInitialData& other);

  /// \brief Moves the `other` implementation without allocation and leaves an
  /// empty object.
  AnyInitialData(AnyInitialData&&) noexcept = default;

  /// \brief Moves the `other` implementation without allocation and leaves an
  /// empty object.
  AnyInitialData& operator=(AnyInitialData&&) noexcept = default;
  /// @}

  /// @{
  /// \name Actions

  /// \brief Initializes a patch level within a gridding algorithm
  ///
  /// \note In cases where a gridding algorithm needs to rebuild a level, the
  /// InitializeData function might be invoked multiple times for the level.
  void InitializeData(PatchLevel& patch_level, const GriddingAlgorithm& grid,
                      int level, Duration time);
  /// @}

private:
  std::unique_ptr<detail::InitialDataStrategy<GriddingAlgorithm>> initial_data_;
};

///////////////////////////////////////////////////////////////////////////////
//                                                              Implementation

///////////////////////////////////////////////////////////////////////////////
//                                                                Constructors

namespace detail {
template <typename T, typename GriddingAlgorithm>
struct InitialDataWrapper : public InitialDataStrategy<GriddingAlgorithm> {
  using PatchLevel = typename GridTraits<GriddingAlgorithm>::PatchLevel;

  InitialDataWrapper(const T& initial_data) : initial_data_{initial_data} {}
  InitialDataWrapper(T&& initial_data)
      : initial_data_{std::move(initial_data)} {}

  void InitializeData(PatchLevel& patch_level, const GriddingAlgorithm& grid,
                      int level, Duration time) override {
    initial_data_.InitializeData(patch_level, grid, level, time);
  }

  std::unique_ptr<InitialDataStrategy<GriddingAlgorithm>>
  Clone() const override {
    return std::make_unique<InitialDataWrapper<T, GriddingAlgorithm>>(
        initial_data_);
  }

  T initial_data_;
};
} // namespace detail

template <typename GriddingAlgorithm>
template <typename T, typename>
AnyInitialData<GriddingAlgorithm>::AnyInitialData(T&& initial_data)
    : initial_data_{
          std::make_unique<detail::InitialDataWrapper<std::decay_t<T>, GriddingAlgorithm>>(
              std::forward<T>(initial_data))} {}

template <typename GriddingAlgorithm>
AnyInitialData<GriddingAlgorithm>::AnyInitialData(const AnyInitialData& other)
    : initial_data_(other.initial_data_->Clone()) {}

template <typename GriddingAlgorithm>
AnyInitialData<GriddingAlgorithm>&
AnyInitialData<GriddingAlgorithm>::operator=(const AnyInitialData& other) {
  AnyInitialData tmp(other);
  return *this = std::move(tmp);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                     Actions

template <typename GriddingAlgorithm>
void AnyInitialData<GriddingAlgorithm>::InitializeData(
    PatchLevel& patch_level, const GriddingAlgorithm& grid, int level,
    Duration time) {
  if (initial_data_) {
    return initial_data_->InitializeData(patch_level, grid, level, time);
  }
}

} // namespace fub

#endif
