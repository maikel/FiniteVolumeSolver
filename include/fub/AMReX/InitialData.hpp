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

#include "fub/core/type_traits.hpp"
#include "fub/AMReX/PatchHierarchy.hpp"

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include <memory>

namespace fub {
namespace amrex {

template <typename GriddingAlgorithm>
struct InitialDataStrategy {
  virtual ~InitialDataStrategy() = default;

  virtual void InitializeData(PatchLevel& patch_level,
                              const GriddingAlgorithm& grid, int level, Duration time) = 0;

  virtual std::unique_ptr<InitialDataStrategy> Clone() const = 0;
};

template <typename T, typename GriddingAlgorithm>
struct InitialDataWrapper : public InitialDataStrategy<GriddingAlgorithm> {
  InitialDataWrapper(const T& initial_data) : initial_data_{initial_data} {}
  InitialDataWrapper(T&& initial_data)
      : initial_data_{std::move(initial_data)} {}

  void InitializeData(PatchLevel& patch_level,
                      const GriddingAlgorithm& grid, int level, Duration time) override {
    initial_data_.InitializeData(patch_level, grid, level, time);
  }

  std::unique_ptr<InitialDataStrategy<GriddingAlgorithm>> Clone() const override {
    return std::make_unique<InitialDataWrapper<T, GriddingAlgorithm>>(initial_data_);
  }

  T initial_data_;
};

template <typename GriddingAlgorithm>
struct AnyInitialData {
  AnyInitialData() = default;

  AnyInitialData(const AnyInitialData& other)
      : initial_data_(other.initial_data_->Clone()) {}
  AnyInitialData& operator=(const AnyInitialData& other) {
    AnyInitialData tmp(other);
    return *this = std::move(tmp);
  }

  AnyInitialData(AnyInitialData&&) noexcept = default;
  AnyInitialData& operator=(AnyInitialData&&) noexcept = default;

  template <typename T>
  AnyInitialData(const T& initial_data)
      : initial_data_{std::make_unique<InitialDataWrapper<T, GriddingAlgorithm>>(initial_data)} {}

  template <typename T>
  AnyInitialData(T&& initial_data)
      : initial_data_{std::make_unique<InitialDataWrapper<remove_cvref_t<T>, GriddingAlgorithm>>(
            std::move(initial_data))} {}

  void InitializeData(PatchLevel& patch_level, const GriddingAlgorithm& grid,
                      int level, Duration time) {
    if (initial_data_) {
      return initial_data_->InitializeData(patch_level, grid, level, time);
    }
  }

  std::unique_ptr<InitialDataStrategy<GriddingAlgorithm>> initial_data_;
};

} // namespace amrex
} // namespace fub

#endif
