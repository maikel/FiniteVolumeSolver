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

#ifndef FUB_AMREX_CUTCELL_INITIAL_DATA_HPP
#define FUB_AMREX_CUTCELL_INITIAL_DATA_HPP

#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/cutcell/PatchHierarchy.hpp"
#include "fub/CartesianCoordinates.hpp"
#include "fub/Equation.hpp"
#include "fub/PatchDataView.hpp"

#include <AMReX.H>

namespace fub {
namespace amrex {
namespace cutcell {

struct InitialDataStrategy {
  virtual ~InitialDataStrategy() = default;

  virtual void InitializeData(PatchLevel& patch_level,
                              const fub::amrex::cutcell::GriddingAlgorithm& grid, int level, Duration time) = 0;

  virtual std::unique_ptr<InitialDataStrategy> Clone() const = 0;
};

template <typename T> struct InitialDataWrapper : public InitialDataStrategy {
  InitialDataWrapper(const T& initial_data) : initial_data_{initial_data} {}
  InitialDataWrapper(T&& initial_data)
      : initial_data_{std::move(initial_data)} {}

  void InitializeData(PatchLevel& patch_level,
                      const fub::amrex::cutcell::GriddingAlgorithm& grid, int level, Duration time) override {
    initial_data_.InitializeData(patch_level, grid, level, time);
  }

  std::unique_ptr<InitialDataStrategy> Clone() const override {
    return std::make_unique<InitialDataWrapper<T>>(initial_data_);
  }

  T initial_data_;
};

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
      : initial_data_{std::make_unique<InitialDataWrapper<T>>(initial_data)} {}

  template <typename T>
  AnyInitialData(T&& initial_data)
      : initial_data_{std::make_unique<InitialDataWrapper<remove_cvref_t<T>>>(
            std::move(initial_data))} {}

  void InitializeData(PatchLevel& patch_level, const fub::amrex::cutcell::GriddingAlgorithm& grid, int level, Duration time) {
    if (initial_data_) {
      return initial_data_->InitializeData(patch_level, grid, level, time);
    }
  }

  std::unique_ptr<InitialDataStrategy> initial_data_;
};

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
