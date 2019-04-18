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

#include "fub/CartesianCoordinates.hpp"
#include "fub/Equation.hpp"
#include "fub/PatchDataView.hpp"
#include "fub/grid/AMReX/PatchHandle.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"

#include <AMReX.H>

namespace fub {
namespace amrex {
struct InitialDataStrategy {
  virtual ~InitialDataStrategy() = default;
  virtual void
  InitializeData(const PatchDataView<double, AMREX_SPACEDIM + 1>& states,
                 const PatchHandle& patch) = 0;
};

template <typename T> struct InitialDataWrapper : public InitialDataStrategy {
  InitialDataWrapper(const T& initial_data) : initial_data_{initial_data} {}
  InitialDataWrapper(T&& initial_data)
      : initial_data_{std::move(initial_data)} {}

  void InitializeData(const PatchDataView<double, AMREX_SPACEDIM + 1>& states,
                      const PatchHandle& patch) override {
    initial_data_.InitializeData(states, patch);
  }

  T initial_data_;
};

struct InitialData {
  template <typename T>
  InitialData(const T& initial_data)
      : initial_data_{std::make_unique<InitialDataWrapper<T>>(initial_data)} {}
  template <typename T>
  InitialData(T&& initial_data)
      : initial_data_{std::make_unique<InitialDataWrapper<remove_cvref_t<T>>>(
            std::move(initial_data))} {}

  void InitializeData(const PatchDataView<double, AMREX_SPACEDIM + 1>& states,
                      const PatchHandle& patch) {
    if (initial_data_) {
      return initial_data_->InitializeData(states, patch);
    }
  }

  std::unique_ptr<InitialDataStrategy> initial_data_;
};

template <typename InitialData, typename Equation> struct AdaptInitialData {
  AdaptInitialData(InitialData data, Equation equation)
      : data_{std::move(data)}, equation_{std::move(equation)} {}

  void InitializeData(const PatchDataView<double, AMREX_SPACEDIM + 1>& states,
                      const PatchHandle& patch) {
    const IndexBox<AMREX_SPACEDIM> cells =
        AsIndexBox(patch.iterator->tilebox());
    View<Complete<Equation>> state_view = Subview(
        MakeView<BasicView<Complete<Equation>>>(states, equation_), cells);
    data_.InitializeData(state_view, patch);
  }

  InitialData data_;
  Equation equation_;
};

} // namespace amrex
} // namespace fub

#endif
