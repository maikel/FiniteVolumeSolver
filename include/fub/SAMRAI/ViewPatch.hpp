// Copyright (c) 2018 Maikel Nadolski
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

#ifndef FUB_SAMRAI_PATCH_VIEW_HPP
#define FUB_SAMRAI_PATCH_VIEW_HPP

#include "fub/Equation.hpp"

#include <SAMRAI/pdat/ArrayData.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/SideData.h>

namespace fub {
namespace samrai {

template <int Rank, typename T>
mdspan<T, Rank> MakeMdSpan(SAMRAI::pdat::ArrayData<T>& array) {
  const SAMRAI::hier::Box& box = array.getBox();
  const int dim = box.getDim().getValue();
  const int depth = array.getDepth();
  if ((depth > 1 && Rank < dim + 1) || (depth == 1 && Rank < dim)) {
    throw std::invalid_argument("Dimension mismatch in make_mdspan.");
  }
  std::array<std::ptrdiff_t, Rank> extents;
  extents.fill(1);
  for (int d = 0; d < dim; ++d) {
    extents[d] = box.numberCells(d);
  }
  if (depth > 1) {
    extents[dim] = depth;
  }
  T* pointer = array.getPointer();
  return MdSpan(pointer, extents);
}

template <int Rank, typename T>
mdspan<const T, Rank> MakeMdSpan(const SAMRAI::pdat::ArrayData<T>& array) {
  const SAMRAI::hier::Box& box = array.getBox();
  const int dim = box.getDim().getValue();
  const int depth = array.getDepth();
  if ((depth > 1 && Rank < dim + 1) || (depth == 1 && Rank < dim)) {
    throw std::invalid_argument("Dimension mismatch in make_mdspan.");
  }
  std::array<std::ptrdiff_t, Rank> extents;
  extents.fill(1);
  for (int d = 0; d < dim; ++d) {
    extents[d] = box.numberCells(d);
  }
  if (depth > 1) {
    extents[dim] = depth;
  }
  const T* pointer = array.getPointer();
  return MdSpan(pointer, extents);
}

template <int Rank> IndexBox<Rank> AsIndexBox(const SAMRAI::hier::Box& box) {
  IndexBox<Rank> index_box{};
  const std::size_t dim = static_cast<std::size_t>(box.getDim().getValue());
  for (std::size_t i = 0; i < dim; ++i) {
    index_box.lower[i] = box.lower(static_cast<SAMRAI::hier::Box::dir_t>(i));
    index_box.upper[i] = box.upper(static_cast<SAMRAI::hier::Box::dir_t>(i));
  }
  return index_box;
}

template <int Rank, typename T>
PatchDataView<T, Rank> MakePatchDataView(SAMRAI::pdat::ArrayData<T>& array) {
  IndexBox<Rank> box = AsIndexBox<Rank>(array.getBox());
  return PatchDataView<T, Rank>(MakeMdSpan(array), box.lower);
}

template <typename State, typename Equation, std::size_t... Is>
auto MakeView(std::index_sequence<Is...>,
              span<SAMRAI::pdat::CellData<double>*> span,
              const Equation& equation) {
  BasicView<State> view;
  constexpr auto depths = Depths<std::decay_t<State>>(equation);
  ((get<Is>(view) =
        MakePatchDataView<std::decay_t<decltype(get<Is>(view))>::Rank()>(
            span[Is]->getArrayData())),
   ...);
  return view;
}

template <typename State, typename Equation>
auto MakeView(span<SAMRAI::pdat::CellData<double>*> span,
              const Equation& equation) {
  return MakeView(std::make_index_sequence<NumVariables<State>::value>(), span,
                  equation);
}

int GetDirection(const SAMRAI::hier::IntVector& directions);

template <typename State, typename Equation, std::size_t... Is>
auto MakeView(std::index_sequence<Is...>,
              span<SAMRAI::pdat::SideData<double>*> span,
              const Equation& equation) {
  BasicView<State> view;
  constexpr auto depths = Depths<std::decay_t<State>>(equation);
  const int dir = GetDirection(span[0]->getDirectionVector());
  ((get<Is>(view) =
        MakePatchDataView<std::decay_t<decltype(get<Is>(view))>::Rank()>(
            span[Is]->getArrayData(dir))),
   ...);
  return view;
}

template <typename State, typename Equation>
auto MakeView(span<SAMRAI::pdat::SideData<double>*> span,
              const Equation& equation) {
  return MakeView(std::make_index_sequence<NumVariables<State>::value>(), span,
                  equation);
}

} // namespace samrai
} // namespace fub

#endif
