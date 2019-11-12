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

#include <SAMRAI/hier/Patch.h>
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
  return mdspan<T, Rank>(pointer, extents);
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
  return mdspan<const T, Rank>(pointer, extents);
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
PatchDataView<const T, Rank>
MakePatchDataView(const SAMRAI::pdat::ArrayData<T>& array) {
  IndexBox<Rank> box = AsIndexBox<Rank>(array.getBox());
  return PatchDataView<const T, Rank>(MakeMdSpan<Rank, T>(array), box.lower);
}

template <int Rank, typename T>
PatchDataView<T, Rank> MakePatchDataView(SAMRAI::pdat::ArrayData<T>& array) {
  IndexBox<Rank> box = AsIndexBox<Rank>(array.getBox());
  return PatchDataView<T, Rank>(MakeMdSpan<Rank, T>(array), box.lower);
}

template <typename PatchData>
std::enable_if_t<std::is_pointer_v<PatchData>, void>
GetPatchData(span<PatchData> patch_datas, SAMRAI::hier::Patch& patch,
             span<const int> data_ids) {
  FUB_ASSERT(data_ids.size() == patch_datas.size());
  std::transform(data_ids.begin(), data_ids.end(), patch_datas.begin(),
                 [&patch](int id) {
                   PatchData pointer_to_data =
                       dynamic_cast<PatchData>(patch.getPatchData(id).get());
                   FUB_ASSERT(pointer_to_data);
                   return pointer_to_data;
                 });
}

template <typename PatchData>
std::enable_if_t<std::is_pointer_v<PatchData>, void>
GetPatchData(span<PatchData> patch_datas, const SAMRAI::hier::Patch& patch,
             span<const int> data_ids) {
  FUB_ASSERT(data_ids.size() == patch_datas.size());
  std::transform(data_ids.begin(), data_ids.end(), patch_datas.begin(),
                 [&patch](int id) {
                   PatchData pointer_to_data =
                       dynamic_cast<PatchData>(patch.getPatchData(id).get());
                   FUB_ASSERT(pointer_to_data);
                   return pointer_to_data;
                 });
}

template <typename State>
BasicView<State> MakeView(span<SAMRAI::pdat::SideData<double>*> span,
                          const typename State::Equation& /* equation */,
                          Direction dir) {
  BasicView<State> view;
  int i = 0;
  const int dir_value = static_cast<int>(dir);
  ForEachVariable(
      [&](auto& variable) {
        using type = typename remove_cvref<decltype(variable)>::type;
        constexpr int Rank = type::rank();
        variable = MakePatchDataView<Rank>(span[i]->getArrayData(dir_value));
        ++i;
      },
      view);
  return view;
}

template <typename State>
BasicView<State> MakeView(span<SAMRAI::pdat::CellData<double>*> span,
                          const typename State::Equation& /* equation */) {
  BasicView<State> view;
  int i = 0;
  static constexpr int Rank = State::Equation::Rank();
  using T = std::conditional_t<std::is_const_v<State>, const double, double>;
  ForEachVariable(
      overloaded{[&](PatchDataView<T, Rank>& variable) {
                   variable = MakePatchDataView<Rank>(span[i]->getArrayData());
                   ++i;
                 },
                 [&](PatchDataView<T, Rank + 1>& variable) {
                   variable =
                       MakePatchDataView<Rank + 1>(span[i]->getArrayData());
                   ++i;
                 }},
      view);
  return view;
}

template <typename State>
BasicView<const State>
MakeView(span<SAMRAI::pdat::CellData<double> const*> span,
         const typename State::Equation& /* equation */) {
  BasicView<const State> view;
  int i = 0;
  static constexpr int Rank = State::Equation::Rank();
  ForEachVariable(
      overloaded{[&](PatchDataView<const double, Rank>& variable) {
                   variable = MakePatchDataView<Rank>(span[i]->getArrayData());
                   ++i;
                 },
                 [&](PatchDataView<const double, Rank + 1>& variable) {
                   variable =
                       MakePatchDataView<Rank + 1>(span[i]->getArrayData());
                   ++i;
                 }},
      view);
  return view;
}

template <typename State, typename Range>
View<State> MakeView(Range&& span, const typename State::Equation& equation,
                     const IndexBox<State::Equation::Rank()>& box) {
  return Subview(MakeView<State>(std::forward<Range>(span), equation), box);
}

template <typename State, typename Range>
View<State> MakeView(Range&& span, const typename State::Equation& equation,
                     Direction dir,
                     const IndexBox<State::Equation::Rank()>& box) {
  return Subview(MakeView<State>(std::forward<Range>(span), equation, dir),
                 box);
}

int GetDirection(const SAMRAI::hier::IntVector& directions);

// template <typename State, typename Equation>
// auto MakeView(span<SAMRAI::pdat::SideData<double>*> span,
//              const Equation& equation) {
//  return MakeView(std::make_index_sequence<NumVariables<State>::value>(),
//  span,
//                  equation);
//}

} // namespace samrai
} // namespace fub

#endif
