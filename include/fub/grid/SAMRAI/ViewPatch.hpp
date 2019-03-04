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

#include <boost/hana/ext/std/array.hpp>

#include <SAMRAI/pdat/ArrayData.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/CellVariable.h>

namespace fub {
namespace samrai {

template <typename MdSpan, typename T>
MdSpan MakeMdSpan(boost::hana::basic_type<MdSpan>,
                  SAMRAI::pdat::ArrayData<T>& array) {
  const SAMRAI::hier::Box& box = array.getBox();
  const int dim = box.getDim().getValue();
  const int depth = array.getDepth();
  constexpr int Rank = MdSpan::rank();
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

template <typename MdSpan, typename T>
MdSpan MakeMdSpan(boost::hana::basic_type<MdSpan>,
                  const SAMRAI::pdat::ArrayData<T>& array) {
  const SAMRAI::hier::Box& box = array.getBox();
  const int dim = box.getDim().getValue();
  const int depth = array.getDepth();
  constexpr int Rank = MdSpan::rank();
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

template <template <typename...> typename T, typename... Args>
T<remove_cvref_t<Args>...> MakeTemplate(template_t<T>,
                                        Args&&... args) {
  return T<remove_cvref_t<Args>...>{std::forward<Args>(args)...};
}

template <typename State, typename Equation,
          int N = boost::hana::length(State::ValueTypes())>
auto MakeView(boost::hana::basic_type<State>,
              nodeduce_t<span<SAMRAI::pdat::CellData<double>*, N>> span,
              const Equation&) {
  constexpr auto types = State::ValueTypes();
  const auto pointers = boost::hana::to_tuple(span);
  return boost::hana::unpack(
      boost::hana::zip(types, pointers), [](auto... args) {
        return State{
            MakeMdSpan(at_c<0>(args), at_c<1>(args)->getArrayData())...};
      });
}

int GetDirection(const SAMRAI::hier::IntVector& directions);

template <typename State, typename Equation,
          int N = boost::hana::length(State::ValueTypes())>
auto MakeView(boost::hana::basic_type<State>,
              nodeduce_t<span<SAMRAI::pdat::SideData<double>*, N>> span,
              const Equation&) {
  constexpr auto types = State::ValueTypes();
  const auto pointers = boost::hana::to_tuple(span);
  const int direction = GetDirection(span[0]->getDirectionVector());
  return boost::hana::unpack(
      boost::hana::zip(types, pointers), [&](auto... args) {
        return State{MakeMdSpan(at_c<0>(args),
                                at_c<1>(args)->getArrayData(direction))...};
      });
}

} // namespace samrai
} // namespace fub

#endif