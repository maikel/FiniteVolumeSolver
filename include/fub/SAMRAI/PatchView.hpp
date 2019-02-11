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

template <typename MdSpan>
MdSpan MakeMdSpan(boost::hana::basic_type<MdSpan>,
                  const SAMRAI::hier::Patch& patch,
                  boost::hana::basic_type<SAMRAI::pdat::CellVariable<double>>,
                  int id) {
  SAMRAI::hier::PatchData* pdata = patch.getPatchData(id).get();
  auto data = static_cast<SAMRAI::pdat::CellData<double>*>(pdata);
  return MakeMdSpan(boost::hana::type_c<MdSpan>, data->getArrayData());
}

template <typename StateView>
auto ViewCellDataOnPatch(const SAMRAI::hier::Patch& patch,
                         std::array<int, StateView::size()> ids) {
  constexpr auto value_types = StateView::value_types();
  auto tids = boost::hana::to<boost::hana::tuple_tag>(ids);
  auto value_type_and_id = boost::hana::zip(value_types, tids);
  return boost::hana::unpack(value_type_and_id, [&](auto... type_and_id) {
    using namespace boost::hana::literals;
    return StateView{
        MakeMdSpan(type_and_id[0_c], patch,
                   boost::hana::type_c<SAMRAI::pdat::CellVariable<double>>,
                   type_and_id[1_c])...};
  });
}

} // namespace samrai
} // namespace fub

#endif