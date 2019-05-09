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

#ifndef FUB_AMREX_CUTCELL_RECONSTRUCTION_HPP
#define FUB_AMREX_CUTCELL_RECONSTRUCTION_HPP

#include "fub/CompleteFromCons.hpp"
#include "fub/CutCellData.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"

namespace fub {
namespace amrex {
namespace cutcell {
template <typename Eq, typename Equation = std::decay_t<Eq>>
void CompleteFromCons(Eq&& eq, const View<Complete<Equation>>& complete_view,
                      const View<const Conservative<Equation>>& cons_view,
                      const CutCellData<Equation::Rank()>& cutcell_data) {
  Complete<Equation> complete(eq);
  Conservative<Equation> cons(eq);
  const auto& flags = cutcell_data.flags;
  ForEachIndex(Box<0>(complete_view), [&](auto... is) {
    if (flags(is...).isRegular() || flags(is...).isSingleValued()) {
      std::array<std::ptrdiff_t, sizeof...(is)> cell{is...};
      Load(cons, cons_view, cell);
      CompleteFromCons(eq, complete, cons);
      Store(complete_view, complete, cell);
    }
  });
}

template <typename Equation> struct Reconstruction {
  static constexpr int Rank = Equation::Rank();

  explicit Reconstruction(const Equation& eq) : equation_{eq} {}

  template <typename Context>
  void CompleteFromCons(Context& context, PatchHandle patch, Direction dir,
                        Duration) {
    const auto tilebox = AsIndexBox<Rank>(patch.iterator->tilebox());
    View<Complete<Equation>> state =
        Subview(MakeView<BasicView<Complete<Equation>>>(context.GetData(patch),
                                                        equation_),
                tilebox);
    View<const Conservative<Equation>> scratch =
        AsCons(Subview(MakeView<BasicView<const Complete<Equation>>>(
                           context.GetScratch(patch, dir), equation_),
                       tilebox));
    ::amrex::FabType type = context.GetCutCellPatchType(patch);
    if (type == ::amrex::FabType::regular) {
      ::fub::CompleteFromCons(equation_, state, scratch);
    } else if (type == ::amrex::FabType::singlevalued) {
      CutCellData<Equation::Rank()> cutcell_data =
          context.GetCutCellData(patch, dir);
      ::fub::amrex::cutcell::CompleteFromCons(equation_, state, scratch,
                                              cutcell_data);
    }
  }

  Equation equation_;
};

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif