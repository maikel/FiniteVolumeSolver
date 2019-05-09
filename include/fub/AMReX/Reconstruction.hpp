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

#ifndef FUB_AMREX_RECONSTRUCTION_HPP
#define FUB_AMREX_RECONSTRUCTION_HPP

#include "fub/CompleteFromCons.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"

namespace fub {
namespace amrex {

template <typename Equation> struct Reconstruction {
  explicit Reconstruction(const Equation& eq) : equation_{eq} {}

  static constexpr int Rank = Equation::Rank();

  template <typename Context>
  void CompleteFromCons(Context& context, PatchHandle patch, Direction dir,
                        Duration) {
    const IndexBox<Rank> tilebox = AsIndexBox<Rank>(patch.iterator->tilebox());

    View<Complete<Equation>> state =
        Subview(MakeView<BasicView<Complete<Equation>>>(context.GetData(patch),
                                                        equation_),
                tilebox);

    View<Conservative<Equation>> scratch =
        Subview(AsCons(MakeView<BasicView<Complete<Equation>>>(
                    context.GetScratch(patch, dir), equation_)),
                tilebox);

    ::fub::CompleteFromCons(equation_, state, AsConst(scratch));
  }

  Equation equation_;
};

} // namespace amrex
} // namespace fub

#endif