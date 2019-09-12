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

#ifndef FUB_AMREX_INITIAL_DATA_CONSTANT_DATA_HPP
#define FUB_AMREX_INITIAL_DATA_CONSTANT_DATA_HPP

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/State.hpp"

namespace fub::amrex {

template <typename Equation> class ConstantData {
public:
  ConstantData(Equation eq, Complete<Equation> state);

  void InitializeData(::amrex::MultiFab& mf,
                      const ::amrex::Geometry& geom) const;

private:
  Equation equation_;
  Complete<Equation> state_;
};

template <typename Equation>
ConstantData<Equation>::ConstantData(Equation eq, Complete<Equation> state)
    : equation_{std::move(eq)}, state_{std::move(state)} {}

template <typename Equation>
void ConstantData<Equation>::InitializeData(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& /* geom */) const {
  ForEachFab(execution::openmp, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    auto view = MakeView<Complete<Equation>>(fab, equation_, mfi.tilebox());
    ForEachIndex(Box<0>(view),
                 [&](auto... is) { Store(view, state_, {is...}); });
  });
}

} // namespace fub::amrex

#endif
