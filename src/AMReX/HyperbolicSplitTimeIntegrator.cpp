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

#include "fub/AMReX/HyperbolicSplitTimeIntegrator.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/ForEach.hpp"

namespace fub::amrex {

void ForwardIntegrator::UpdateConservatively(::amrex::MultiFab& dest,
                                             const ::amrex::MultiFab& src,
                                             const ::amrex::MultiFab& fluxes,
                                             const ::amrex::Geometry& geom,
                                             Duration dt, Direction dir) {
  const int n_cons = fluxes.nComp();
  const double dx = geom.CellSize(int(dir));
  const double lambda = dt.count() / dx;
  for (::amrex::MFIter mfi(dest, true); mfi.isValid(); ++mfi) {
    ::amrex::FArrayBox& next = dest[mfi];
    const ::amrex::FArrayBox& prev = src[mfi];
    const ::amrex::FArrayBox& flux = fluxes[mfi];
    const ::amrex::Box& box = mfi.tilebox();
    IndexBox<AMREX_SPACEDIM> ibox =
        Grow(AsIndexBox<AMREX_SPACEDIM>(box), dir, {1, 1});
    for (int c = 0; c < n_cons; ++c) {
      ForEachIndex(ibox, [c, lambda, dir, &next, &prev,
                          &flux](AMREX_D_DECL(std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k)) {
        const ::amrex::IntVect iv{AMREX_D_DECL(int(i), int(j), int(k))};
        const ::amrex::IntVect left = iv;
        ::amrex::IntVect right = iv;
        right.shift(int(dir), 1);
        next(iv, c) = prev(iv, c) + lambda * (flux(left, c) - flux(right, c));
      });
    }
  }
}
} // namespace fub::amrex
