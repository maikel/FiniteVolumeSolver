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

#include "fub/AMReX/boundary_condition/TransmissiveBoundary.hpp"
#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

namespace fub::amrex {

namespace {
int GetSign(int side) { return (side == 0) - (side == 1); }

::amrex::IntVect MapToSrc(const ::amrex::IntVect& dest,
                          const ::amrex::Geometry& geom, int side,
                          Direction dir) {
  const int boundary = (side == 0) * geom.Domain().smallEnd(int(dir)) +
                       (side == 1) * geom.Domain().bigEnd(int(dir));
  const int distance = dest[int(dir)] - boundary;
  const int sign = int((distance > 0) - (distance < 0));
  ::amrex::IntVect src{dest};
  src[int(dir)] -= 2 * distance - sign;
  return src;
}

} // namespace

void TransmissiveBoundary::FillBoundary(::amrex::MultiFab& mf,
                                        const ::amrex::Geometry& geom, Duration,
                                        const GriddingAlgorithm&) {
  const int ngrow = mf.nGrow(int(dir));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
  for (::amrex::MFIter mfi(mf, true); mfi.isValid(); ++mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    for (const ::amrex::Box& boundary : boundaries) {
      ::amrex::Box shifted =
          ::amrex::shift(boundary, int(dir), GetSign(side) * ngrow);
      if (!geom.Domain().contains(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = fab.box() & boundary;
      if (!box_to_fill.isEmpty()) {
        const int ncomp = fab.nComp();
        for (int c = 0; c < ncomp; ++c) {
          ForEachIndex(box_to_fill, [this, c, &fab, &geom](auto... is) {
            ::amrex::IntVect dest{int(is)...};
            ::amrex::IntVect src = MapToSrc(dest, geom, side, dir);
            fab(dest, c) = fab(src, c);
          });
        }
      }
    }
  }
}

} // namespace fub::amrex
