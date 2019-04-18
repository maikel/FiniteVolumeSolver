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

#include "fub/grid/AMReX/BoundaryCondition.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"

namespace fub {
namespace amrex {

BoundaryCondition::BoundaryCondition(Function f, const ::amrex::Geometry& geom,
                                     int level)
    : function_{f}, geom_{geom}, level_num_{level} {}

void BoundaryCondition::FillBoundary(::amrex::MultiFab& mf, int, int,
                                     double time_point, int) {
  if (geom_.isAllPeriodic())
    return;

  //! create a grown domain box containing valid + periodic cells
  const ::amrex::Box& domain = geom_.Domain();
  ::amrex::Box gdomain = ::amrex::convert(domain, mf.boxArray().ixType());
  const ::amrex::IntVect& ngrow = mf.nGrowVect();
  std::array<int, AMREX_SPACEDIM + 1> gcw{};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    gcw[static_cast<std::size_t>(i)] = ngrow[i];
    if (geom_.isPeriodic(i)) {
      gdomain.grow(i, ngrow[i]);
    }
  }
  for (::amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    ::amrex::Box box = mfi.growntilebox();
    if (!gdomain.contains(box)) {
      PatchHandle patch{level_num_, &mfi};
      PatchDataView<double, AMREX_SPACEDIM + 1> data =
          MakePatchDataView(mf[mfi]);
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (ngrow[dir] && box.smallEnd(dir) < gdomain.smallEnd(dir)) {
          function_(data, patch, Location{static_cast<std::size_t>(dir), 0}, ngrow[dir],
                    Duration(time_point));
        }
        if (ngrow[dir] && gdomain.bigEnd(dir) < box.bigEnd(dir)) {
          function_(data, patch, Location{static_cast<std::size_t>(dir), 1}, ngrow[dir],
                    Duration(time_point));
        }
      }
    }
  }
}

} // namespace amrex
} // namespace fub
