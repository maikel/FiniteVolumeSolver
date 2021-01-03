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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_CONSTANT_BOUNDARY_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_CONSTANT_BOUNDARY_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/Direction.hpp"
#include "fub/ForEach.hpp"

#include <AMReX_MultiFab.H>

namespace fub::amrex {

/// \ingroup BoundaryCondition
///
template <typename Equation> struct ConstantBoundary {
  Direction dir;
  int side;
  Equation equation;
  Complete<Equation> state;

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid,
                    int level) {
    const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    FillBoundary(mf, geom);
  }

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir) {
    if (dir == this->dir) {
      FillBoundary(mf, gridding, level);
    }
  }

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom) {
    int idir = static_cast<int>(dir);
    if (!geom.isPeriodic(idir)) {
      const int ngrow = mf.nGrow(idir);
      ::amrex::Box fully_grown_box = geom.growNonPeriodicDomain(ngrow);
      ::amrex::Box grown_dir_box = ::amrex::grow(fully_grown_box, idir, -ngrow);
      ::amrex::BoxList boundaries = ::amrex::complementIn(
          fully_grown_box, ::amrex::BoxList{grown_dir_box});
      if (boundaries.isEmpty()) {
        return;
      }
      ForEachFab(fub::execution::openmp, mf, [&](const ::amrex::MFIter& mfi) {
        ::amrex::FArrayBox& fab = mf[mfi];
        for (const ::amrex::Box& boundary : boundaries) {
          ::amrex::Box shifted =
              ::amrex::shift(boundary, idir, GetSign(side) * ngrow);
          if (!geom.Domain().intersects(shifted)) {
            continue;
          }
          ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
          if (!box_to_fill.isEmpty()) {
            auto states =
                MakeView<Complete<Equation>>(fab, equation, box_to_fill);
            ForEachIndex(Box<0>(states),
                         [&](auto... is) { Store(states, state, {is...}); });
          }
        }
      });
    }
  }
};

} // namespace fub::amrex

#endif