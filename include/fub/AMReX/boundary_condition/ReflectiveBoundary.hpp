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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_REFLECTIVE_BOUNDARY_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_REFLECTIVE_BOUNDARY_HPP

#include "fub/AMReX/BoundaryCondition.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/Execution.hpp"
#include "fub/boundary_condition/ReflectiveBoundary.hpp"

namespace fub::amrex {

template <typename Tag, typename Equation> class ReflectiveBoundary {
public:
  ReflectiveBoundary(const Equation& equation, Direction dir, int side) : ReflectiveBoundary(Tag(), equation, dir, side) {}
  ReflectiveBoundary(Tag, const Equation& equation, Direction dir, int side);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration dt, const GriddingAlgorithm&);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration dt, const GriddingAlgorithm& grid, Direction dir);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom);

private:
  Local<Tag, ::fub::ReflectiveBoundary<Equation>> implementation_;
  Direction dir_;
  int side_;
};

template <typename Equation>
ReflectiveBoundary(const Equation&, Direction, int) -> ReflectiveBoundary<execution::SequentialTag, Equation>;

template <typename Tag, typename Equation>
ReflectiveBoundary<Tag, Equation>::ReflectiveBoundary(Tag,
                                                      const Equation& equation,
                                                      Direction dir, int side)
    : implementation_{::fub::ReflectiveBoundary<Equation>(equation)}, dir_{dir},
      side_{side} {}

template <typename Tag, typename Equation>
void ReflectiveBoundary<Tag, Equation>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom) {
  const int ngrow = mf.nGrow(int(dir_));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }
  ForEachFab(Tag(), mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    for (const ::amrex::Box& boundary : boundaries) {
      ::amrex::Box shifted =
          ::amrex::shift(boundary, int(dir_), GetSign(side_) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<Equation>>(
            fab, implementation_->GetEquation(), mfi.growntilebox());
        auto box = AsIndexBox<Equation::Rank()>(box_to_fill);
        implementation_->FillBoundary(states, box, dir_, side_);
      }
    }
  });
}

template <typename Tag, typename Equation>
void ReflectiveBoundary<Tag, Equation>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom, Duration,
    const GriddingAlgorithm&) {
  FillBoundary(mf, geom);
}

template <typename Tag, typename Equation>
void ReflectiveBoundary<Tag, Equation>::FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration, const GriddingAlgorithm&, Direction dir) {
                      if (dir == dir_) {
                      FillBoundary(mf, geom);
                      }
                    }

template <typename Tag, typename Equation>
ReflectiveBoundary(Tag, const Equation&, Direction, int)
    ->ReflectiveBoundary<Tag, Equation>;

} // namespace fub::amrex

#endif
