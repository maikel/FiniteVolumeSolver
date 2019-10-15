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

#ifndef FUB_AMREX_CUTCELL_BOUNDARY_CONDITION_REFLECTIVE_BOUNDARY_HPP
#define FUB_AMREX_CUTCELL_BOUNDARY_CONDITION_REFLECTIVE_BOUNDARY_HPP

#include "fub/AMReX/BoundaryCondition.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/Execution.hpp"
#include "fub/boundary_condition/ReflectiveBoundary.hpp"

namespace fub::amrex::cutcell {

template <typename Tag, typename Equation> class ReflectiveBoundary {
public:
  ReflectiveBoundary(Tag, const Equation& equation, Direction dir, int side);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration dt, const GriddingAlgorithm& grid);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    const GriddingAlgorithm& grid);

private:
  Local<Tag, Equation> equation_;
  Direction dir_;
  int side_;
};

template <typename GriddingAlgorithm>
int FindLevel_(const ::amrex::Geometry& geom,
               const GriddingAlgorithm& gridding) {
  for (int level = 0; level < gridding.GetPatchHierarchy().GetNumberOfLevels();
       ++level) {
    if (geom.Domain() ==
        gridding.GetPatchHierarchy().GetGeometry(level).Domain()) {
      return level;
    }
  }
  return -1;
}

template <typename Tag, typename Equation>
ReflectiveBoundary<Tag, Equation>::ReflectiveBoundary(Tag,
                                                      const Equation& equation,
                                                      Direction dir, int side)
    : equation_{Local<Tag, Equation>{equation}}, dir_{dir}, side_{side} {}

template <typename Tag, typename Equation>
void ReflectiveBoundary<Tag, Equation>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
    const GriddingAlgorithm& grid) {
  const int ngrow = mf.nGrow(int(dir_));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }
  static constexpr int Rank = Equation::Rank();
  const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(dir_);
  const int level_num = FindLevel_(geom, grid);
  const ::amrex::MultiFab& alphas =
      grid.GetPatchHierarchy().GetEmbeddedBoundary(level_num)->getVolFrac();
  ForEachFab(Tag(), mf, [&](const ::amrex::MFIter& mfi) {
    Complete<Equation> state(*equation_);
    Complete<Equation> reflected(*equation_);
    ::amrex::FArrayBox& fab = mf[mfi];
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    for (const ::amrex::Box& boundary : boundaries) {
      ::amrex::Box shifted =
          ::amrex::shift(boundary, int(dir_), GetSign(side_) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      if (!box_to_fill.isEmpty()) {
        auto states =
            MakeView<Complete<Equation>>(fab, *equation_, mfi.growntilebox());
        auto box = AsIndexBox<Equation::Rank()>(box_to_fill);
        ForEachIndex(Box<0>(states), [&](auto... is) {
          std::array<std::ptrdiff_t, Rank> dest{is...};
          std::array<std::ptrdiff_t, Rank> src =
              ReflectIndex(dest, box, dir_, side_);
          ::amrex::IntVect iv{
              AMREX_D_DECL(int(src[0]), int(src[1]), int(src[2]))};
          if (alpha(iv) > 0.0) {
            Load(state, states, src);
            FUB_ASSERT(state.density > 0.0);
            Reflect(reflected, state, unit, *equation_);
            Store(states, reflected, dest);
          }
        });
      }
    }
  });
}

template <typename Tag, typename Equation>
void ReflectiveBoundary<Tag, Equation>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom, Duration,
    const GriddingAlgorithm& grid) {
  return FillBoundary(mf, geom, grid);
}

template <typename Tag, typename Equation>
ReflectiveBoundary(Tag, const Equation&, Direction, int)
    ->ReflectiveBoundary<Tag, Equation>;

} // namespace fub::amrex::cutcell

#endif
