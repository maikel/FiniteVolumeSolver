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

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/Execution.hpp"
#include "fub/boundary_condition/ReflectiveBoundary.hpp"

namespace fub::amrex::cutcell {

/// \ingroup BoundaryCondition
///
template <typename Tag, typename Equation> class ReflectiveBoundary {
public:
  ReflectiveBoundary(Tag, const Equation& equation, Direction dir, int side);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);

private:
  Local<Tag, Equation> equation_;
  Direction dir_;
  int side_;
};

template <typename Tag, typename Equation>
ReflectiveBoundary<Tag, Equation>::ReflectiveBoundary(Tag,
                                                      const Equation& equation,
                                                      Direction dir, int side)
    : equation_{Local<Tag, Equation>{equation}}, dir_{dir}, side_{side} {}

template <typename Tag, typename Equation>
void ReflectiveBoundary<Tag, Equation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& gridding, int level,
    Direction dir) {
  if (dir == dir_) {
    FillBoundary(mf, gridding, level);
  }
}

template <typename Tag, typename Equation>
void ReflectiveBoundary<Tag, Equation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& grid, int level) {
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  const int ngrow = mf.nGrow(int(dir_));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }
  static constexpr int Rank = Equation::Rank();
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(dir_);
  const ::amrex::MultiFab& alphas =
      grid.GetPatchHierarchy().GetEmbeddedBoundary(level)->getVolFrac();
  ForEachFab(Tag(), mf, [&](const ::amrex::MFIter& mfi) {
    Complete<Equation> state(*equation_);
    Complete<Equation> reflected(*equation_);
    Complete<Equation> zeros(*equation_);
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
        auto box = AsIndexBox<Rank>(box_to_fill);
        ForEachIndex(box, [&](auto... is) {
          std::array<std::ptrdiff_t, sRank> dest{is...};
          std::array<std::ptrdiff_t, sRank> src =
              ReflectIndex(dest, box, dir_, side_);
          ::amrex::IntVect iv{
              AMREX_D_DECL(int(src[0]), int(src[1]), int(src[2]))};
          ::amrex::IntVect dest_iv{
              AMREX_D_DECL(int(dest[0]), int(dest[1]), int(dest[2]))};
          if (alpha(dest_iv) > 0.0 && alpha(iv) > 0.0) {
            Load(state, states, src);
            FUB_ASSERT(state.density > 0.0);
            Reflect(reflected, state, unit, *equation_);
            Store(states, reflected, dest);
          } else {
            Store(states, zeros, dest);
          }
        });
      }
    }
  });
}

template <typename Tag, typename Equation>
ReflectiveBoundary(Tag, const Equation&, Direction, int)
    ->ReflectiveBoundary<Tag, Equation>;

} // namespace fub::amrex::cutcell

#endif
