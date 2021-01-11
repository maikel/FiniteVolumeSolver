// Copyright (c) 2020 Maikel Nadolski
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

#ifndef FUB_AMREX_CUTCELL_BOUNDARY_CONDITION_REFLECTIVE_BOUNDARY2_HPP
#define FUB_AMREX_CUTCELL_BOUNDARY_CONDITION_REFLECTIVE_BOUNDARY2_HPP

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/Execution.hpp"
#include "fub/boundary_condition/ReflectiveBoundary.hpp"

namespace fub::amrex::cutcell {

/// \ingroup BoundaryCondition
///
/// \brief This boundary condition provides wall boundary which acts only on a
/// specified subset of ghost cells.
template <typename Tag, typename Equation> class ReflectiveBoundary2 {
public:
  /// \brief Constructs the boundary condition with respective execution tag.
  ReflectiveBoundary2(Tag, const Equation& equation, Direction dir, int side,
                      const ::amrex::Box& boundary_section);

  /// \brief Delegates the construction to the tag constructor.
  ReflectiveBoundary2(const Equation& equation, Direction dir, int side,
                      const ::amrex::Box& boundary_section)
      : ReflectiveBoundary2(Tag(), equation, dir, side, boundary_section) {}

  /// \brief Fill the boundary section with reflected states. The reflected
  /// state is taken from a mirrored index by given direction and side.
  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  /// \brief Conditionally fill the boundary section with reflected states, if
  /// dir == dir_. The reflected state is taken from a mirrored index by given
  /// direction and side.
  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);

private:
  /// \brief The equation defines how to reflect states.
  ///
  /// This needs to be Local'ized because the implementation might use OpenMP.
  Local<Tag, Equation> equation_;
  Direction dir_;
  int side_;
  ::amrex::Box boundary_section_;
};

template <typename Tag, typename Equation>
ReflectiveBoundary2<Tag, Equation>::ReflectiveBoundary2(
    Tag, const Equation& equation, Direction dir, int side,
    const ::amrex::Box& boundary_section)
    : equation_{Local<Tag, Equation>{equation}}, dir_{dir}, side_{side},
      boundary_section_{boundary_section} {}

template <typename Tag, typename Equation>
void ReflectiveBoundary2<Tag, Equation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& gridding, int level,
    Direction dir) {
  if (dir == dir_) {
    FillBoundary(mf, gridding, level);
  }
}

template <typename Tag, typename Equation>
void ReflectiveBoundary2<Tag, Equation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& grid, int level) {
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
    ::amrex::Box box_to_fill = mfi.growntilebox() & boundary_section_;
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
  });
}

template <typename Tag, typename Equation>
ReflectiveBoundary2(Tag, const Equation&, Direction, int, const ::amrex::Box&)
    -> ReflectiveBoundary2<Tag, Equation>;

template <typename Equation>
ReflectiveBoundary2(const Equation&, Direction, int, const ::amrex::Box&)
    -> ReflectiveBoundary2<execution::SequentialTag, Equation>;

} // namespace fub::amrex::cutcell

#endif
