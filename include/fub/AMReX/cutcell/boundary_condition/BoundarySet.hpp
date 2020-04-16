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

#ifndef FUB_AMREX_CUTCELL_BOUNDARY_CONDITION_BOUNDARY_SET_HPP
#define FUB_AMREX_CUTCELL_BOUNDARY_CONDITION_BOUNDARY_SET_HPP

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"

#include <AMReX_MultiFab.H>

#include <vector>

namespace fub::amrex::cutcell {

/// \ingroup BoundaryCondition
///
/// \brief This class manages a list of boundary conditions which get applied in
/// the order in which they are stored.
struct BoundarySet {
  /// \brief Apply each condition on the specified MultiFab `mf`.
  ///
  /// \param[in,out] mf The MultiFab that will be filled.
  ///
  /// \param[in] gridding The GriddingAlgorithm will serve as a context.
  ///
  /// \param[in] level The level number which is associated with the specified
  /// MultiFab `mf`.
  ///
  /// \throw Exceptions will be propogated to the user if any of the stored
  /// boundary conditions throws an exception.
  ///
  /// \pre `mf` is defined with the same (BoxArray, DistributionMapping) as
  /// defined in the gridding algorithm at refinement level `level`.
  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  /// \brief Apply each condition on the specified MultiFab `mf` but only in
  /// direction `dir`.
  ///
  /// \param[in,out] mf The MultiFab that will be filled.
  ///
  /// \param[in] gridding The GriddingAlgorithm will serve as a context.
  ///
  /// \param[in] level The level number which is associated with the specified
  /// MultiFab `mf`.
  ///
  /// \param[in] dir The direction in which the physical boundary shall be
  /// filled.
  ///
  /// \throw Exceptions will be propogated to the user if any of the stored
  /// boundary conditions throws an exception.
  ///
  /// \pre `mf` is defined with the same (BoxArray, DistributionMapping) as
  /// defined in the gridding algorithm at refinement level `level`.
  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);

  /// \brief This vector stores all boundary conditions.
  std::vector<AnyBoundaryCondition> conditions;
};

} // namespace fub::amrex::cutcell

#endif
