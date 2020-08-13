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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_REF_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_REF_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"

#include <AMReX_PhysBCFunct.H>

namespace fub::amrex {

/// \ingroup BoundaryCondition
///
/// \brief This class references a BoundaryCondition object and acts as an
/// adapter such that it is enabled to be used as an AMReX boundary condition
template <typename BC, typename GriddingAlgorithm>
class BoundaryConditionRef {
public:
  /// @{
  /// \name Constructors

  /// \brief This adapter class is not default constructible.
  BoundaryConditionRef() = delete;

  /// \brief Constructs a reference to `condition` and stores a pointer to
  /// `grid` at refinement level `level`.
  ///
  /// \param[in] condition The condition that will be invoked
  ///
  /// \param[in] grid The reference gridding algorithm can be used if a Geometry
  /// object or boundary conditions are needed
  ///
  /// \param[in] level The refinement level for which this boundary conditions
  /// will be applied to.
  ///
  /// \throw Nothing.
  BoundaryConditionRef(BC& condition, const GriddingAlgorithm& grid, int level);
  /// @}

  /// \name Actions

  /// \brief This function fills the ghost layer of the specified MultiFab `mf`.
  ///
  /// \throw Throws any exception that is thrown by the reference boundary
  /// condition.
  void operator() (::amrex::MultiFab& mf, int dcomp, int ncomp,
                   ::amrex::IntVect const& nghost, ::amrex::Real time,
                   int bccomp);

  BC* pointer;
  const GriddingAlgorithm* gridding;
  int level;
};

template <typename BC, typename GriddingAlgorithm>
BoundaryConditionRef<BC, GriddingAlgorithm>::BoundaryConditionRef(
    BC& condition, const GriddingAlgorithm& grid, int lvl)
    : pointer{&condition}, gridding{&grid}, level{lvl} {}

template <typename BC, typename GriddingAlgorithm>
void BoundaryConditionRef<BC, GriddingAlgorithm>::operator() (::amrex::MultiFab& mf, int, int,
                    ::amrex::IntVect const&, ::amrex::Real,
                    int) {
 if (pointer && gridding) {
   pointer->FillBoundary(mf, *gridding, level);
 }
}

} // namespace fub::amrex

#endif
