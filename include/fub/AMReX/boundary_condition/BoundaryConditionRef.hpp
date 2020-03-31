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

namespace fub::amrex {

/// \ingroup BoundaryCondition
/// \brief This class references a BoundaryCondition object and acts as an
/// adapter such that it is enabled to be used as an AMReX boundary condition
template <typename BC, typename GriddingAlgorithm>
class BoundaryConditionRef : public ::amrex::PhysBCFunctBase {
public:
  /// @{
  /// \name Constructors

  /// \brief This adapter class is not default constructible.
  BoundaryConditionRef() = delete;

  /// \brief Constructs a reference to `condition` and stores a pointer to
  /// `grid` at refinement level `level`.
  BoundaryConditionRef(BC& condition, const GriddingAlgorithm& grid, int level);

  /// \brief This class is neither copy- nor move-assignable.
  BoundaryConditionRef(const BoundaryConditionRef&) = delete;
  BoundaryConditionRef& operator=(const BoundaryConditionRef&) = delete;

  BoundaryConditionRef(BoundaryConditionRef&&) = delete;
  BoundaryConditionRef& operator=(BoundaryConditionRef&&) = delete;
  /// @}

  /// \name Actions
  /// \brief This function fills the ghost layer of the specified MultiFab `mf`.
  void FillBoundary(::amrex::MultiFab& mf, int, int, const ::amrex::IntVect&,
                    double time_point, int) override;
};

} // namespace fub::amrex

#endif
