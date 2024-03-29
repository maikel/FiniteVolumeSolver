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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_TRANSMISSIVE_BOUNDARY_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_TRANSMISSIVE_BOUNDARY_HPP

#include "fub/Direction.hpp"

#include "fub/AMReX/GriddingAlgorithm.hpp"

namespace fub::amrex {

/// \ingroup BoundaryCondition
///
/// \brief This class copies the inner grid values to the boundary.
///
/// In case of the Euler equations this is equivalent to model a supersonic
/// boundary since all signals leave domain.
///
/// This condition also fills the physical boundary in one direction at only one
/// side. You can construct an object as in the following example:
///
/// ```cpp
/// // ...
/// fub::amrex:TransmissiveBoundary boundary{fub::Direction::X, 0};
/// // ...
/// ```
struct TransmissiveBoundary {
  Direction dir;
  int side;

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid,
                    int level, Direction dir);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom);
};

} // namespace fub::amrex

#endif
