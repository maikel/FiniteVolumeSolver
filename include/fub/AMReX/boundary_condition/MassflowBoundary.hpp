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

#ifndef FUB_AMREX_MASSFLOW_BOUNDARY_HPP
#define FUB_AMREX_MASSFLOW_BOUNDARY_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/Direction.hpp"
#include "fub/equations/IdealGasMix.hpp"

#include <AMReX.H>

namespace fub::amrex {

/// \ingroup BoundaryCondition
///
/// \brief This boundary models an inflow boundary with constant mean mass flow.
template <int Rank> class MassflowBoundary {
public:
  MassflowBoundary(const IdealGasMix<Rank>& eq,
                   const ::amrex::Box& coarse_inner_box,
                   double required_massflow, double surface_area, Direction dir,
                   int side);

  MassflowBoundary(const std::string& name, const IdealGasMix<Rank>& eq,
                   const ::amrex::Box& coarse_inner_box,
                   double required_massflow, double surface_area, Direction dir,
                   int side);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid,
                    int level, Direction dir);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
                    const ::amrex::Geometry& geom,
                    const Complete<IdealGasMix<AMREX_SPACEDIM>>& state);

private:
  IdealGasMix<Rank> equation_;
  ::amrex::Box coarse_inner_box_;
  double required_massflow_;
  double surface_area_;
  Direction dir_;
  int side_;
};

} // namespace fub::amrex

#endif // !FUB_AMREX_CUTCELL_ISENTROPIC_PRESSURE_BOUNDARY_HPP
