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

#include "fub/AMReX/CartesianGridGeometry.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"

namespace fub::amrex {

CartesianGridGeometry::CartesianGridGeometry(const ProgramOptions& options) {
  std::array<int, 3> cells{32, 32, 32};
  cells = GetOptionOr(options, "cell_dimensions", cells);
  std::copy_n(cells.data(), AMREX_SPACEDIM, cell_dimensions.data());
  coordinates = GetOptionOr(options, "coordinates", coordinates);
  std::array<int, 3> p{0, 0, 0};
  p = GetOptionOr(options, "periodicity", p);
  std::copy_n(p.data(), AMREX_SPACEDIM, periodicity.data());
}

::amrex::Box BoxWhichContains(const ::amrex::RealBox& xbox,
                              const ::amrex::Geometry& geom) {
  ::amrex::Box domain = geom.Domain();
  ::amrex::IntVect lo = domain.smallEnd();
  ::amrex::IntVect up = domain.bigEnd();
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    for (int i = domain.smallEnd(d); i < domain.bigEnd(d); ++i) {
      const double x = geom.CellCenter(i, d);
      if (x < xbox.lo(d)) {
        lo[d] = std::max(lo[d], i);
      }
      if (x > xbox.hi(d)) {
        up[d] = std::min(up[d], i);
      }
    }
  }
  return ::amrex::Box{lo, up};
}

::amrex::Geometry
GetCoarseGeometry(const CartesianGridGeometry& grid_geometry) {
  ::amrex::Box domain{{},
                      {AMREX_D_DECL(grid_geometry.cell_dimensions[0] - 1,
                                    grid_geometry.cell_dimensions[1] - 1,
                                    grid_geometry.cell_dimensions[2] - 1)}};
  return ::amrex::Geometry(domain, &grid_geometry.coordinates, -1,
                           grid_geometry.periodicity.data());
}

::amrex::RealBox DomainAroundCenter(const ::amrex::RealArray& x,
                                    const ::amrex::RealArray& rx) {
  return ::amrex::RealBox{
      {AMREX_D_DECL(x[0] - rx[0], x[1] - rx[1], x[2] - rx[2])},
      {AMREX_D_DECL(x[0] + rx[0], x[1] + rx[1], x[2] + rx[2])}};
}

} // namespace fub::amrex
