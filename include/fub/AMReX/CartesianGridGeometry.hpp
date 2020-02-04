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

#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include <AMReX_Geometry.H>

#include <array>
#include <fmt/ranges.h>

#ifndef FUB_AMREX_CARTESIAN_GRID_GEOMETRY_HPP
#define FUB_AMREX_CARTESIAN_GRID_GEOMETRY_HPP

namespace fub {
namespace amrex {

struct CartesianGridGeometry {
  CartesianGridGeometry() = default;
  CartesianGridGeometry(const ProgramOptions& options);

  template <typename Log> void Print(Log& log);

  std::array<int, AMREX_SPACEDIM> cell_dimensions{AMREX_D_DECL(32, 32, 32)};
  ::amrex::RealBox coordinates{
      std::array<double, AMREX_SPACEDIM>{},
      std::array<double, AMREX_SPACEDIM>{AMREX_D_DECL(1.0, 1.0, 1.0)}};
  std::array<int, AMREX_SPACEDIM> periodicity{};
};

template <typename Log> void CartesianGridGeometry::Print(Log& log) {
  std::array<double, AMREX_SPACEDIM> xlo{};
  std::array<double, AMREX_SPACEDIM> xup{};
  std::copy_n(coordinates.lo(), AMREX_SPACEDIM, xlo.data());
  std::copy_n(coordinates.hi(), AMREX_SPACEDIM, xup.data());
  BOOST_LOG(log) << fmt::format(" - cell_dimensions = {{{}}}",
                                fmt::join(cell_dimensions, ", "));
  BOOST_LOG(log) << fmt::format(" - x_lower = {{{}}}", fmt::join(xlo, ", "));
  BOOST_LOG(log) << fmt::format(" - x_upper = {{{}}}", fmt::join(xup, ", "));
  BOOST_LOG(log) << fmt::format(" - periodicity = {{{}}}",
                                fmt::join(periodicity, ", "));
}

::amrex::Geometry GetCoarseGeometry(const CartesianGridGeometry& grid_geometry);

::amrex::RealBox DomainAroundCenter(const ::amrex::RealArray& x,
                                    const ::amrex::RealArray& rx);

::amrex::Box BoxWhichContains(const ::amrex::RealBox& xbox,
                              const ::amrex::Geometry& geom);

} // namespace amrex
} // namespace fub

#endif
