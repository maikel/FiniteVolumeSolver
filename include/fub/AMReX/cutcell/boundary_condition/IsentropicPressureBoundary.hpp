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

#ifndef FUB_AMREX_CUTCELL_ISENTROPIC_PRESSURE_BOUNDARY_HPP
#define FUB_AMREX_CUTCELL_ISENTROPIC_PRESSURE_BOUNDARY_HPP

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/Direction.hpp"
#include "fub/equations/IdealGasMix.hpp"

#include <boost/log/sources/severity_channel_logger.hpp>
#include <boost/log/trivial.hpp>

#include <AMReX.H>

namespace fub::amrex::cutcell {

class IsentropicPressureBoundary {
public:
  IsentropicPressureBoundary(const std::string& name,
                             const IdealGasMix<AMREX_SPACEDIM>& eq,
                             const ::amrex::Box& coarse_inner_box,
                             double outer_pressure, Direction dir, int side);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration dt, const GriddingAlgorithm& grid);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
                    const ::amrex::Geometry& geom,
                    const Complete<IdealGasMix<AMREX_SPACEDIM>>& state);

private:
  boost::log::sources::severity_channel_logger<
      boost::log::trivial::severity_level>
      log_;
  IdealGasMix<AMREX_SPACEDIM> equation_;
  ::amrex::Box coarse_inner_box_;
  double outer_pressure_;
  Direction dir_;
  int side_;
};

} // namespace fub::amrex::cutcell

#endif // !FUB_AMREX_CUTCELL_ISENTROPIC_PRESSURE_BOUNDARY_HPP