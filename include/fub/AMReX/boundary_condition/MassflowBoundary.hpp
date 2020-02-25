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

#include <boost/log/common.hpp>
#include <boost/log/core.hpp>
#include <boost/log/sources/channel_logger.hpp>

#include <AMReX.H>

namespace fub::amrex {
struct MassflowBoundaryOptions {
  MassflowBoundaryOptions() = default;
  MassflowBoundaryOptions(const ProgramOptions& options);

  std::string channel_name{"MassflowBoundary"};
  ::amrex::Box coarse_inner_box{};
  double required_massflow = 0.0;
  double surface_area = 0.0;
  Direction dir = Direction::X;
  int side = 0;
};

template <int Rank> class MassflowBoundary {
public:
  MassflowBoundary(const IdealGasMix<Rank>& eq,
                   const MassflowBoundaryOptions& options);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration dt, const GriddingAlgorithm& grid);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration dt, const GriddingAlgorithm& grid, Direction dir);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    const Complete<IdealGasMix<Rank>>& state);

private:
  IdealGasMix<Rank> equation_;
  MassflowBoundaryOptions options_;
};

extern template class MassflowBoundary<1>;
extern template class MassflowBoundary<2>;
extern template class MassflowBoundary<3>;

} // namespace fub::amrex

#endif // !FUB_AMREX_MASSFLOW_BOUNDARY_HPP
