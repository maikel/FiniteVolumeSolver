// Copyright (c) 2021 Christian Zenker
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

#ifndef FUB_AMREX_CUTCELL_MASSFLOW_BOUNDARY_PERFECTGAS_JIRASEK_HPP
#define FUB_AMREX_CUTCELL_MASSFLOW_BOUNDARY_PERFECTGAS_JIRASEK_HPP

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"

#include "fub/Direction.hpp"
#include "fub/equations/PerfectGas.hpp"

#include <AMReX.H>

namespace fub::amrex::cutcell {
/// \ingroup BoundaryCondition
///
struct MassflowBoundary_PerfectGas_JirasekOptions {
  MassflowBoundary_PerfectGas_JirasekOptions() = default;
  MassflowBoundary_PerfectGas_JirasekOptions(const ProgramOptions& options);

  void Print(SeverityLogger& log) const;
  
  std::string channel_name{"MassflowBoundary_PerfectGas_Jirasek"};
  ::amrex::Box coarse_inner_box{};
  double required_massflow = 0.0;
  double surface_area = 1.0;
  Duration ramp_time{0.0};
  Direction dir = Direction::X;
  int side = 0;
};

/// \ingroup BoundaryCondition
///
class MassflowBoundary_PerfectGas_Jirasek {
public:
  MassflowBoundary_PerfectGas_Jirasek(
      const PerfectGas<AMREX_SPACEDIM>& eq,
      const MassflowBoundary_PerfectGas_JirasekOptions& options);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir) {
    if (dir == options_.dir) {
      FillBoundary(mf, gridding, level);
    }
  }

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
                    const ::amrex::Geometry& geom,
                    const Complete<PerfectGas<AMREX_SPACEDIM>>& state);

private:
  PerfectGas<AMREX_SPACEDIM> equation_;
  MassflowBoundary_PerfectGas_JirasekOptions options_;
};

} // namespace fub::amrex::cutcell

#endif // !FUB_AMREX_CUTCELL_MASSFLOW_BOUNDARY_PERFECTGAS_JIRASEK_HPP
