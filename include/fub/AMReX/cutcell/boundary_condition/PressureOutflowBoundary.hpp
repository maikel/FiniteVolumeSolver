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

#ifndef FUB_AMREX_CUTCELL_PRESSURE_OUTFLOW_BOUNDARY_HPP
#define FUB_AMREX_CUTCELL_PRESSURE_OUTFLOW_BOUNDARY_HPP

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/Direction.hpp"
#include "fub/equations/PerfectGas.hpp"
#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace fub::amrex::cutcell {
/// \ingroup BoundaryCondition
///
struct PressureOutflowOptions {
  PressureOutflowOptions() = default;
  PressureOutflowOptions(const ProgramOptions& options);

  void Print(SeverityLogger& log);

  double outer_pressure = 101325.0;
  Direction direction = Direction::X;
  int side = 0;
};

/// \ingroup BoundaryCondition
///
class PressureOutflowBoundary {
public:
  PressureOutflowBoundary(const PerfectGas<AMREX_SPACEDIM>& eq,
                          const PressureOutflowOptions& options);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid,
                    int level, Direction dir) {
    if (dir == options_.direction) {
      FillBoundary(mf, grid, level);
    }
  }

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
                    const ::amrex::Geometry& geom);

private:
  PerfectGas<AMREX_SPACEDIM> equation_;
  PressureOutflowOptions options_;
};

} // namespace fub::amrex::cutcell

#endif // !FUB_AMREX_CUTCELL_PRESSURE_OUTFLOW_BOUNDARY_HPP
