// Copyright (c) 2021 Maikel Nadolski
// Copyright (c) 2021 Christian Zenker
// Copyright (c) 2021 Rupert Klein
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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_TURBINE_PLENUM_BOUNDARY_CONDITION_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_TURBINE_PLENUM_BOUNDARY_CONDITION_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/equations/PerfectGasMix.hpp"
#include "fub/equations/perfect_gas_mix/PlenaControl.hpp"

namespace fub::amrex {

class TurbinePlenumBoundaryCondition {
public:
  using ControlState = fub::perfect_gas_mix::gt::ControlState;

  TurbinePlenumBoundaryCondition(
      const PerfectGasMix<1>& equation,
      std::shared_ptr<const ControlState> control_state)
      : equation_{equation}, control_state_{control_state} {}

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid,
                    int level, Direction dir) {
    if (dir == Direction::X) {
      FillBoundary(mf, grid, level);
    }
  }

private:
  PerfectGasMix<1> equation_;
  std::shared_ptr<const ControlState> control_state_;
};

} // namespace fub::amrex

#endif