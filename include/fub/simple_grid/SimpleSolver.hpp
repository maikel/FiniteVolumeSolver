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

#ifndef FUB_SIMPLE_SOLVER_HPP
#define FUB_SIMPLE_SOLVER_HPP

#include "fub/simple_grid/SimpleGrid.hpp"

namespace fub {
template <typename Integrator, typename FluxMethod> class SimpleSolver {
public:
  using Equation = typename Integrator::Equation;

  SimpleSolver(Integrator integrator, FluxMethod flux_method)
      : integrator_{integrator}, flux_method_{flux_method} {}

  template <typename BoundaryCondition>
  double ComputeStableDt(SimpleGrid<Equation>& grid,
                         BoundaryCondition& boundary_condition) {
    StateView<const double, Equation> states = grid.states();
    StateView<double, Equation> scratch = grid.scratch();
    const int fill_width = grid.ghost_cell_width();
    CopyToInnerRegion(scratch, states, fill_width);
    const double t = grid.time_point();
    boundary_condition.SetPhysicalBoundaryConditions(scratch, t, fill_width);
    return flux_method_.ComputeStableDtOnSpans(scratch, grid.dx());
  }

  template <typename BoundaryCondition>
  void AdvanceGrid(SimpleGrid<Equation>& grid,
                   BoundaryCondition& boundary_condition, double dt) {
    StateView<double, Equation> states = grid.states();
    StateView<double, Equation> scratch = grid.scratch();

    // Copy states into the inner of scratch
    const int fill_width = grid.ghost_cell_width();
    CopyToInnerRegion(scratch, states, fill_width);
    const double t = grid.time_point();
    boundary_condition.SetPhysicalBoundaryConditions(scratch, t, fill_width);

    // Compute fluxes with states from scratch
    const double dx = grid.dx();
    ConsView<double, Equation> fluxes = grid.fluxes();
    flux_method_.ComputeNumericFluxesOnSpans(fluxes, scratch, dt, dx);

    // Advance time on store results in states
    auto inner = ViewInnerRegion(scratch, fill_width);
    integrator_.AdvanceTimeOnSpans(states, inner, fluxes, dt, dx);

    grid.time_point(t + dt);
  }

private:
  Integrator integrator_;
  FluxMethod flux_method_;
};

} // namespace fub

#endif