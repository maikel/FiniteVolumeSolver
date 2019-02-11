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

#include "fub/simple_grid/SimpleSolver.hpp"
#include "fub/equations/Advection.hpp"
#include "fub/simple_grid/SimpleGrid.hpp"
#include "fub/simple_grid/boundary_condition/PeriodicBoundary.hpp"

#include "fub/flux_method/GodunovMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"
#include "fub/solver/HyperbolicSplitTimeIntegrator.hpp"

#include "fub/simple_grid/boundary_condition/PeriodicBoundary.hpp"

#include <benchmark/benchmark.h>

#include <cmath>
#include <iostream>

using namespace fub;

static double wave_packet(double x) {
  return std::exp(-100 * x * x) * std::sin(40.0 * M_PI * x);
}

static void BM_AdvanceState(benchmark::State& state) {
  std::array<double, 3> velocity{1.0};
  Advection equation(velocity);

  const std::ptrdiff_t n_cells = state.range(0);
  const int gcw = 2;
  const double dx = 1.0 / n_cells;
  SimpleGrid grid(equation, n_cells, gcw, dx);

  auto mass = grid.states().mass;
  for (std::ptrdiff_t i = 0; i < n_cells; ++i) {
    const double x = 0.5 * dx + i * dx;
    mass(i, 0, 0) = wave_packet(x - 0.5);
  }

  DimensionalSplitHyperbolicTimeIntegrator integrator(equation);
  MusclHancockMethod flux_method(equation, GodunovMethod(equation));
  SimpleSolver solver(integrator, flux_method);
  simple::PeriodicBoundary boundary(equation);

  for (auto _ : state) {
    solver.AdvanceGrid(grid, boundary, 0.9 * dx);
  }
}
BENCHMARK(BM_AdvanceState)->Range(8, 8 << 10);

BENCHMARK_MAIN();