#include "fub/equations/ShallowWater.hpp"
#include "fub/simple_grid/SimpleGrid.hpp"
#include "fub/simple_grid/SimpleSolver.hpp"
#include "fub/simple_grid/boundary_condition/PeriodicBoundary.hpp"

#include "fub/flux_method/GodunovMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"
#include "fub/solver/DimensionalSplitHyperbolicTimeIntegrator.hpp"

#include "fub/simple_grid/boundary_condition/PeriodicBoundary.hpp"

#include <benchmark/benchmark.h>

#include <cmath>
#include <iostream>

using namespace fub;

static double glocke(double x) { return std::exp(-100 * x * x); }

static void BM_AdvanceState(benchmark::State& state) {
  ShallowWater equation{};

  const std::ptrdiff_t n_cells = state.range(0);
  const int gcw = 2;
  const double dx = 1.0 / n_cells;
  SimpleGrid grid(equation, n_cells, gcw, dx);

  auto heigth = grid.states().heigth;
  auto momentum = grid.states().momentum;
  for (std::ptrdiff_t i = 0; i < n_cells; ++i) {
    const double x = 0.5 * dx + i * dx;
    heigth(i, 0, 0) = 1.0 + glocke(x - 0.5);
    momentum(i, 0, 0, 0) = 0.0;
    momentum(i, 0, 0, 1) = 0.0;
  }

  DimensionalSplitHyperbolicTimeIntegrator integrator(equation);
  // MusclHancockMethod flux_method(equation, GodunovMethod(equation));
  GodunovMethod flux_method(equation);
  SimpleSolver solver(integrator, flux_method);
  simple::PeriodicBoundary boundary(equation);

  for (auto _ : state) {
    for (int i = 0; i < 100000; ++i) 
      solver.AdvanceGrid(grid, boundary, 1e-6);
    }
  }
}
BENCHMARK(BM_AdvanceState)->Range(8, 8 << 10);

BENCHMARK_MAIN();