#include "fub/simple_grid/SimpleSolver.hpp"
#include "fub/equations/Burgers.hpp"
#include "fub/simple_grid/SimpleGrid.hpp"
#include "fub/simple_grid/boundary_condition/PeriodicBoundary.hpp"

#include "fub/flux_method/GodunovMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"
#include "fub/solver/HyperbolicSplitIntegrator.hpp"

#include "fub/simple_grid/boundary_condition/PeriodicBoundary.hpp"

#include <benchmark/benchmark.h>

#include <cmath>
#include <iostream>

using namespace fub;

static double wave_packet(double x) {
  return std::exp(-100 * x * x) * std::sin(40.0 * M_PI * x);
}

static void BM_AdvanceState(benchmark::State& state) {
  Burgers equation{};

  const std::ptrdiff_t n_cells = state.range(0);
  const int gcw = 2;
  const double dx = 1.0 / n_cells;
  SimpleGrid grid(equation, n_cells, gcw, dx);

  auto mass = grid.states().mass;
  for (std::ptrdiff_t i = 0; i < n_cells; ++i) {
    const double x = 0.5 * dx + i * dx;
    mass(i, 0, 0) = wave_packet(x - 0.5);
  }

  HyperbolicSplitIntegrator integrator(equation);
  MusclHancockMethod flux_method(equation, GodunovMethod(equation));
  // GodunovMethod flux_method(equation);
  SimpleSolver solver(integrator, flux_method);
  simple::PeriodicBoundary boundary(equation);

  for (auto _ : state) {
    solver.AdvanceGrid(grid, boundary, 0.9 * dx);
  }
}
BENCHMARK(BM_AdvanceState)->Range(8, 8 << 10);

BENCHMARK_MAIN();