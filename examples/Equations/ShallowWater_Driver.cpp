#include "fub/equations/ShallowWater.hpp"
#include "fub/simple_grid/SimpleGrid.hpp"
#include "fub/simple_grid/SimpleSolver.hpp"
#include "fub/simple_grid/boundary_condition/PeriodicBoundary.hpp"

#include "fub/flux_method/GodunovMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"
#include "fub/solver/DimensionalSplitHyperbolicTimeIntegrator.hpp"

#include "fub/simple_grid/boundary_condition/PeriodicBoundary.hpp"

#include <cmath>
#include <iostream>

using namespace fub;

static double glocke(double x) { return std::exp(-100 * x * x); }

void print(const SimpleGrid<ShallowWater>& grid) {
  auto heigth = grid.states().heigth;
  auto momentum = grid.states().momentum;
  std::cout << "# time point = " << grid.time_point() << '\n';
  for (std::ptrdiff_t i = 0; i < heigth.extent(0); ++i) {
    const double x = 0.5 * grid.dx() + i * grid.dx();
    std::cout << x << "\t" << heigth(i, 0, 0) << '\n';
  }
  std::cout << "\n\n";
}

int main() {
  ShallowWater equation{};

  const std::ptrdiff_t n_cells = 200;
  const int gcw = 1;
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

  while (grid.time_point() < 1.0) {
    double dt = solver.ComputeStableDt(grid, boundary);
    solver.AdvanceGrid(grid, boundary, 0.5 * dt);
    print(grid);
  }
}