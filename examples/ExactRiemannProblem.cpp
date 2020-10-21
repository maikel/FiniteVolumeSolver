#include "fub/Solver.hpp"

int main()
{
  fub::PerfectGas<1> equation;
  fub::ExactRiemannSolver<fub::PerfectGas<1>> solver(equation);

  using Complete = fub::Complete<fub::PerfectGas<1>>;
  
  const double rhoL = 1.0;
  const fub::Array<double, 1, 1> uL{0.0};
  const double pL = 1.0;

  const double rhoR = 0.125;
  const fub::Array<double, 1, 1> uR{0.0};
  const double pR = 0.1;

  Complete left = equation.CompleteFromPrim(rhoL, uL, pL);
  Complete right = equation.CompleteFromPrim(rhoR, uR, pR);
  const auto [pM, uM] = solver.ComputeMiddleState(left, right, fub::Direction::X);
  // std::cout << fmt::format("{:12f}\t{:12f}\n\n", pM, uM);
  Complete solution;

  const int nx = 3200;
  const double x0 = -0.5;
  const double xE = +0.5;
  const double dx = (xE - x0) / nx;
  
  for (int i = 0; i < nx; ++i) {
    const double x = x0 + 0.5 * dx + i * dx;
    solver.SolveRiemannProblem(solution, left, right, pM, uM, fub::Duration(0.2), x, fub::Direction::X);

    const double rho = solution.density;
    const double pressure = solution.pressure;

    std::cout << fmt::format("{:12f}\t{:12f}\n", rho, pressure);
  }
}