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

#include "fub/equations/PerfectGas.hpp"
#include <iostream>

int main() {
  using namespace fub;
  PerfectGas<1> equation{};
  ExactRiemannSolver<PerfectGas<1>> solver{equation};

  using Complete = PerfectGas<1>::Complete;
  using Conservative = PerfectGas<1>::Conservative;
  auto from_prim = [&equation](Complete& state) {
    state.energy = state.pressure * equation.gamma_minus_1_inv;
    state.speed_of_sound =
        std::sqrt(equation.gamma * state.pressure / state.density);
  };

  Complete left;
  Complete right;
  Complete solution;

  // (i) SOD problem.

  left.density = 1.0;
  left.momentum = 0.0;
  left.pressure = 1.0;
  from_prim(left);

  right.density = 0.125;
  right.momentum = 0.0;
  right.pressure = 0.1;
  from_prim(right);

  double pM;
  double uM;
  std::tie(pM, uM) = solver.ComputeMiddleState(left, right, Direction::X);
  std::cout << "p*: " << pM << ", u*: " << uM << '\n';

  solver.SolveRiemannProblem(solution, left, right, Direction::X);
  std::cout << "solution: rho=" << solution.density
            << ", u=" << solution.momentum / solution.density
            << ", p=" << solution.pressure << "\n\n";

  // (ii) 123 problem.

  left.density = 1.0;
  left.momentum = -2.0;
  left.pressure = 0.4;
  from_prim(left);

  right.density = 1.0;
  right.momentum = +2.0;
  right.pressure = 0.4;
  from_prim(right);

  std::tie(pM, uM) = solver.ComputeMiddleState(left, right, Direction::X);
  std::cout << "p*: " << pM << ", u*: " << uM << '\n';

  solver.SolveRiemannProblem(solution, left, right, Direction::X);
  std::cout << "solution: rho=" << solution.density
            << ", u=" << solution.momentum / solution.density
            << ", p=" << solution.pressure << "\n\n";

  // (iii) severe SOD problem

  left.density = 1.0;
  left.momentum = 0.0;
  left.pressure = 1000.0;
  from_prim(left);

  right.density = 1.0;
  right.momentum = 0.0;
  right.pressure = 0.01;
  from_prim(right);

  std::tie(pM, uM) = solver.ComputeMiddleState(left, right, Direction::X);
  std::cout << "p*: " << pM << ", u*: " << uM << '\n';

  solver.SolveRiemannProblem(solution, left, right, Direction::X);
  std::cout << "solution: rho=" << solution.density
            << ", u=" << solution.momentum / solution.density
            << ", p=" << solution.pressure << "\n\n";

  // (iv) Woodward and Colella

  left.density = 1.0;
  left.momentum = 0.0;
  left.pressure = 0.1;
  from_prim(left);

  right.density = 1.0;
  right.momentum = 0.0;
  right.pressure = 100.0;
  from_prim(right);

  std::tie(pM, uM) = solver.ComputeMiddleState(left, right, Direction::X);
  std::cout << "p*: " << pM << ", u*: " << uM << '\n';

  solver.SolveRiemannProblem(solution, left, right, Direction::X);
  std::cout << "solution: rho=" << solution.density
            << ", u=" << solution.momentum / solution.density
            << ", p=" << solution.pressure << "\n\n";

  left.density =  2.52667;
  left.momentum = 1204.31;
  left.pressure = 312160;
  from_prim(left);
  
  right.density = 2.28972;
  right.momentum = 1150.82;
  right.pressure = 278308;
  from_prim(right);

  std::tie(pM, uM) = solver.ComputeMiddleState(left, right, Direction::X);
  std::cout << "p*: " << pM << ", u*: " << uM << '\n';

  solver.SolveRiemannProblem(solution, left, right, Direction::X);
  std::cout << "solution: rho=" << solution.density
            << ", u=" << solution.momentum / solution.density
            << ", p=" << solution.pressure << "\n\n";

  Conservative f_left;
  equation.Flux(f_left, solution, Direction::X);

  // Linear Shock Tube problematic states
  left.density = 2.28972;
  left.momentum = 1150.82;
  left.pressure = 278308;
  from_prim(left);

  right.density = 2.06243;
  right.momentum = 1084.96;
  right.pressure = 248331;
  from_prim(right);

  std::tie(pM, uM) = solver.ComputeMiddleState(left, right, Direction::X);
  std::cout << "p*: " << pM << ", u*: " << uM << '\n';

  solver.SolveRiemannProblem(solution, left, right, Direction::X);
  std::cout << "solution: rho=" << solution.density
            << ", u=" << solution.momentum / solution.density
            << ", p=" << solution.pressure << "\n\n";

  Conservative f_right;
  equation.Flux(f_right, solution, Direction::X);

  std::cout << f_left.density << " - " << f_right.density << '\n';
  std::cout << f_left.momentum << " - " << f_right.momentum << '\n';
  std::cout << f_left.energy << " - " << f_right.energy << '\n';
}