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

#ifndef FUB_IDEAL_GAS_FLUX_GODUNOV_METHOD_HPP
#define FUB_IDEAL_GAS_FLUX_GODUNOV_METHOD_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/ForEach.hpp"
#include "fub/core/span.hpp"
#include "fub/flux_method/FluxMethod.hpp"

#include <numeric>

namespace fub {

template <typename EquationT, typename RiemannSolverT,
          int Width = kDefaultChunkSize>
class Godunov {
public:
  using Equation = EquationT;
  using Complete = typename Equation::Complete;
  using Conservative = typename Equation::Conservative;
  using CompleteArray = ::fub::CompleteArray<Equation, Width>;
  using ConservativeArray = ::fub::ConservativeArray<Equation, Width>;

  static constexpr int ChunkSize = Width;

  Godunov(const Equation& equation) : equation_{equation} {}

  void ComputeNumericFlux(Conservative& numeric_flux,
                          span<const Complete, 2> states, Duration /* dt */,
                          double /* dx */, Direction dir) {
    riemann_solver_.SolveRiemannProblem(riemann_solution_, states[0], states[1],
                                        dir);
    equation_.Flux(numeric_flux, riemann_solution_, dir);
  }

  void ComputeNumericFlux(ConservativeArray& numeric_flux,
                          span<const CompleteArray, 2> states,
                          Duration /* dt */, double /* dx */, Direction dir) {
    riemann_solver_.SolveRiemannProblem(riemann_solution_arr_, states[0],
                                        states[1], dir);
    equation_.Flux(numeric_flux, riemann_solution_arr_, dir);
  }

  double ComputeStableDt(span<const Complete, 2> states, double dx,
                         Direction dir) {
    auto signals = riemann_solver_.ComputeSignals(states[0], states[1], dir);
    const double s_max =
        std::accumulate(signals.begin(), signals.end(), 0.0,
                        [](double x, double y) { return std::max(x, y); });
    return dx / s_max;
  }

  static constexpr int GetStencilWidth() noexcept { return 1; }

  const Equation& GetEquation() const noexcept { return equation_; }

private:
  Equation equation_;
  RiemannSolverT riemann_solver_{equation_};
  Complete riemann_solution_{equation_};
  CompleteArray riemann_solution_arr_{equation_};
};

template <typename Equation, typename RPSolver = ExactRiemannSolver<Equation>,
          int Width = kDefaultChunkSize>
struct GodunovMethod : public FluxMethod<Godunov<Equation, RPSolver, Width>> {
  using FluxMethod<Godunov<Equation, RPSolver, Width>>::FluxMethod;
};

template <typename Equation>
GodunovMethod(const Equation& eq)->GodunovMethod<Equation>;

} // namespace fub

#endif