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

template <typename EquationT,
          typename RiemannSolverT = ExactRiemannSolver<EquationT>>
class Godunov {
public:
  using Equation = EquationT;
  using RiemannSolver = RiemannSolverT;
  using Complete = typename Equation::Complete;
  using Conservative = typename Equation::Conservative;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;

  Godunov(const Equation& equation) : equation_{equation} {}
  Godunov(const Equation& equation, const RiemannSolver& solver)
      : equation_{equation}, riemann_solver_{solver} {}

  static constexpr int GetStencilWidth() noexcept;

  const Equation& GetEquation() const noexcept;
  Equation& GetEquation() noexcept;

  const RiemannSolver& GetRiemannSolver() const noexcept;
  RiemannSolver& GetRiemannSolver() noexcept;

  void ComputeNumericFlux(Conservative& numeric_flux,
                          span<const Complete, 2> states, Duration /* dt */,
                          double /* dx */, Direction dir);

  void ComputeNumericFlux(ConservativeArray& numeric_flux,
                          span<const CompleteArray, 2> states,
                          Duration /* dt */, double /* dx */, Direction dir);

  void ComputeNumericFlux(ConservativeArray& numeric_flux,
                          Array1d face_fraction,
                          span<const CompleteArray, 2> states,
                          span<const Array1d, 2> volume_fraction,
                          Duration /* dt */, double /* dx */, Direction dir);

  double ComputeStableDt(span<const Complete, 2> states, double dx,
                         Direction dir);

  Array1d ComputeStableDt(span<const CompleteArray, 2> states, double dx,
                          Direction dir);

  Array1d ComputeStableDt(span<const CompleteArray, 2> states,
                          Array1d face_fraction,
                          span<const Array1d, 2> volume_fraction, double dx,
                          Direction dir);

private:
  Equation equation_;
  RiemannSolver riemann_solver_{equation_};
  Complete riemann_solution_{equation_};
  CompleteArray riemann_solution_arr_{equation_};
};

/// \ingroup flux-method
template <typename Equation, typename RPSolver = ExactRiemannSolver<Equation>>
struct GodunovMethod : public FluxMethod<Godunov<Equation, RPSolver>> {
  using FluxMethod<Godunov<Equation, RPSolver>>::FluxMethod;
};

template <typename Equation>
GodunovMethod(const Equation& eq)->GodunovMethod<Equation>;

///////////////////////////////////////////////////////////////////////////////
// Implementation

template <typename Equation, typename RiemannSolver>
constexpr int Godunov<Equation, RiemannSolver>::GetStencilWidth() noexcept {
  return 1;
}

template <typename Equation, typename RiemannSolver>
const Equation& Godunov<Equation, RiemannSolver>::GetEquation() const noexcept {
  return equation_;
}

template <typename Equation, typename RiemannSolver>
Equation& Godunov<Equation, RiemannSolver>::GetEquation() noexcept {
  return equation_;
}

template <typename Equation, typename RiemannSolver>
const RiemannSolver& Godunov<Equation, RiemannSolver>::GetRiemannSolver() const
    noexcept {
  return riemann_solver_;
}

template <typename Equation, typename RiemannSolver>
RiemannSolver& Godunov<Equation, RiemannSolver>::GetRiemannSolver() noexcept {
  return riemann_solver_;
}

template <typename Equation, typename RiemannSolver>
void Godunov<Equation, RiemannSolver>::ComputeNumericFlux(
    Conservative& numeric_flux, span<const Complete, 2> states,
    Duration /* dt */, double /* dx */, Direction dir) {
  riemann_solver_.SolveRiemannProblem(riemann_solution_, states[0], states[1],
                                      dir);
  Flux(equation_, numeric_flux, riemann_solution_, dir);
}

template <typename Equation, typename RiemannSolver>
void Godunov<Equation, RiemannSolver>::ComputeNumericFlux(
    ConservativeArray& numeric_flux, span<const CompleteArray, 2> states,
    Duration /* dt */, double /* dx */, Direction dir) {
  riemann_solver_.SolveRiemannProblem(riemann_solution_arr_, states[0],
                                      states[1], dir);
  Flux(equation_, numeric_flux, riemann_solution_arr_, dir);
}

template <typename Equation, typename RiemannSolver>
void Godunov<Equation, RiemannSolver>::ComputeNumericFlux(
    ConservativeArray& numeric_flux, Array1d face_fraction,
    span<const CompleteArray, 2> states,
    span<const Array1d, 2> /* volume_fraction */, Duration /* dt */,
    double /* dx */, Direction dir) {
  MaskArray mask = face_fraction > 0.0;
  riemann_solver_.SolveRiemannProblem(riemann_solution_arr_, states[0],
                                      states[1], mask, dir);
  Flux(equation_, numeric_flux, riemann_solution_arr_, dir);
  const Array1d zero = Array1d::Zero();
  ForEachComponent([&](auto&& nf, Array1d f) { nf = mask.select(f, zero); },
                   numeric_flux, numeric_flux);
}

template <typename Equation, typename RiemannSolver>
double Godunov<Equation, RiemannSolver>::ComputeStableDt(
    span<const Complete, 2> states, double dx, Direction dir) {
  auto signals = riemann_solver_.ComputeSignals(states[0], states[1], dir);
  const double s_max = std::accumulate(
      signals.begin(), signals.end(), 0.0,
      [](double x, double y) { return std::max(x, std::abs(y)); });
  return dx / s_max;
}

template <typename Equation, typename RiemannSolver>
Array1d Godunov<Equation, RiemannSolver>::ComputeStableDt(
    span<const CompleteArray, 2> states, double dx, Direction dir) {
  auto signals = riemann_solver_.ComputeSignals(states[0], states[1], dir);
  Array1d zero = Array1d::Zero();
  const Array1d s_max =
      std::accumulate(signals.begin(), signals.end(), zero,
                      [](const Array1d& x, const Array1d& y) -> Array1d {
                        return x.max(y.abs());
                      });
  return Array1d(dx) / s_max;
}

template <typename Equation, typename RiemannSolver>
Array1d Godunov<Equation, RiemannSolver>::ComputeStableDt(
    span<const CompleteArray, 2> states, Array1d face_fraction,
    span<const Array1d, 2>, double dx, Direction dir) {
  std::array<Complete, 2> state{};
  Array1d dts = Array1d::Constant(std::numeric_limits<double>::infinity());
  for (int i = 0; i < face_fraction.size(); ++i) {
    if (face_fraction[i] > 0.0) {
      Load(state[0], states[0], i);
      Load(state[1], states[1], i);
      dts[i] = ComputeStableDt(state, dx, dir);
    }
  }
  return dts;
}

} // namespace fub

#endif
