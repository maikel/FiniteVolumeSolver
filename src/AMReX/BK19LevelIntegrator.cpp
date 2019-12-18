// Copyright (c) 2019 Maikel Nadolski
// Copyright (c) 2019 Stefan Vater
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

#include "fub/AMReX/BK19LevelIntegrator.hpp"

namespace fub::amrex {
static constexpr int Rank = AMREX_SPACEDIM;
using Equation = CompressibleAdvection<Rank>;

namespace {
double EvaluateRhsU_(const Equation& equation, const Complete& state, Duration time_step_size) {
  const double f = equation.f;
  const double dt = time_step_size.count();
  const double U = state.velocity[0];
  const double V = state.velocity[1];
  const double dt_times_f = dt * f;
  const double dt_times_f_square = dt_times_f * dt_times_f;
  return (U + dt * f * V) / (1 + dt_times_f_square);
}

double EvaluateRhsV_(const Equation& equation, const Complete& state, Duration time_step_size) {
  const double f = equation.f(x);
  const double dt = time_step_size.count();
  const double U = state.velocity[0];
  const double V = state.velocity[1];
  const double dt_times_f = dt * f;
  const double dt_times_f_square = dt_times_f * dt_times_f;
  return (V - dt * f * U) / (1 + dt_times_f_square);
}

double EvaluateRhsW_(const Equation& equation, const Complete& state, Duration time_step_size) {
  const double dt = time_step_size.count();
  const double U = state.velocity[0];
  const double V = state.velocity[1];
  const double dt_times_f = dt * f;
  const double dt_times_f_square = dt_times_f * dt_times_f;
  return (V - dt * f * U) / (1 + dt_times_f_square);
}

}
}