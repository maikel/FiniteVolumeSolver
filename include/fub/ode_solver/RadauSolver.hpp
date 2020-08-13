// Copyright (c) 2018 Maikel Nadolski
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

#ifndef FUB_ODE_SOLVER_RADAU_SOLVER_HPP
#define FUB_ODE_SOLVER_RADAU_SOLVER_HPP

#include "fub/ode_solver/OdeSolverFactory.hpp"
#include <stdexcept>
#include <vector>

namespace fub {
namespace radau_error {
struct Exception : std::runtime_error {
  using runtime_error::runtime_error;
};

struct InputIsNotConsistent : Exception {
  InputIsNotConsistent() : Exception("RadauSolver: Input is not consistent.") {}
};

struct LargerNMaxIsNeeded : Exception {
  LargerNMaxIsNeeded() : Exception("RadauSolver: Larger NMAX is needed.") {}
};

struct StepSizeBecomesTooSmall : Exception {
  StepSizeBecomesTooSmall()
      : Exception("RadauSolver: Step size becomes too small.") {}
};

struct MatrixIsRepeatedlySingular : Exception {
  MatrixIsRepeatedlySingular()
      : Exception("RadauSolver: Matrix is repeatedly singular.") {}
};

struct RhsEvaluatesToNaN : Exception {
  RhsEvaluatesToNaN() : Exception("RadauSolver: RHS evaluates to NaN.") {}
};
} // namespace radau_error

class RadauSolver : public OdeSolver {
public:
  double relative_tolerance{1e-12};
  double absolute_tolerance{1e-13};

  double initial_time_step_size{1e-6};

  // Recompute Jacobian less often (default 0.001)
  double jacobian_threshhold{0.01};

  RadauSolver() = default;

  static std::ptrdiff_t GetWorkSpaceSize(int size);
  static std::ptrdiff_t GetIWorkSpaceSize(int size);

private:
  std::unique_ptr<OdeSolver> clone() const override {
    return std::make_unique<RadauSolver>();
  }

  void integrate(system_type system, span<double> y_0, double t, double dt,
                 feedback_type* feedback,
                 jacobian_type* jacobian) const override;

  mutable std::vector<double> work_space_{};
  mutable std::vector<int> iwork_space_{};
};

static RegisterSpecificFactory register_radau_solver_{
    "RADAU",
    [](const OdeSolverOptions&) { return std::make_unique<RadauSolver>(); }};

} // namespace fub

#endif
