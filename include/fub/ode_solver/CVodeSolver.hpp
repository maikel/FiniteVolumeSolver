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

#ifndef FUB_ODE_SOLVER_CVODE_SOLVER_HPP
#define FUB_ODE_SOLVER_CVODE_SOLVER_HPP

#include "fub/config.hpp"
#include "fub/ode_solver/OdeSolverFactory.hpp"
#include <functional>
#include <vector>

namespace fub {

struct CVodeSolverOptions : public OdeSolverOptions {
  explicit CVodeSolverOptions(int nvars)
      : OdeSolverOptions{OdeSolverType::CVode}, n_variables{nvars} {}

  int n_variables;
};

class CVodeSolver : public OdeSolver {
public:
  using root_function_type =
      std::function<int(span<double>, span<const double>, double)>;

  CVodeSolver(int size);

  void SetRootFunction(root_function_type root_function, int size) {
    root_function_ = std::move(root_function);
    root_size_ = size;
  }

private:
  std::unique_ptr<OdeSolver> clone() const override;

  void integrate(system_type system, span<double> y_0, double t, double dt,
                 feedback_type* feedback,
                 jacobian_type* jacobian) const override;

  struct CVodeInternalData_;
  struct ReleaseCVodeMemory {
    void operator()(CVodeInternalData_* p) const noexcept;
  };

  std::unique_ptr<CVodeInternalData_, ReleaseCVodeMemory> data_;
  root_function_type root_function_{};
  int root_size_{};
  double relative_tolerance_{1e-10};
};

#ifdef FUB_WITH_SUNDIALS
static RegisterSpecificFactory register_cvode_solver_{
    "CVode", [](const OdeSolverOptions& base) {
      FUB_ASSERT(base.type == OdeSolverType::CVode);
      const CVodeSolverOptions* opts =
          static_cast<const CVodeSolverOptions*>(&base);
      return std::make_unique<CVodeSolver>(opts->n_variables);
    }};
#endif

} // namespace fub

#endif