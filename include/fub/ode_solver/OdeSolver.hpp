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

#ifndef FUB_ODE_SOLVER_ODE_SOLVER_HPP
#define FUB_ODE_SOLVER_ODE_SOLVER_HPP

#include "fub/core/function_ref.hpp"
#include "fub/core/span.hpp"

#include <memory>

namespace fub {

struct OdeSolver {
  using system_type =
      function_ref<void(span<double>, span<const double>, double)>;

  using feedback_type = function_ref<void(span<const double>, double)>;

  using jacobian_type =
      function_ref<void(span<double>, span<const double>, double)>;

  virtual ~OdeSolver() = default;

  void Integrate(system_type rhs, span<double> y_0, double t, double dt) const {
    integrate(rhs, y_0, t, dt, nullptr, nullptr);
  }

  void Integrate(system_type rhs, span<double> y_0, double t, double dt,
                 feedback_type feedback) const {
    integrate(rhs, y_0, t, dt, &feedback, nullptr);
  }

  void Integrate(system_type rhs, span<double> y_0, double t, double dt,
                 jacobian_type jacobian) const {
    integrate(rhs, y_0, t, dt, nullptr, &jacobian);
  }

  void Integrate(system_type rhs, span<double> y_0, double t, double dt,
                 feedback_type feedback, jacobian_type jacobian) const {
    integrate(rhs, y_0, t, dt, &feedback, &jacobian);
  }

  std::unique_ptr<OdeSolver> Clone() const { return clone(); }

private:
  virtual void integrate(system_type system, span<double> y_0, double t,
                         double dt, feedback_type* feedback,
                         jacobian_type* jacobian) const = 0;

  virtual std::unique_ptr<OdeSolver> clone() const = 0;
};

} // namespace fub

#endif