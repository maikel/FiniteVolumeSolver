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

#ifndef FUB_NEWTON_ITERATION_HPP
#define FUB_NEWTON_ITERATION_HPP

#include "fub/ext/Eigen.hpp"

namespace fub {

template <typename Function, typename Derivative>
double NewtonIteration(Function&& f, Derivative&& Df, double x0,
                        double tolerance = 1e-7, int max_iterations = 5000) {
  double x = x0;
  double fx = f(x);
  int counter = 0;
  double err = std::abs(fx);
  while (err > tolerance && counter < max_iterations) {
    double dfx = Df(x);
    double offset = fx / dfx;
    double next = x - offset;
    x = ((next - next) != (next - next) || next < 0.0) ? 0.5 * x : next;
    fx = f(x);
    err = std::abs(fx);
    counter += 1;
  }
  if (counter == max_iterations && fx > tolerance) {
    std::ostringstream out;
    out << "Newton Iteration did not converge. (x0 = " << x0 << ", x = " << x
        << ", fx = " << fx << ", counter = " << counter << ", err = " << err
        << ")";
    throw std::runtime_error(out.str());
  }
  return x;
}

} // namespace fub

#endif