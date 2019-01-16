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

#include "fub/geometry/Halfspace.hpp"

#include "fub/core/assert.hpp"

namespace fub {
Halfspace::Halfspace(const Coordinates& normal, double offset)
    : normal_{normal}, offset_{offset} {
  double norm2 = 0.0;
  for (double n : normal_) {
    norm2 += n * n;
  }
  FUB_ASSERT(norm2 > 0.0);
  const double norm = std::sqrt(norm2);
  for (double& n : normal_) {
    n /= norm;
  }
}

double Halfspace::ComputeDistanceTo(const Coordinates& x) const {
  FUB_ASSERT(x.size() == normal_.size());
  const int size = x.size();
  double proj = -offset_;
  for (int i = 0; i < size; ++i) {
    proj += x[i] * normal_[i];
  }
  return proj;
}
} // namespace fub