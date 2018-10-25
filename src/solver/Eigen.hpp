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

#include "Eigen/Eigen"
#include "fub/core/span.hpp"

namespace fub {

template <typename T, std::size_t N>
using Array = Eigen::Array<T, N, 1>;

template <typename T, std::ptrdiff_t N>
Eigen::Map<const Eigen::Array<T, N, 1>> asEigen(span<const T, N> array) {
  return Eigen::Map<const Eigen::Array<T, N, 1>>(array.data(), array.size());
}

template <typename T, std::ptrdiff_t N>
Eigen::Map<Eigen::Array<T, N, 1>> asEigen(span<T, N> array) {
  return Eigen::Map<Eigen::Array<T, N, 1>>(array.data(), array.size());
}

template <typename T, int N>
span(Eigen::Array<T, N, 1>& array) -> span<T, N>;

template <typename T, int N>
span(const Eigen::Array<T, N, 1>& array) -> span<const T, N>;

}