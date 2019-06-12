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

#ifndef FUB_EXECUTION_HPP
#define FUB_EXECUTION_HPP

#include "fub/ext/omp.hpp"
#include <optional>

namespace fub {
namespace execution {

struct SequentialTag {};
inline constexpr SequentialTag seq{};

struct SimdTag {};
inline constexpr SimdTag simd{};

struct OpenMpTag {};
inline constexpr OpenMpTag openmp{};

struct OpenMpSimdTag : OpenMpTag, SimdTag {};
inline constexpr OpenMpSimdTag openmp_simd{};

} // namespace execution

namespace detail {
template <typename Tag, typename T> struct LocalType {
  using type = std::optional<T>;
};
template <typename T> struct LocalType<execution::OpenMpTag, T> {
  using type = OmpLocal<T>;
};
template <typename T> struct LocalType<execution::OpenMpSimdTag, T> {
  using type = OmpLocal<T>;
};
} // namespace detail
template <typename Tag, typename T>
using Local = typename detail::LocalType<Tag, T>::type;

} // namespace fub

#endif
