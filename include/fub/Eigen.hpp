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

#ifndef FUB_CORE_EIGEN_HPP
#define FUB_CORE_EIGEN_HPP

#include <unsupported/Eigen/CXX11/Tensor>

#include <fub/core/mdspan.hpp>

#include <type_traits>

namespace fub {

constexpr const int kChunkSize = 8;

template <typename T, int N>
using Array =
    std::conditional_t<kChunkSize == 1, Eigen::Array<T, 1, N>,
                       Eigen::Array<T, kChunkSize, N, Eigen::ColMajor>>;

using Array1d = Array<double, 1>;
using Scalar = Array1d;
using Array2d = Array<double, 2>;
using Array3d = Array<double, 3>;
using ArrayXd = Array<double, Eigen::Dynamic>;

template <typename T, int Rank>
using TensorSpan = DynamicMdSpan<T, Rank, layout_left>;

template <typename T, typename Extents, typename Layout, typename... Indices>
Array<std::remove_cv_t<T>, 1> Load(basic_mdspan<T, Extents, Layout> mdspan,
                                   Indices... indices) {
  FUB_ASSERT(mdspan.stride(0) == 1);
  std::array<std::ptrdiff_t, sizeof...(Indices)> is{
      static_cast<std::ptrdiff_t>(indices)...};
  if (mdspan.extent(0) - is[0] < kChunkSize) {
    return LoadN(mdspan, mdspan.extent(0) - is[0], indices...);
  }
  Eigen::Map<const Array<std::remove_cv_t<T>, 1>> mapped(&mdspan(indices...),
                                                         kChunkSize);
  return mapped;
}

template <typename T, typename Extents, typename Layout, typename... Indices>
Array<std::remove_cv_t<T>, 1> LoadN(basic_mdspan<T, Extents, Layout> mdspan,
                                    int n, Indices... indices) {
  FUB_ASSERT(0 < n && n < kChunkSize);
  Array<std::remove_cv_t<T>, 1> array;
  auto out = std::copy_n(&mdspan(indices...), n, array.data());
  std::fill_n(out, kChunkSize - n, array[0]);
  return array;
}

template <typename T, typename... Indices>
void Store(basic_mdspan<T, extents<dynamic_extent>, layout_stride> mdspan,
           const Array<std::remove_cv_t<T>, 1>& chunk, Indices... offsets) {
  std::array<std::ptrdiff_t, sizeof...(Indices)> is{
      static_cast<std::ptrdiff_t>(offsets)...};
  if (mdspan.extent(0) - is[0] < kChunkSize) {
    StoreN(mdspan, mdspan.extent(0) - is[0], chunk, offsets...);
    return;
  }
  Eigen::Map<Array<std::remove_cv_t<T>, 1>> mapped(&mdspan(offsets...),
                                                   kChunkSize);
  mapped = chunk;
}

template <typename F, typename First, typename... Arrays>
F ForEachCol(F function, First&& first, Arrays&&... arrays) {
  const int cols = first.cols();
  for (int c = 0; c < cols; ++c) {
    function(first.col(c), arrays.col(c)...);
  }
  return function;
}

template <typename T, typename... Indices>
void StoreN(basic_mdspan<T, extents<dynamic_extent>, layout_stride> mdspan,
            int limited, const Array<std::remove_cv_t<T>, 1>& chunk,
            Indices... offsets) {
  std::copy_n(chunk.data(), limited, &mdspan(offsets...));
}

} // namespace fub

#endif