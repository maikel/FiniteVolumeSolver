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

#include <Eigen/Eigen>

#include "fub/PatchDataView.hpp"
#include "fub/Direction.hpp"

#include <type_traits>

namespace fub {

constexpr const int kDefaultChunkSize = 8;

template <typename T, int N, int M = kDefaultChunkSize>
using Array = std::conditional_t<M == 1, Eigen::Array<T, 1, N>,
                                 Eigen::Array<T, M, N, Eigen::ColMajor>>;

using Array1d = Array<double, 1>;
// using Scalar = Array1d;
using Array2d = Array<double, 2>;
using Array3d = Array<double, 3>;
using ArrayXd = Array<double, Eigen::Dynamic>;

template <int N, typename T, typename Extents, typename Layout,
          typename... Indices,
          typename = std::enable_if_t<(std::is_integral_v<Indices> && ...)>>
Eigen::Array<std::remove_cv_t<T>, N, 1>
LoadN(constant<N>, basic_mdspan<T, Extents, Layout> mdspan, int size,
      std::ptrdiff_t i0, Indices... indices) {
  FUB_ASSERT(0 < size && size < N);
  Eigen::Array<std::remove_cv_t<T>, N, 1> array;
  for (int i = 0; i < size; ++i) {
    array[i] = mdspan(i0 + i, indices...);
  }
  return array;
}

template <std::size_t N>
std::array<std::ptrdiff_t, N> Shift(const std::array<std::ptrdiff_t, N>& idx,
                                    Direction dir, std::ptrdiff_t shift) {
  auto shifted(idx);
  shifted[static_cast<std::size_t>(dir)] += shift;
  return shifted;
}

template <int Rank>
Eigen::Matrix<double, Rank, 1> UnitVector(Direction dir) noexcept {
  Eigen::Matrix<double, Rank, 1> unit = Eigen::Matrix<double, Rank, 1>::Zero();
  Eigen::Index dir_v = static_cast<Eigen::Index>(dir);
  unit[dir_v] = 1.0;
  return unit;
}

inline Eigen::Matrix<double, 2, 2> MakeRotation(const Eigen::Matrix<double, 2, 1>& a,
                                         const Eigen::Matrix<double, 2, 1>& b) {
  if (a == -b) {
    return -Eigen::Matrix<double, 2, 2>::Identity();
  }
  const double s = a[0]*b[1] - a[1]*b[0];
  const double c = a.dot(b);
  Eigen::Matrix<double, 2, 2> rotation;
  rotation << c, -s, s, c;
  return rotation;
}

template <int N, typename T, typename Extents, typename Layout,
          typename... Indices>
Eigen::Array<std::remove_cv_t<T>, N, 1>
LoadN(constant<N>, basic_mdspan<T, Extents, Layout> mdspan, int size,
      const std::array<std::ptrdiff_t, Extents::rank()>& index) {
  FUB_ASSERT(0 < size && size < N);
  Eigen::Array<std::remove_cv_t<T>, N, 1> array{};
  for (int i = 0; i < size; ++i) {
    array[i] = mdspan(Shift(index, Direction::X, i));
  }
  return array;
}

template <int N, typename T, typename Extents, typename Layout,
          typename... Indices>
Eigen::Array<std::remove_cv_t<T>, N, 1>
Load(constant<N> n, basic_mdspan<T, Extents, Layout> mdspan,
     Indices... indices) {
  if (mdspan.stride(0) == 1) {
    return Eigen::Map<const Eigen::Array<std::remove_cv_t<T>, N, 1>>(
        &mdspan(indices...), n);
  } else {
    using Stride = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
    return Eigen::Map<const Eigen::Array<std::remove_cv_t<T>, N, 1>,
                      Eigen::Unaligned, Stride>(&mdspan(indices...), n, 1,
                                                Stride(1, mdspan.stride(0)));
  }
}

template <typename T, typename Extents, typename Layout, typename... Indices>
auto Load(basic_mdspan<T, Extents, Layout> mdspan, Indices... indices) {
  return Load(constant<kDefaultChunkSize>(), mdspan, indices...);
}

template <typename T, int N, typename Layout, typename... Indices>
void StoreN(mdspan<T, 1, Layout> mdspan, int n,
            const Eigen::Array<T, N, 1>& chunk, std::ptrdiff_t offset) {
  for (int i = 0; i < n; ++i) {
    mdspan(offset + i) = chunk[i];
  }
}

template <typename T, int Rank, typename Layout, int N, int Options,
          typename... Indices>
void Store(const PatchDataView<T, Rank, Layout>& view,
           const Eigen::Array<T, N, 1, Options>& chunk,
           const nodeduce_t<std::array<std::ptrdiff_t, Rank>>& index) {
  FUB_ASSERT(view.Extent(0) >= N);
  if (view.stride(0) == 1) {
    Eigen::Map<Eigen::Array<T, N, 1>> mapped(&view(index), N);
    mapped = chunk;
  } else {
    using Stride = Eigen::Stride<1, Eigen::Dynamic>;
    Eigen::Map<Eigen::Array<T, N, 1>, Eigen::Unaligned, Stride> mapped(
        &view(index), N, Stride(1, view.stride(0)));
    mapped = chunk;
  }
  return;
}

template <typename T, int Rank, typename L>
void Store(const PatchDataView<T, Rank, L>& view, nodeduce_t<const T&> value,
           const nodeduce_t<std::array<std::ptrdiff_t, Rank>>& index) {
  view(index) = value;
}

} // namespace fub

#endif