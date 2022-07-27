// Copyright (c) 2022 Maikel Nadolski
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

#include "fub/cutcell_method/AnyMdLimiter.hpp"
#include "fub/ext/Eigen.hpp"
#include "fub/core/assert.hpp"
#include "fub/CutCellData.hpp"

#include <limits>

namespace fub
{
namespace {
template <int Rows, int Cols>
Eigen::Matrix<double, 2, Cols>
TakeRows(const Eigen::Matrix<double, Rows, Cols>& A,
         const std::array<int, 2>& W_k) {
  Eigen::Matrix<double, 2, Cols> A_k;
  A_k.row(0) = A.row(W_k[0]);
  A_k.row(1) = A.row(W_k[1]);
  return A_k;
}

bool IsNonSingular(const Eigen::Matrix<double, Eigen::Dynamic, 2>& A, const std::array<int, 2>& W_k)
{
  Eigen::Matrix<double, 2, 2> A_k = TakeRows(A, W_k);
  return A_k.determinant() != 0.0;
}
} // namespace

Eigen::Vector2d SolveLinearOptimizationProblem(
    int maxiter, const Eigen::Vector2d& x_0, const Eigen::Vector2d& c,
    const Eigen::Matrix<double, Eigen::Dynamic, 2>& A,
    const Eigen::VectorXd& b) {
  // Create initial working set
  const int m = A.rows();
  std::array<int, 2> W_k{0, 1};
  Eigen::Matrix<double, 2, 2> A_k = TakeRows(A, W_k);
  Eigen::Vector2d x_k = x_0;

  for (int i = 0; i < maxiter; ++i) {
    // Step 1: Calculate the Lagrange multipliers λ_k in R^2 by solving
    //
    //                      (A_k)T λ_l = c
    //
    FUB_ASSERT(A_k.determinant() != 0.0);
    Eigen::Vector2d lambda_k = A_k.transpose().inverse() * c;

    // Step 2: If λ_k >= 0, STOP. In this case, the point x_k is optimal.
    if ((lambda_k.array() >= 0.0).all()) {
      return x_k;
    }

    // Step 3: Otherwise select q such that (λ_k)q < 0. This constraint will be
    // removed from W_k.
    int q = lambda_k[0] < lambda_k[1] ? 0 : 1;
    int qq = lambda_k[0] < lambda_k[1] ? 1 : 0;

    // Step 4: Calculate the descent direction pk from A_k p_k = e_q, e_q = q-th
    // coordinate vector.
    Eigen::Vector2d p_k = A_k.inverse() * Eigen::Vector2d::Unit(q);

    // Step 5: Choose the step length α_k to be taken along p_k:
    double alpha_k = std::numeric_limits<double>::max();
    int min_i = std::numeric_limits<int>::max();
    for (int i = 0; i < m; ++i) {
      const double ap = A.row(i).transpose().dot(p_k);
      if (ap < 0.0) {
        const double Aix = A.row(i).transpose().dot(x_k);
        const double y_i = (b[i] - Aix) / ap;
        if (y_i < alpha_k && IsNonSingular(A, {W_k[qq], i})) {
          min_i = i;
          alpha_k = y_i;
        }
      }
    }

    // Step 6: Update working set
    x_k += alpha_k * p_k;
    FUB_ASSERT(min_i != W_k[q]);
    W_k[q] = min_i;
    A_k = TakeRows(A, W_k);
  }

  throw std::runtime_error("Convergence Error: Max iterations reached while "
                           "solving a linear optimization probem.");
}

// template <int Rank>
// void BasicHGridReconstruction<Rank>::LimitGradients(
//     const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
//     StridedDataView<const double, Rank> u,
//     StridedDataView<const char, Rank> needs_limiter,
//     const CutCellData<Rank>& geom, const Coordinates<Rank>& dx) const {
//   ForEachIndex(u.Box(), [&](auto... is) {
//     if (needs_limiter(is...) && geom.volume_fractions(is...) > 0.0) {
//       limiter_.LimitGradientsAtIndex(grad_u, u, geom, {is...}, dx);
//     }
//   });
// }

template <int Rank>
bool IsConnected(const CutCellData<Rank>& geom, const Index<Rank>& i,
                 const Index<Rank>& j) {
  if (i != j && geom.volume_fractions(j) > 0.0) {
    Coordinates<Rank> diff;
    for (int d = 0; d < Rank; ++d) {
      diff[d] = static_cast<double>(j[d] - i[d]);
    }
    if (diff.squaredNorm() != 1.0) {
      return false;
    }
    const Coordinates<Rank> xI = GetVolumeCentroid(geom, i);
    const Coordinates<Rank> xJ = diff + GetVolumeCentroid(geom, j);
    const Coordinates<Rank> xB = GetBoundaryCentroid(geom, i);
    const Coordinates<Rank> n = GetBoundaryNormal(geom, i);
    return n.dot(xJ - xI) >= n.dot(xB);
  } else {
    return false;
  }
}


template <int Rank>
void NoMdLimiter<Rank>::LimitGradientsAtIndex(
    const std::array<StridedDataView<double, Rank>, Rank>&,
    StridedDataView<const double, Rank>, const CutCellData<Rank>&,
    const Index<Rank>&, const Coordinates<Rank>&) const {}

template struct NoMdLimiter<2>;

template <int Rank>
void UpwindMdLimiter<Rank>::LimitGradientsAtIndex(
    const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
    StridedDataView<const double, Rank>, const CutCellData<Rank>&,
    const Index<Rank>& index, const Coordinates<Rank>&) const {
  grad_u[0](index) = 0.0;
  grad_u[1](index) = 0.0;
}

template struct UpwindMdLimiter<2>;

template <int Rank>
void LinearOptimizationLimiter<Rank>::LimitGradientsAtIndex(
    const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
    StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
    const Index<Rank>& index, const Coordinates<Rank>& dx) const {
  if constexpr (Rank == 2) {
    IndexBox<Rank> neighborhood = Neighborhood<Rank>(index, 1);
    double u_M = u(index);
    double Dx = grad_u[0](index);
    double Dy = grad_u[1](index);
    Eigen::Vector2d D(Dx, Dy);
    const Coordinates<Rank> x_M = GetAbsoluteVolumeCentroid(geom, index, dx);
    // we want to minimze cost_vector.dot(limiter)
    Eigen::Vector2d cost_vector = -D.array().abs();

    // to build up positivity constraints we count how many neighbors we have
    std::array<Index<Rank>, 9> neighbors;
    int count_neighbors = 0;
    ForEachIndex(neighborhood, [&](int i, int j) {
      if (IsConnected(geom, index, {i, j})) {
        neighbors[count_neighbors] = Index<Rank>{i, j};
        count_neighbors += 1;
      }
    });

    const int n = 4 + 2 * count_neighbors;
    Eigen::Matrix<double, Eigen::Dynamic, 2> A(n, 2);
    Eigen::Matrix<double, Eigen::Dynamic, 1> b(n);

    // The first 4 constraints are

    // 0 <= limiter_x
    A.row(0) = Eigen::Vector2d(1, 0);
    b[0] = 0;
    // 0 <= limiter_y
    A.row(1) = Eigen::Vector2d(0, 1);
    b[1] = 0;

    // limiter_x <= 1   <=>   -1 <= -limiter_x
    A.row(2) = Eigen::Vector2d(-1, 0);
    b[2] = -1;
    // limiter_y <= 1   <=>   -1 <= -limiter_y
    A.row(3) = Eigen::Vector2d(0, -1);
    b[3] = -1;

    int r = 4;
    for (const Index<Rank>& i : neighbors) {
      if (r == n) {
        break;
      }
      double u_i = u(i);
      const Coordinates<Rank> x_i = GetAbsoluteVolumeCentroid(geom, i, dx);
      Eigen::Vector2d a = (x_i - x_M).array() * D.array();
      if (u_i < u_M) {
        A.row(r) = a;
        b[r] = u_i - u_M;
        A.row(r + 1) = -a;
        b[r + 1] = 0.0;
      } else {
        A.row(r) = -a;
        b[r] = -(u_i - u_M);
        A.row(r + 1) = a;
        b[r + 1] = 0.0;
      }
      r += 2;
    }

    const int max_iter = 100;
    Eigen::Vector2d limiter_initial_guess = Eigen::Vector2d::Zero();
    Eigen::Vector2d limiter = SolveLinearOptimizationProblem(
        max_iter, limiter_initial_guess, cost_vector, A, b);

    // Now we do the limiting!
    grad_u[0](index) = limiter[0] * D[0];
    grad_u[1](index) = limiter[1] * D[1];
  }
}

}