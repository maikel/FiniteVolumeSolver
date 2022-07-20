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

#ifndef FUB_CUTCELL_METHOD_MY_STABILISATION_HPP
#define FUB_CUTCELL_METHOD_MY_STABILISATION_HPP

#include "fub/CutCellData.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"
#include "fub/StateRow.hpp"
#include "fub/core/span.hpp"
#include "fub/core/tuple.hpp"
#include "fub/equations/RequireMassflow.hpp"

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>

#include <algorithm>

namespace fub {

template <int Rank>
IndexBox<Rank> Neighborhood(const Index<Rank>& i, int width) {
  Index<Rank> lower = i;
  Index<Rank> upper = i;
  for (int i = 0; i < Rank; ++i) {
    lower[i] -= width;
    upper[i] += width + 1;
  }
  return {lower, upper};
}

void MyStab_ComputeStableFluxComponents(
    const PatchDataView<double, 2, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, 2, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, 2, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, 2, layout_stride>& regular_fluxes,
    const PatchDataView<const double, 2, layout_stride>& boundary_fluxes,
    const CutCellData<2>& geom, Duration dt, double dx, Direction dir);

void MyStab_ComputeStableFluxComponents(
    const PatchDataView<double, 3, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, 3, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, 3, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, 3, layout_stride>& regular_fluxes,
    const PatchDataView<const double, 3, layout_stride>& boundary_fluxes,
    const CutCellData<3>& geom, Duration dt, double dx, Direction dir);

template <int Rank> using Coordinates = Eigen::Matrix<double, Rank, 1>;

inline Coordinates<2>
ComputeReflectedCoordinates(const Coordinates<2>& offset,
                            const Coordinates<2>& boundary_normal) {
  Coordinates<2> reflected =
      offset - 2.0 * offset.dot(boundary_normal) * boundary_normal;
  return reflected;
}

template <int Rank>
Index<Rank> RelativeCellIndex(const Coordinates<Rank>& x,
                              const Coordinates<Rank>& dx) {
  Index<Rank> i{};
  Coordinates<Rank> rel_x = x / dx;
  for (int d = 0; d < Rank; ++d) {
    i[d] = std::floor(rel_x[d]);
  }
  return i;
}

template <typename State, std::ptrdiff_t Rank>
void ApplyGradient(
    State& u, span<const State, Rank> grad,
    nodeduce_t<const Eigen::Matrix<double, Rank, 1>&> x) noexcept {
  if constexpr (Rank == 2) {
    ForEachComponent(
        [&x](double& u, auto&&... grad_x) {
          Coordinates<Rank> grad{grad_x...};
          u = grad.dot(x);
        },
        u, grad[0], grad[1]);
  } else if constexpr (Rank == 3) {
    ForEachComponent(
        [&x](double& u, auto&&... grad_x) {
          Coordinates<Rank> grad{grad_x...};
          u = grad.dot(x);
        },
        u, grad[0], grad[1], grad[2]);
  }
}

template <int Rank> struct BasicHGridReconstruction {
  void ComputeGradients(span<double, 2> gradient, span<const double, 4> states,
                        span<const Coordinates<Rank>, 4> x);

  void ComputeGradients(span<double, 2> gradient, span<const double, 5> states,
                        span<const Coordinates<Rank>, 5> x);

  void
  LimitGradients(const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
                 StridedDataView<const double, Rank> u,
                 StridedDataView<const char, Rank> needs_limiter,
                 const CutCellData<Rank>& geom,
                 const Coordinates<Rank>& dx) const;

private:
  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const;
};

extern template struct BasicHGridReconstruction<2>;

template <typename Equation>
struct ConservativeHGridReconstruction
    : BasicHGridReconstruction<Equation::Rank()> {
  using Base = BasicHGridReconstruction<Equation::Rank()>;

  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;

  using K = CGAL::Exact_predicates_exact_constructions_kernel;
  using Point_2 = K::Point_2;
  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  static constexpr int Rank = Equation::Rank();

  struct IntegrationPoints {
    static constexpr int kMaxSources = 9;

    std::array<Index<Rank>, kMaxSources> index{};
    std::array<double, kMaxSources> volume{};
    std::array<Coordinates<Rank>, kMaxSources> xM{};

    Index<Rank> iB{};
    Coordinates<Rank> xB{};
  };

  Point_2 ShiftFrom(const Coordinates<2>& x, double xshift, double yshift) {
    Point_2 p{x[0] + xshift, x[1] + yshift};
    return p;
  }

  Polygon_2 GetPolygon(const CutCellData<2>& geom, const Index<2>& index,
                       const Coordinates<2>& dx) {
    std::vector<Point_2> points{};
    if (geom.volume_fractions(index) == 1.0) {
      const Coordinates<2> xC = GetAbsoluteVolumeCentroid(geom, index, dx);
      points.push_back(ShiftFrom(xC, -0.5 * dx[0], -0.5 * dx[1]));
      points.push_back(ShiftFrom(xC, -0.5 * dx[0], +0.5 * dx[1]));
      points.push_back(ShiftFrom(xC, +0.5 * dx[0], -0.5 * dx[1]));
      points.push_back(ShiftFrom(xC, +0.5 * dx[0], +0.5 * dx[1]));
    } else if (geom.volume_fractions(index) > 0.0) {
      const auto ix = 0;
      const auto iy = 1;
      const Index<2> fL = index;
      const Index<2> fR = Shift(index, Direction::X, 1);
      const double betaL = geom.face_fractions[ix](fL);
      const double betaR = geom.face_fractions[ix](fR);
      const double betaLy = geom.face_fractions[iy](fL);
      const double betaRy =
          geom.face_fractions[iy](Shift(index, Direction::Y, 1));
      const Coordinates<2> xC{(double(index[0]) + 0.5) * dx[0] , (double(index[1]) + 0.5) * dx[1]};
      const Coordinates<2> xB = GetAbsoluteBoundaryCentroid(geom, index, dx);
      const Coordinates<2> xN = GetBoundaryNormal(geom, index);

      const double dBeta_x = betaR - betaL;
      const double dBeta_y = betaRy - betaLy;

      Eigen::Vector2d xLL = xB;
      Eigen::Vector2d xLR = xB;

      xLL[iy] += 0.5 * dBeta_x * dx[iy];
      xLL[ix] -= 0.5 * dBeta_y * dx[ix];

      xLR[iy] -= 0.5 * dBeta_x * dx[iy];
      xLR[ix] += 0.5 * dBeta_y * dx[ix];

      points.push_back(Point_2{xLL[0], xLL[1]});
      points.push_back(Point_2{xLR[0], xLR[1]});

      auto e_x = Coordinates<2>::Unit(0);
      auto e_y = Coordinates<2>::Unit(1);
      std::vector<Coordinates<2>> xs{};
      xs.push_back(xC - 0.5 * dx);
      xs.push_back(xC + 0.5 * dx);
      xs.push_back(xC + 0.5 * dx[0] * e_x - 0.5 * dx[1] * e_y);
      xs.push_back(xC - 0.5 * dx[0] * e_x + 0.5 * dx[1] * e_y);
      for (Coordinates<2> x : xs) {
        if (xN.dot(xB - x) < 0) {
          points.emplace_back(x[0], x[1]);
        }
      }
    }
    Polygon_2 polygon{};
    CGAL::convex_hull_2(points.begin(), points.end(),
                        std::back_inserter(polygon));
    return polygon;
  }

  Polygon_2 GetHGridPolygonAtBoundary(const CutCellData<2>& geom,
                                      const Index<2>& index,
                                      const Eigen::Vector2d& dx,
                                      Direction dir) {
    FUB_ASSERT(dir == Direction::X || dir == Direction::Y);
    const auto ix = static_cast<std::size_t>(dir);
    const auto iy = 1 - ix;
    const Index<Rank> fL = index;
    const Index<Rank> fR = Shift(index, dir, 1);
    const double betaL = geom.face_fractions[ix](fL);
    const double betaR = geom.face_fractions[ix](fR);
    const double betaLy = geom.face_fractions[iy](fL);
    const double betaRy =
        geom.face_fractions[iy](Shift(index, Direction(iy), 1));
    const Coordinates<Rank> xB = GetAbsoluteBoundaryCentroid(geom, index, dx);
    const Coordinates<Rank> xN = GetBoundaryNormal(geom, index);

    const double dBeta_x = betaR - betaL;
    const double dBeta_y = betaRy - betaLy;
    const double alpha = betaL > betaR
                             ? 0.5 + geom.boundary_centeroids(index, ix)
                             : 0.5 - geom.boundary_centeroids(index, ix);

    const int sign = (dBeta_x > 0) - (dBeta_x < 0);
    Eigen::Vector2d xR = xB;
    xR[ix] -= sign * (1.0 - alpha) * dx[ix];

    Eigen::Vector2d xRL = xR;
    Eigen::Vector2d xRR = xR;
    Eigen::Vector2d xLL = xB;
    Eigen::Vector2d xLR = xB;

    xRR[iy] -= 0.5 * dBeta_x * dx[iy];
    xLR[iy] -= 0.5 * dBeta_x * dx[iy];

    xLL[iy] += 0.5 * dBeta_x * dx[iy];
    xRL[iy] += 0.5 * dBeta_x * dx[iy];

    xLL[ix] -= 0.5 * dBeta_y * dx[ix];
    xLR[ix] += 0.5 * dBeta_y * dx[ix];

    Eigen::Vector2d mRL = xRL - 2.0 * xN.dot(xRL - xLL) * xN;
    Eigen::Vector2d mRR = xRR - 2.0 * xN.dot(xRR - xLR) * xN;

    std::vector<Point_2> points{};
    points.push_back(Point_2{xLL[0], xLL[1]});
    points.push_back(Point_2{xLR[0], xLR[1]});
    points.push_back(Point_2{mRL[0], mRL[1]});
    points.push_back(Point_2{mRR[0], mRR[1]});

    Polygon_2 polygon{};
    CGAL::convex_hull_2(points.begin(), points.end(),
                        std::back_inserter(polygon));

    return polygon;
  }

  IntegrationPoints GetIntegrationPoints(const Index<Rank>& index,
                                         const CutCellData<Rank>& geom,
                                         const Coordinates<Rank>& dx,
                                         Direction dir) {
    // Fill aux data for the cut-cell itself
    IntegrationPoints integration{};
    integration.xB = GetAbsoluteBoundaryCentroid(geom, index, dx);
    integration.iB = index;
    const Polygon_2 hgrid_polygon =
        GetHGridPolygonAtBoundary(geom, index, dx, dir);
    int count = 0;
    ForEachIndex(Neighborhood<Rank>(index, 1), [&](auto... is) {
      Index<Rank> i{is...};
      const Polygon_2 neighbor_polygon = GetPolygon(geom, i, dx);
      if (neighbor_polygon.is_empty()) {
        return;
      }
      if (CGAL::do_intersect(hgrid_polygon, neighbor_polygon)) {
        FUB_ASSERT(count < IntegrationPoints::kMaxSources);
        std::vector<Polygon_with_holes_2> intersections{};
        CGAL::intersection(hgrid_polygon, neighbor_polygon,
                           std::back_inserter(intersections));
        FUB_ASSERT(intersections.size() == 1);
        integration.index[count] = i;
        Eigen::Vector2d xM = Eigen::Vector2d::Zero();
        for (const Point_2& p : intersections[0].outer_boundary()) {
          const double px = p[0].exact().convert_to<double>();
          const double py = p[1].exact().convert_to<double>();
          xM += Eigen::Vector2d{px, py};
        }
        xM /= intersections[0].outer_boundary().size();
        integration.xM[count] = xM;
        integration.volume[count] = intersections[0]
                                        .outer_boundary()
                                        .area()
                                        .exact()
                                        .convert_to<double>();
        count += 1;
      }
    });
    return integration;
  }

  bool IntegrateInteriorCellState(Conservative& integral,
                                  Conservative& integral_gradient,
                                  const View<const Conservative>& states,
                                  const View<const Conservative>& gradient_x,
                                  const View<const Conservative>& gradient_y,
                                  const View<const Conservative>& gradient_z,
                                  const CutCellData<Rank>& geom,
                                  const Index<Rank>& index,
                                  const Coordinates<Rank>& dx, Direction dir) {
    const auto d = static_cast<std::size_t>(dir);
    const Index<Rank> fL = index;
    const Index<Rank> fR = Shift(fL, dir, 1);
    const double betaL = geom.face_fractions[d](fL);
    const double betaLds = geom.doubly_shielded_fractions[d](fL);
    const double betaLsr = geom.shielded_right_fractions[d](fL);
    const double betaR = geom.face_fractions[d](fR);
    const double betaRds = geom.doubly_shielded_fractions[d](fR);
    const double betaRsl = geom.shielded_left_fractions[d](fR);
    if (betaL > betaR && betaLds < betaLsr) {
      //           Here we have the boundary from right
      //
      //        ------------- -------------
      //       |             |             |
      //       |             |             |
      //       |     iL      |      i      |
      //       |             |            /|
      //       |             |          /  |
      //       |             |        /    |
      //        ------------- -------------
      const double alpha = 0.5 + geom.boundary_centeroids(index, d);
      const Index<Rank> iL = Shift(index, dir, -1);
      Load(gradient_[0], gradient_x, iL);
      Load(gradient_[1], gradient_y, iL);
      Load(gradient_[2], gradient_z, iL);
      const Coordinates<Rank> xL = GetAbsoluteVolumeCentroid(geom, iL, dx);
      const Coordinates<Rank> xLssR =
          GetAbsoluteSinglyShieldedFromRightVolumeCentroid(
              geom, fL, Side::Lower, dir, dx);
      span<const Conservative, Rank> grads{gradient_.data(), Rank};
      ApplyGradient(gradient_dir_, grads, xLssR - xL);
      Load(state_, states, iL);
      state_ += gradient_dir_;

      ForEachVariable(
          [&](auto&& y, auto&& dy, const auto& xL, const auto& dxL) {
            y = (1.0 - alpha) * (xL + alpha * dx[d] * dxL);
            dy = (1.0 - alpha) * dxL;
          },
          integral, integral_gradient, state_, grads[d]);

      const Coordinates<Rank> xC = GetAbsoluteVolumeCentroid(geom, index, dx);
      const Coordinates<Rank> xCssR =
          xLssR +
          0.5 * (1.0 - alpha) * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
      Load(gradient_[0], gradient_x, index);
      Load(gradient_[1], gradient_y, index);
      Load(gradient_[2], gradient_z, index);
      ApplyGradient(gradient_dir_, grads, xCssR - xC);
      Load(state_, states, index);
      state_ += gradient_dir_;

      ForEachVariable(
          [&](auto&& y, auto&& dy, const auto& xR, const auto& dxR) {
            y += alpha * xR;
            dy += alpha * dxR;
          },
          integral, integral_gradient, state_, grads[d]);

      return true;

    } else if (betaL < betaR && betaRds < betaRsl) {

      const double alpha = 0.5 - geom.boundary_centeroids(index, d);
      const Index<Rank> iR = Shift(index, dir, +1);
      Load(gradient_[0], gradient_x, iR);
      Load(gradient_[1], gradient_y, iR);
      Load(gradient_[2], gradient_z, iR);
      const Coordinates<Rank> xR = GetAbsoluteVolumeCentroid(geom, iR, dx);
      const Coordinates<Rank> xRssR =
          GetAbsoluteSinglyShieldedFromLeftVolumeCentroid(geom, fR, Side::Upper,
                                                          dir, dx);
      span<const Conservative, Rank> grads{gradient_.data(), Rank};
      ApplyGradient(gradient_dir_, grads, xRssR - xR);
      Load(state_, states, iR);
      state_ += gradient_dir_;

      ForEachVariable(
          [&](auto&& y, auto&& dy, const auto& xR, const auto& dxR) {
            y = (1.0 - alpha) * (xR + alpha * dx[d] * dxR);
            dy = (1.0 - alpha) * dxR;
          },
          integral, integral_gradient, state_, grads[d]);

      const Coordinates<Rank> xC = GetAbsoluteVolumeCentroid(geom, index, dx);
      const Coordinates<Rank> xCssR =
          xRssR -
          0.5 * (1.0 - alpha) * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
      Load(gradient_[0], gradient_x, index);
      Load(gradient_[1], gradient_y, index);
      Load(gradient_[2], gradient_z, index);
      ApplyGradient(gradient_dir_, grads, xCssR - xC);
      Load(state_, states, index);
      state_ += gradient_dir_;

      ForEachVariable(
          [&](auto&& y, auto&& dy, const auto& xL, const auto& dxL) {
            y += alpha * xL;
            dy += alpha * dxL;
          },
          integral, integral_gradient, state_, grads[d]);

      return true;
    }

    return false;
  }

  bool IntegrateCellState(Conservative& integral,
                          Conservative& integral_gradient,
                          const View<const Conservative>& states,
                          const View<const Conservative>& gradient_x,
                          const View<const Conservative>& gradient_y,
                          const View<const Conservative>& gradient_z,
                          const CutCellData<Rank>& geom,
                          const IntegrationPoints& integration,
                          const Coordinates<Rank>& dx, Direction dir) {
    const double total_volume = std::accumulate(integration.volume.begin(),
                                                integration.volume.end(), 0);
    for (std::size_t i = 0;
         i < IntegrationPoints::kMaxSources && integration.volume[i]; ++i) {
      FUB_ASSERT(integration.volume[i] > 0);
      const Index<Rank> index = integration.index[i];
      if (Contains(Box<0>(states), index) &&
          geom.volume_fractions(index) > 0.0) {
        Load(gradient_[0], gradient_x, index);
        Load(gradient_[1], gradient_y, index);
        Load(gradient_[2], gradient_z, index);
        span<const Conservative, Rank> grads{gradient_.data(), Rank};
        const Coordinates<Rank> xC = GetAbsoluteVolumeCentroid(geom, index, dx);
        const Coordinates<Rank> xM = integration.xM[i];
        const Coordinates<Rank> dx = xM - xC;
        ApplyGradient(gradient_dir_, grads, dx);
        Load(state_, states, index);
        state_ += gradient_dir_;
        const Coordinates<Rank> e_d =
            Eigen::Matrix<double, Rank, 1>::Unit(int(dir));
        const Coordinates<Rank> xN = GetBoundaryNormal(geom, index);
        const Coordinates<Rank> e_r = e_d - 2.0 * xN.dot(e_d) * xN;
        ApplyGradient(gradient_dir_, grads, e_r);
        const double volume = integration.volume[i];
        ForEachVariable(
            [volume, total_volume](auto&& u, auto&& grad_u, auto&& u_0,
                                   auto&& grad_u_0) {
              u += volume / total_volume * u_0;
              grad_u += volume / total_volume * grad_u_0;
            },
            integral, integral_gradient, state_, gradient_dir_);
      } else {
        return false;
      }
    }
    return true;
  }

  ConservativeHGridReconstruction(const Equation& equation)
      : equation_(equation) {}

  void ComputeGradients(span<Conservative, 2> gradient,
                        span<const Conservative, 4> states,
                        span<const Coordinates<Rank>, 4> x) {
    ForEachComponent(
        [&](double& grad_x, double& grad_y, double uM, double u1, double u2,
            double u3) {
          std::array<double, 2> grad{0.0, 0.0};
          const std::array<double, 4> quantities{uM, u1, u2, u3};
          Base::ComputeGradients(grad, quantities, x);
          grad_x = grad[0];
          grad_y = grad[1];
        },
        gradient[0], gradient[1], states[0], states[1], states[2], states[3]);
  }

  void ComputeGradients(span<Conservative, 2> gradient,
                        span<const Conservative, 5> states,
                        span<const Coordinates<Rank>, 5> x) {
    ForEachComponent(
        [&](double& grad_x, double& grad_y, double uM, double u1, double u2,
            double u3, double u4) {
          std::array<double, 2> grad{0.0, 0.0};
          const std::array<double, 5> quantities{uM, u1, u2, u3, u4};
          Base::ComputeGradients(grad, quantities, x);
          grad_x = grad[0];
          grad_y = grad[1];
        },
        gradient[0], gradient[1], states[0], states[1], states[2], states[3],
        states[4]);
  }

  void ComputeGradients(const View<Conservative>& gradient_x,
                        const View<Conservative>& gradient_y,
                        const View<Conservative>& gradient_z,
                        const View<const Conservative>& states,
                        const StridedDataView<const char, Rank>& /* flags */,
                        const CutCellData<Rank>& geom,
                        const Coordinates<Rank>& dx) {
    Conservative zero{equation_};
    if constexpr (Rank == 2) {
      std::array<Conservative, 2> gradient{equation_};
      std::array<Conservative, 5> u;
      u.fill(Conservative{equation_});
      const IndexBox<Rank> box =
          Shrink(Shrink(Box<0>(gradient_x), Direction::X, {1, 1}), Direction::Y,
                 {1, 1});
      ForEachIndex(box, [&](int i, int j) {
        ////////////////////////////////////////////////
        // All regular case
        if (geom.volume_fractions(i, j) == 1.0 &&
            geom.volume_fractions(i + 1, j) == 1.0 &&
            geom.volume_fractions(i - 1, j) == 1.0 &&
            geom.volume_fractions(i, j + 1) == 1.0 &&
            geom.volume_fractions(i, j - 1) == 1.0) {
          Load(u[0], AsCons(states), {i - 1, j});
          Load(u[1], AsCons(states), {i + 1, j});
          ForEachComponent(
              [&](double& gradient, double qL, double qR) {
                return gradient = 0.5 * (qR - qL) / dx[0];
              },
              gradient[0], u[0], u[1]);
          Load(u[2], AsCons(states), {i, j - 1});
          Load(u[3], AsCons(states), {i, j + 1});
          ForEachComponent(
              [&](double& gradient, double qL, double qR) {
                return gradient = 0.5 * (qR - qL) / dx[1];
              },
              gradient[1], u[2], u[3]);
          //////////////////////////////////////////////////
          // Cut-Cell case
        } else if (geom.volume_fractions(i, j) == 0.0) {
          SetZero(gradient[0]);
          SetZero(gradient[1]);
        } else {
          std::array<Index<Rank + 1>, 4> edges{
              Index<Rank + 1>{i, j, 0}, Index<Rank + 1>{i + 1, j, 0},
              Index<Rank + 1>{i, j, 1}, Index<Rank + 1>{i, j + 1, 1}};
          std::array<double, 4> betas{};
          std::transform(edges.begin(), edges.end(), betas.begin(),
                         [&geom](const Index<Rank + 1>& ijd) {
                           return geom.face_fractions[ijd[2]](ijd[0], ijd[1]);
                         });
          std::array<Index<Rank>, 4> neighbors{
              Index<Rank>{i - 1, j}, Index<Rank>{i + 1, j},
              Index<Rank>{i, j - 1}, Index<Rank>{i, j + 1}};
          std::array<double, 4> alphas{};
          std::transform(neighbors.begin(), neighbors.end(), alphas.begin(),
                         [&geom](const Index<Rank>& ij) {
                           return geom.volume_fractions(ij);
                         });
          std::array<int, 4> is{0, 1, 2, 3};
          std::sort(is.begin(), is.end(), [&](int i, int j) {
            return betas[i] >= betas[j] && alphas[i] >= alphas[j];
          });
          if (betas[is[2]] <= 0.0) {
            // Set neighbors[is[2]] to corner between neighbors[is[0]] and
            // neighbors[is[1]]
            FUB_ASSERT(betas[is[0]] > 0.0 && betas[is[1]] > 0.0);
            FUB_ASSERT(edges[is[0]][2] != edges[is[1]][2]);
            Index<Rank> corner{};
            corner[edges[is[0]][2]] = neighbors[is[0]][edges[is[0]][2]];
            corner[edges[is[1]][2]] = neighbors[is[1]][edges[is[1]][2]];
            neighbors[is[2]] = corner;
            Load(u[0], AsCons(states), {i, j});
            Load(u[1], AsCons(states), neighbors[is[0]]);
            Load(u[2], AsCons(states), neighbors[is[1]]);
            Load(u[3], AsCons(states), neighbors[is[2]]);
            std::array<Coordinates<Rank>, 4> xM;
            xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
            xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
            xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
            xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
            ComputeGradients(
                gradient,
                fub::span<const Conservative>(u).template subspan<0, 4>(), xM);
            // Load(u[0], AsCons(states), {i, j});
            // Load(u[1], AsCons(states), neighbors[is[0]]);
            // Load(u[2], AsCons(states), neighbors[is[1]]);
            // u[3] = u[1];
            // u[4] = u[2];
            // const Coordinates<Rank> xB =
            //     GetAbsoluteBoundaryCentroid(geom, {i, j}, dx);
            // const Coordinates<Rank> n = GetBoundaryNormal(geom, {i, j});
            // std::array<Coordinates<Rank>, 5> xM;
            // xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
            // xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
            // xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
            // xM[3] = xB + ComputeReflectedCoordinates(xM[1] - xB, n);
            // xM[4] = xB + ComputeReflectedCoordinates(xM[2] - xB, n);
            // ComputeGradients(gradient, span(u), xM);
          } else if (betas[is[3]] > 0.0) {
            FUB_ASSERT(betas[is[0]] > 0.0 && betas[is[2]] > 0.0 &&
                       betas[is[1]] > 0.0);
            Load(u[0], AsCons(states), {i, j});
            Load(u[1], AsCons(states), neighbors[is[0]]);
            Load(u[2], AsCons(states), neighbors[is[1]]);
            Load(u[3], AsCons(states), neighbors[is[2]]);
            Load(u[4], AsCons(states), neighbors[is[3]]);
            std::array<Coordinates<Rank>, 5> xM;
            xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
            xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
            xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
            xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
            xM[4] = GetAbsoluteVolumeCentroid(geom, neighbors[is[3]], dx);
            ComputeGradients(gradient, u, xM);
          } else if (betas[is[2]] > 0.0) {
            FUB_ASSERT(betas[is[0]] > 0.0 && betas[is[1]] > 0.0);
            Load(u[0], AsCons(states), {i, j});
            Load(u[1], AsCons(states), neighbors[is[0]]);
            Load(u[2], AsCons(states), neighbors[is[1]]);
            Load(u[3], AsCons(states), neighbors[is[2]]);
            std::array<Coordinates<Rank>, 4> xM;
            xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
            xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
            xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
            xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
            ComputeGradients(
                gradient,
                fub::span<const Conservative>(u).template subspan<0, 4>(), xM);
          }
        }
        Store(gradient_x, gradient[0], {i, j});
        Store(gradient_y, gradient[1], {i, j});
        Store(gradient_z, zero, {i, j});
      });
      // ForEachComponent(
      //     [&](const StridedDataView<double, 2>& u_x,
      //         const StridedDataView<double, 2>& u_y,
      //         const StridedDataView<const double, 2>& u) {
      //       std::array<StridedDataView<double, 2>, 2>
      //       grad_u{u_x.Subview(box),
      //                                                        u_y.Subview(box)};
      //       Base::LimitGradients(grad_u, u.Subview(box), flags, geom, dx);
      //     },
      //     gradient_x, gradient_y, states);
    }
  }

  void SetZero(Conservative& cons) {
    ForEachComponent([](auto&& u) { u = 0.0; }, cons);
  }

  void ReconstructSinglyShieldedStencil(
      span<Complete, 2> h_grid_singly_shielded,
      span<Conservative, 2> h_grid_singly_shielded_gradients,
      span<const Complete, 2> h_grid_eb,
      span<const Conservative, 2> h_grid_eb_gradients,
      const View<const Complete>& states,
      const View<const Conservative>& gradient_x,
      const View<const Conservative>& gradient_y,
      const View<const Conservative>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, const Coordinates<Rank>& dx, Direction dir) {
    const int d = static_cast<std::size_t>(dir);
    const Index<Rank> face_L = cell;
    const Index<Rank> face_R = Shift(face_L, dir, 1);
    const double betaL = cutcell_data.face_fractions[d](face_L);
    const double betaR = cutcell_data.face_fractions[d](face_R);

    const Index<Rank + 1> cell_d = EmbedIndex<Rank>(cell, dir);

    std::array<View<const Conservative>, 3> gradients{gradient_x, gradient_y,
                                                      gradient_z};
    if (betaL < betaR) {
      const double alpha = 0.5 - cutcell_data.boundary_centeroids(cell_d);

      const Index<Rank> iR = Shift(cell, dir, +1);
      Load(gradient_[0], gradient_x, iR);
      Load(gradient_[1], gradient_y, iR);
      Load(gradient_[2], gradient_z, iR);
      const Coordinates<Rank> xR =
          GetAbsoluteVolumeCentroid(cutcell_data, iR, dx);
      const Coordinates<Rank> xB =
          GetAbsoluteBoundaryCentroid(cutcell_data, cell, dx);
      const Coordinates<Rank> xRssR =
          xB + (alpha + 0.5) * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
      span<const Conservative, Rank> grads{gradient_.data(), Rank};
      ApplyGradient(gradient_dir_, grads, xRssR - xR);
      Load(h_grid_singly_shielded[1], states, iR);
      ForEachComponent([](double& x, double dx) { x += dx; },
                       AsCons(h_grid_singly_shielded[1]), gradient_dir_);
      CompleteFromCons(equation_, h_grid_singly_shielded[1],
                       h_grid_singly_shielded[1]);
      Load(h_grid_singly_shielded_gradients[1], gradients[d], iR);

      Load(state_, AsCons(states), cell);
      Load(gradient_[0], gradient_x, cell);
      Load(gradient_[1], gradient_y, cell);
      Load(gradient_[2], gradient_z, cell);
      const Coordinates<Rank> xC =
          GetAbsoluteVolumeCentroid(cutcell_data, cell, dx);
      const Coordinates<Rank> xCssR =
          xB + 0.5 * alpha * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
      ApplyGradient(gradient_dir_, grads, xCssR - xC);
      ForEachComponent([](double& x, double dx) { x += dx; }, AsCons(state_),
                       gradient_dir_);
      ForEachComponent(
          [&](double& Q_l, double& dQ_l, double Q_b, double dQ_b, double Q_i,
              double dQ_i) {
            Q_l = (1 - alpha) * (Q_b + alpha * dx[d] * dQ_b) + alpha * Q_i;
            dQ_l = (1 - alpha) * dQ_b + alpha * dQ_i;
          },
          AsCons(h_grid_singly_shielded[0]),
          h_grid_singly_shielded_gradients[0], AsCons(h_grid_eb[0]),
          h_grid_eb_gradients[0], state_, gradient_[d]);
      CompleteFromCons(equation_, h_grid_singly_shielded[0],
                       h_grid_singly_shielded[0]);
    } else if (betaR < betaL) {
      const double alpha = 0.5 + cutcell_data.boundary_centeroids(cell_d);

      const Index<Rank> iL = Shift(cell, dir, -1);
      Load(gradient_[0], gradient_x, iL);
      Load(gradient_[1], gradient_y, iL);
      Load(gradient_[2], gradient_z, iL);
      const Coordinates<Rank> xL =
          GetAbsoluteVolumeCentroid(cutcell_data, iL, dx);
      const Coordinates<Rank> xB =
          GetAbsoluteBoundaryCentroid(cutcell_data, cell, dx);
      const Coordinates<Rank> xLssR =
          xB - (alpha + 0.5) * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
      span<const Conservative, Rank> grads{gradient_.data(), Rank};
      ApplyGradient(gradient_dir_, grads, xLssR - xL);
      Load(h_grid_singly_shielded[0], states, iL);
      ForEachComponent([](double& x, double dx) { x += dx; },
                       AsCons(h_grid_singly_shielded[0]), gradient_dir_);
      CompleteFromCons(equation_, h_grid_singly_shielded[0],
                       h_grid_singly_shielded[0]);
      Load(h_grid_singly_shielded_gradients[0], gradients[d], iL);

      Load(state_, AsCons(states), cell);
      Load(gradient_[0], gradient_x, cell);
      Load(gradient_[1], gradient_y, cell);
      Load(gradient_[2], gradient_z, cell);
      const Coordinates<Rank> xC =
          GetAbsoluteVolumeCentroid(cutcell_data, cell, dx);
      const Coordinates<Rank> xCssR =
          xB - 0.5 * alpha * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
      ApplyGradient(gradient_dir_, grads, xCssR - xC);
      ForEachComponent([](double& x, double dx) { x += dx; }, AsCons(state_),
                       gradient_dir_);
      ForEachComponent(
          [&](double& Q_r, double& dQ_r, double Q_b, double dQ_b, double Q_i,
              double dQ_i) {
            Q_r = (1 - alpha) * (Q_b + alpha * dx[d] * dQ_b) + alpha * Q_i;
            dQ_r = (1 - alpha) * dQ_b + alpha * dQ_i;
          },
          AsCons(h_grid_singly_shielded[1]),
          h_grid_singly_shielded_gradients[1], AsCons(h_grid_eb[1]),
          h_grid_eb_gradients[1], state_, gradient_[d]);
      CompleteFromCons(equation_, h_grid_singly_shielded[1],
                       h_grid_singly_shielded[1]);
    }
  }

  // Here, we reconstruct a Riemman Problem, where a specified mass flux is
  // required on the boundary and the interior state is known.
  //
  // Using the exact Riemann Problem, we compute the missing state behind the
  // boundary such that the required mass flux is obtained on the cut-cell
  // boundary.
  template <typename ReconstructionMethod>
  void FindRiemannProblemForRequiredMassFlux(
      ReconstructionMethod& reconstruction,
      span<Complete, 2> h_grid_embedded_boundary,
      span<Conservative, 2> h_grid_embedded_boundary_slopes,
      double required_massflux, const View<const Complete>& states,
      const View<const Conservative>& gradient_x,
      const View<const Conservative>& gradient_y,
      const View<const Conservative>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration dt, Eigen::Matrix<double, Rank, 1> dx, Direction dir) {
    const Coordinates<Rank> normal = GetBoundaryNormal(cutcell_data, cell);
    const int dir_v = static_cast<int>(dir);
    const Coordinates<Rank> unit_vector =
        ((normal[dir_v] >= 0) - (normal[dir_v] < 0)) *
        Eigen::Matrix<double, Rank, 1>::Unit(dir_v);
    SetZero(interior_state_);
    SetZero(interior_gradient_);
    IntegrateInteriorCellState(interior_state_, interior_gradient_, states,
                               gradient_x, gradient_y, gradient_z, cutcell_data,
                               cell, dx, dir);
    RequireMassflow_SolveExactRiemannProblem require_massflow{};

    const Side boundary_side = normal[dir_v] >= 0 ? Side::Lower : Side::Upper;
    reconstruction.Reconstruct(reconstructed_interior_state_, interior_state_,
                               interior_gradient_, dt, dx[dir_v], dir,
                               boundary_side);

    if (boundary_side == Side::Lower) {
      require_massflow(equation_, boundary_state_, interior_state_,
                       required_massflux, 1.0, dir);
      require_massflow(equation_, reconstructed_boundary_state_,
                       reconstructed_interior_state_, required_massflux, 1.0,
                       dir);
    } else {
      Reflect(reflected_interior_state_, interior_state_, unit_vector,
              equation_);
      require_massflow(equation_, boundary_state_, reflected_interior_state_,
                       required_massflux, 1.0, dir);
      Reflect(boundary_state_, boundary_state_, unit_vector, equation_);

      Reflect(reflected_interior_state_, reconstructed_interior_state_,
              unit_vector, equation_);
      require_massflow(equation_, reconstructed_boundary_state_,
                       reflected_interior_state_, required_massflux, 1.0, dir);
      Reflect(reconstructed_boundary_state_, reconstructed_boundary_state_,
              unit_vector, equation_);
    }

    ForEachComponent(
        [&](double& gradient, double x_rec, double x) {
          gradient = (x_rec - x) * 0.5 / dx[dir_v];
        },
        boundary_gradient_, reconstructed_boundary_state_, boundary_state_);

    const std::size_t d = static_cast<std::size_t>(dir);
    const Index<Rank> face_L = cell;
    const Index<Rank> face_R = Shift(face_L, dir, 1);
    const double betaL = cutcell_data.face_fractions[d](face_L);
    const double betaR = cutcell_data.face_fractions[d](face_R);
    if (betaL < betaR) {
      ForEachComponent([](auto&& u) { u *= -1.0; }, boundary_gradient_);

      CompleteFromCons(equation_, h_grid_embedded_boundary[0], boundary_state_);
      h_grid_embedded_boundary_slopes[0] = boundary_gradient_;
      CompleteFromCons(equation_, h_grid_embedded_boundary[1], interior_state_);
      h_grid_embedded_boundary_slopes[1] = interior_gradient_;

    } else if (betaR < betaL) {
      ForEachComponent([](auto&& u) { u *= -1.0; }, interior_gradient_);

      CompleteFromCons(equation_, h_grid_embedded_boundary[1], boundary_state_);
      h_grid_embedded_boundary_slopes[1] = boundary_gradient_;
      CompleteFromCons(equation_, h_grid_embedded_boundary[0], interior_state_);
      h_grid_embedded_boundary_slopes[0] = interior_gradient_;
    }
  }

  void ReconstructEmbeddedBoundaryStencil(
      span<Complete, 2> h_grid_embedded_boundary,
      span<Conservative, 2> h_grid_embedded_boundary_slopes,
      const View<const Complete>& states,
      const View<const Conservative>& gradient_x,
      const View<const Conservative>& gradient_y,
      const View<const Conservative>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, Eigen::Matrix<double, Rank, 1> dx, Direction dir) {
    IntegrationPoints boundary_aux_data =
        GetIntegrationPoints(cell, cutcell_data, dx, dir);
    const Coordinates<Rank> normal = GetBoundaryNormal(cutcell_data, cell);
    SetZero(boundary_state_);
    SetZero(boundary_gradient_);
    if (!IntegrateCellState(boundary_state_, boundary_gradient_, states,
                            gradient_x, gradient_y, gradient_z, cutcell_data,
                            boundary_aux_data, dx, dir)) {
      Load(boundary_state_, AsCons(states), cell);
      SetZero(boundary_gradient_);
    }
    SetZero(interior_state_);
    SetZero(interior_gradient_);
    if (!IntegrateInteriorCellState(interior_state_, interior_gradient_, states,
                                    gradient_x, gradient_y, gradient_z,
                                    cutcell_data, cell, dx, dir)) {
      interior_state_ = boundary_state_;
      interior_gradient_ = boundary_gradient_;
    }
    const std::size_t d = static_cast<std::size_t>(dir);
    const Index<Rank> face_L = cell;
    const Index<Rank> face_R = Shift(face_L, dir, 1);
    const double betaL = cutcell_data.face_fractions[d](face_L);
    const double betaR = cutcell_data.face_fractions[d](face_R);

    Reflect(boundary_state_, boundary_state_, normal, equation_);
    Reflect(boundary_gradient_, boundary_gradient_, normal, equation_);

    if (betaL < betaR) {
      ForEachComponent([](auto&& u) { u *= -1.0; }, boundary_gradient_);

      CompleteFromCons(equation_, h_grid_embedded_boundary[0], boundary_state_);
      h_grid_embedded_boundary_slopes[0] = boundary_gradient_;
      CompleteFromCons(equation_, h_grid_embedded_boundary[1], interior_state_);
      h_grid_embedded_boundary_slopes[1] = interior_gradient_;

    } else if (betaR < betaL) {
      ForEachComponent([](auto&& u) { u *= -1.0; }, interior_gradient_);

      CompleteFromCons(equation_, h_grid_embedded_boundary[1], boundary_state_);
      h_grid_embedded_boundary_slopes[1] = boundary_gradient_;
      CompleteFromCons(equation_, h_grid_embedded_boundary[0], interior_state_);
      h_grid_embedded_boundary_slopes[0] = interior_gradient_;
    }
  }

  void ReconstructRegularStencil(span<Complete, 2> h_grid_regular,
                                 span<Conservative, 2> h_grid_regular_gradients,
                                 const View<const Complete>& states,
                                 const View<const Conservative>& gradient_x,
                                 const View<const Conservative>& gradient_y,
                                 const View<const Conservative>& gradient_z,
                                 const CutCellData<Rank>& geom,
                                 const Index<Rank>& face, Duration /*dt*/,
                                 Eigen::Matrix<double, Rank, 1> dx,
                                 Direction dir) {
    Index<Rank> iL = LeftTo(face, dir);
    Index<Rank> iR = RightTo(face, dir);

    const Coordinates<Rank> face_xM =
        GetAbsoluteUnshieldedCentroid(geom, face, dx, dir);
    const Coordinates<Rank> xL_us = Shift(face_xM, dir, -0.5 * dx[int(dir)]);
    const Coordinates<Rank> xR_us = Shift(face_xM, dir, +0.5 * dx[int(dir)]);
    const Coordinates<Rank> xL = GetAbsoluteVolumeCentroid(geom, iL, dx);
    const Coordinates<Rank> xR = GetAbsoluteVolumeCentroid(geom, iR, dx);
    Load(state_, AsCons(states), iL);
    Load(gradient_[0], gradient_x, iL);
    Load(gradient_[1], gradient_y, iL);
    Load(gradient_[2], gradient_z, iL);

    span<const Conservative, Rank> grads{gradient_.data(), Rank};

    const Coordinates<Rank> delta_xL = xL_us - xL;
    ApplyGradient(gradient_dir_, grads, delta_xL);
    state_ += gradient_dir_;
    CompleteFromCons(equation_, h_grid_regular[0], state_);
    h_grid_regular_gradients[0] = gradient_[int(dir)];

    Load(state_, AsCons(states), iR);
    Load(gradient_[0], gradient_x, iR);
    Load(gradient_[1], gradient_y, iR);
    Load(gradient_[2], gradient_z, iR);

    const Coordinates<Rank> delta_xR = xR_us - xR;
    ApplyGradient(gradient_dir_, grads, delta_xR);
    state_ += gradient_dir_;
    CompleteFromCons(equation_, h_grid_regular[1], state_);
    h_grid_regular_gradients[1] = gradient_[int(dir)];
  }

  Equation equation_;
  std::array<Conservative, 3> gradient_{Conservative(equation_),
                                        Conservative(equation_),
                                        Conservative(equation_)};
  std::array<Conservative, 4> stencil{
      Conservative(equation_), Conservative(equation_), Conservative(equation_),
      Conservative(equation_)};
  Conservative limited_slope_{equation_};
  Conservative state_{equation_};
  Conservative gradient_dir_{equation_};
  Conservative reconstructed_boundary_state_{equation_};
  Conservative boundary_state_{equation_};
  Conservative boundary_gradient_{equation_};
  Conservative reconstructed_interior_state_{equation_};
  Conservative reflected_interior_state_{equation_};
  Conservative interior_state_{equation_};
  Conservative interior_gradient_{equation_};
};

template <typename EquationT, typename FluxMethod,
          typename HGridReconstruction =
              ConservativeHGridReconstruction<EquationT>>
class MyCutCellMethod : public FluxMethod {
public:
  using Equation = EquationT;

  // Typedefs
  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;

  // Static Variables

  static constexpr int Rank = Equation::Rank();
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  static constexpr int StencilWidth = FluxMethod::GetStencilWidth();
  static_assert(StencilWidth > 0);
  static constexpr std::size_t StencilSize =
      static_cast<std::size_t>(2 * StencilWidth);

  template <typename T> using DataView = PatchDataView<T, Rank, layout_stride>;

  // Constructors

  /// Constructs a CutCell method from a given base flux method.
  ///
  /// This constructor uses a default constructed riemann problem solver.
  explicit MyCutCellMethod(const Equation& equation);
  MyCutCellMethod(const Equation& equation, const FluxMethod& flux_method);

  using FluxMethod::ComputeStableDt;

  /// \todo compute stable dt inside of cutcells, i.e. in the reflection with
  /// their boundary state.
  double ComputeStableDt(const View<const Complete>& states,
                         const CutCellData<Rank>& cutcell_data, double dx,
                         Direction dir);

  void ComputeRegularFluxes(const View<Conservative>& regular_fluxes,
                            const View<const Complete>& states,
                            const CutCellData<Rank>& cutcell_data, Duration dt,
                            double dx, Direction dir);

  void ComputeCutCellFluxes(
      const View<Conservative>& stabilised_fluxes,
      const View<Conservative>& shielded_left_fluxes,
      const View<Conservative>& shielded_right_fluxes,
      const View<Conservative>& doubly_shielded_fluxes,
      const View<Conservative>& regular_fluxes,
      const View<Conservative>& boundary_fluxes,
      const PatchDataView<double, Rank + 1>& boundary_massflows,
      const View<const Conservative>& gradient_x,
      const View<const Conservative>& gradient_y,
      const View<const Conservative>& gradient_z,
      const View<const Complete>& states, const CutCellData<Rank>& geom,
      Duration dt, const Eigen::Matrix<double, Rank, 1>& dx, Direction dir);

  void ComputeGradients(const View<Conservative>& gradient_x,
                        const View<Conservative>& gradient_y,
                        const View<Conservative>& gradient_z,
                        const View<const Conservative>& states,
                        const StridedDataView<const char, Rank>& flags,
                        const CutCellData<Rank>& geom,
                        const Coordinates<Rank>& dx) {
    h_grid_reconstruction_.ComputeGradients(gradient_x, gradient_y, gradient_z,
                                            states, flags, geom, dx);
  }

  static constexpr int GetStencilWidth() noexcept {
    return FluxMethod::GetStencilWidth() + 1;
  }

private:
  bool IsBoundaryMassflowRequired(
      const PatchDataView<double, 3>& /* boundary_massflows */,
      const Coordinates<2>& /* normal */, Index<2> /* cell */,
      Direction /* dir */) {
    // const int ix = static_cast<int>(dir);
    // const int iy = 1 - ix;
    // const double m_y = boundary_massflows(cell, iy);
    // if (normal[ix] > 0 && !std::isnan(m_y)) {
    //   boundary_massflows(cell, ix) = -normal[iy] / normal[ix] * m_y;
    //   return true;
    // }
    return false;
  }

  Equation equation_;
  HGridReconstruction h_grid_reconstruction_;

  std::array<Complete, 2> h_grid_eb_{};
  std::array<Conservative, 2> h_grid_eb_gradients_{};
  std::array<Complete, 2> h_grid_singly_shielded_{};
  std::array<Conservative, 2> h_grid_singly_shielded_gradients_{};
  std::array<Complete, 2> h_grid_regular_{};
  std::array<Conservative, 2> h_grid_regular_gradients_{};

  Conservative boundary_flux_{equation_};
  Conservative singly_shielded_flux_{equation_};
  Conservative regular_flux_{equation_};

  std::array<CompleteArray, StencilSize> stencil_array_{};
  ConservativeArray numeric_flux_array_{equation_};
};

template <typename Equation, typename FluxMethod>
MyCutCellMethod(const Equation&, const FluxMethod&)
    -> MyCutCellMethod<Equation, FluxMethod>;

// IMPLEMENTATION

template <typename Equation, typename FluxMethod, typename HGridReconstruction>
MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::MyCutCellMethod(
    const Equation& eq)
    : MyCutCellMethod(eq, FluxMethod(eq)) {}

template <typename Equation, typename FluxMethod, typename HGridReconstruction>
MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::MyCutCellMethod(
    const Equation& eq, const FluxMethod& flux_method)
    : FluxMethod(flux_method), equation_(eq), h_grid_reconstruction_(eq) {
  h_grid_eb_.fill(Complete(equation_));
  h_grid_eb_gradients_.fill(Conservative(equation_));
  h_grid_singly_shielded_.fill(Complete(equation_));
  h_grid_singly_shielded_gradients_.fill(Conservative(equation_));
  h_grid_regular_.fill(Complete(equation_));
  h_grid_regular_gradients_.fill(Conservative(equation_));
  stencil_array_.fill(CompleteArray(equation_));
}

/// \todo compute stable dt inside of cutcells, i.e. in the reflection with
/// their boundary state.
template <typename Equation, typename FluxMethod, typename HGridReconstruction>
double
MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::ComputeStableDt(
    const View<const Complete>& states, const CutCellData<Rank>& geom,
    double dx, Direction dir) {
  double min_dt = std::numeric_limits<double>::infinity();
  static constexpr int kWidth = StencilWidth;
  IndexBox<Rank> cellbox = Box<0>(states);
  IndexBox<Rank> fluxbox = Shrink(cellbox, dir, {kWidth, kWidth - 1});
  View<const Complete> base = Subview(states, cellbox);
  using ArrayView = PatchDataView<const double, Rank, layout_stride>;
  ArrayView volumes = geom.volume_fractions.Subview(cellbox);
  const int d = static_cast<int>(dir);
  ArrayView faces = geom.face_fractions[d].Subview(fluxbox);
  std::array<View<const Complete>, StencilSize> stencil_views;
  std::array<ArrayView, StencilSize> stencil_volumes;
  for (std::size_t i = 0; i < StencilSize; ++i) {
    stencil_views[i] =
        Shrink(base, dir,
               {static_cast<std::ptrdiff_t>(i),
                static_cast<std::ptrdiff_t>(StencilSize - i) - 1});
    stencil_volumes[i] = volumes.Subview(
        Shrink(cellbox, dir,
               {static_cast<std::ptrdiff_t>(i),
                static_cast<std::ptrdiff_t>(StencilSize - i) - 1}));
  }
  std::tuple views = std::tuple_cat(std::tuple(faces), AsTuple(stencil_volumes),
                                    AsTuple(stencil_views));
  ForEachRow(views, [this, &min_dt, dx, dir](span<const double> faces,
                                             auto... rest) {
    std::tuple args{rest...};
    std::array<span<const double>, StencilSize> volumes =
        AsArray(Take<StencilSize>(args));
    std::array states = std::apply(
        [](const auto&... xs)
            -> std::array<ViewPointer<const Complete>, StencilSize> {
          return {Begin(xs)...};
        },
        Drop<StencilSize>(args));
    std::array<Array1d, StencilSize> alphas;
    alphas.fill(Array1d::Zero());
    Array1d betas = Array1d::Zero();
    int n = static_cast<int>(faces.size());
    while (n >= kDefaultChunkSize) {
      betas = Array1d::Map(faces.data());
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Load(stencil_array_[i], states[i]);
        alphas[i] = Array1d::Map(volumes[i].data());
      }
      Array1d dts =
          FluxMethod::ComputeStableDt(stencil_array_, betas, alphas, dx, dir);
      min_dt = std::min(min_dt, dts.minCoeff());
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Advance(states[i], kDefaultChunkSize);
        volumes[i] = volumes[i].subspan(kDefaultChunkSize);
      }
      faces = faces.subspan(kDefaultChunkSize);
      n = static_cast<int>(faces.size());
    }
    std::copy_n(faces.data(), n, betas.data());
    std::fill_n(betas.data() + n, kDefaultChunkSize - n, 0.0);
    for (std::size_t i = 0; i < StencilSize; ++i) {
      LoadN(stencil_array_[i], states[i], n);
      std::copy_n(volumes[i].data(), n, alphas[i].data());
      std::fill_n(alphas[i].data() + n, kDefaultChunkSize - n, 0.0);
    }
    Array1d dts =
        FluxMethod::ComputeStableDt(stencil_array_, betas, alphas, dx, dir);
    min_dt = std::min(min_dt, dts.minCoeff());
  });
  return min_dt;
}

template <typename Equation, typename FluxMethod, typename HGridReconstruction>
void MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::
    ComputeRegularFluxes(const View<Conservative>& fluxes,
                         const View<const Complete>& states,
                         const CutCellData<Rank>& cutcell_data, Duration dt,
                         double dx, Direction dir) {
  IndexBox<Rank> fluxbox = Box<0>(fluxes);
  static constexpr int kWidth = FluxMethod::GetStencilWidth();
  IndexBox<Rank> cellbox = Grow(fluxbox, dir, {kWidth, kWidth - 1});
  View<const Complete> base = Subview(states, cellbox);
  using ArrayView = PatchDataView<const double, Rank, layout_stride>;
  ArrayView volumes = cutcell_data.volume_fractions.Subview(cellbox);
  const int d = static_cast<int>(dir);
  ArrayView faces = cutcell_data.face_fractions[d].Subview(fluxbox);
  std::array<View<const Complete>, StencilSize> stencil_views;
  std::array<ArrayView, StencilSize> stencil_volumes;
  for (std::size_t i = 0; i < StencilSize; ++i) {
    stencil_views[i] =
        Shrink(base, dir,
               {static_cast<std::ptrdiff_t>(i),
                static_cast<std::ptrdiff_t>(StencilSize - i) - 1});
    stencil_volumes[i] = volumes.Subview(
        Shrink(cellbox, dir,
               {static_cast<std::ptrdiff_t>(i),
                static_cast<std::ptrdiff_t>(StencilSize - i) - 1}));
  }
  std::tuple views =
      std::tuple_cat(std::tuple(fluxes, faces), AsTuple(stencil_volumes),
                     AsTuple(stencil_views));
  ForEachRow(views, [this, dt, dx, dir](const Row<Conservative>& fluxes,
                                        span<const double> faces,
                                        auto... rest) {
    ViewPointer fit = Begin(fluxes);
    ViewPointer fend = End(fluxes);
    std::tuple args{rest...};
    std::array<span<const double>, StencilSize> volumes =
        AsArray(Take<StencilSize>(args));
    std::array states = std::apply(
        [](const auto&... xs)
            -> std::array<ViewPointer<const Complete>, StencilSize> {
          return {Begin(xs)...};
        },
        Drop<StencilSize>(args));
    std::array<Array1d, StencilSize> alphas;
    alphas.fill(Array1d::Zero());
    Array1d betas = Array1d::Zero();
    int n = static_cast<int>(get<0>(fend) - get<0>(fit));
    while (n >= kDefaultChunkSize) {
      betas = Array1d::Map(faces.data());
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Load(stencil_array_[i], states[i]);
        alphas[i] = Array1d::Map(volumes[i].data());
      }
      FluxMethod::ComputeNumericFlux(numeric_flux_array_, betas, stencil_array_,
                                     alphas, dt, dx, dir);
      for (int i = 0; i < betas.size(); ++i) {
        ForEachComponent(
            [&](auto&& flux [[maybe_unused]]) {
              FUB_ASSERT(betas[i] == 0.0 ||
                         (betas[i] > 0.0 && !std::isnan(flux[i])));
            },
            numeric_flux_array_);
      }
      Store(fit, numeric_flux_array_);
      Advance(fit, kDefaultChunkSize);
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Advance(states[i], kDefaultChunkSize);
        volumes[i] = volumes[i].subspan(kDefaultChunkSize);
      }
      faces = faces.subspan(kDefaultChunkSize);
      n = static_cast<int>(get<0>(fend) - get<0>(fit));
    }
    std::copy_n(faces.data(), n, betas.data());
    std::fill_n(betas.data() + n, kDefaultChunkSize - n, 0.0);
    for (std::size_t i = 0; i < StencilSize; ++i) {
      LoadN(stencil_array_[i], states[i], n);
      std::copy_n(volumes[i].data(), n, alphas[i].data());
      std::fill_n(alphas[i].data() + n, kDefaultChunkSize - n, 0.0);
    }
    FluxMethod::ComputeNumericFlux(numeric_flux_array_, betas, stencil_array_,
                                   alphas, dt, dx, dir);
    StoreN(fit, numeric_flux_array_, n);
  });
}

template <typename Equation, typename FluxMethod, typename HGridReconstruction>
void MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::
    ComputeCutCellFluxes(
        const View<Conservative>& stabilised_fluxes,
        const View<Conservative>& shielded_left_fluxes,
        const View<Conservative>& shielded_right_fluxes,
        const View<Conservative>& /* doubly_shielded_fluxes */,
        const View<Conservative>& regular_fluxes,
        const View<Conservative>& boundary_fluxes,
        const PatchDataView<double, Rank + 1>& boundary_massflows,
        const View<const Conservative>& gradient_x,
        const View<const Conservative>& gradient_y,
        const View<const Conservative>& gradient_z,
        const View<const Complete>& states, const CutCellData<Rank>& geom,
        Duration dt, const Eigen::Matrix<double, Rank, 1>& dx, Direction dir) {

  // Iterate through each EB cell and then
  //
  //   (i) Compute h grid centered at EB and keep the slopes for the reflected
  //   states
  //  (ii) Compute boundary flux using the underlying flux method and its h grid
  // (iii) Compute h grid centered at singly shielded fluxes where boundary
  // states and slopes are interpolated from step (i)
  //  (iv) Compute singly shielded fluxes using the underlying flux method and
  //  its h grid
  const auto d = static_cast<std::size_t>(dir);
  const PatchDataView<const double, Rank>& betas = geom.face_fractions[d];
  const PatchDataView<const double, Rank>& betaUs =
      geom.unshielded_fractions[d];
  ForEachIndex(Shrink(Box<0>(states), dir, {1, 1}), [&](auto... is) {
    Index<Rank> cell{is...};
    Index<Rank> faceL = cell;
    Index<Rank> faceR = Shift(faceL, dir, 1);
    const double betaL = betas(faceL);
    const double betaR = betas(faceR);
    const double betaUsL = betaUs(faceL);
    const double betaUsR = betaUs(faceR);

    if (betaL == betaR) {
      return;
    }

    if (IsBoundaryMassflowRequired(boundary_massflows,
                                   GetBoundaryNormal(geom, cell), cell, dir)) {
      h_grid_reconstruction_.FindRiemannProblemForRequiredMassFlux(
          FluxMethod::GetReconstruction(), h_grid_eb_, h_grid_eb_gradients_,
          boundary_massflows(cell, d), states, gradient_x, gradient_y,
          gradient_z, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_,
                                     h_grid_eb_gradients_, dt, dx[d], dir);
      Store(boundary_fluxes, boundary_flux_, cell);
      for (int i = 0; i < Rank; ++i) {
        boundary_massflows(cell, i) =
            std::numeric_limits<double>::signaling_NaN();
      }
    } else {
      h_grid_reconstruction_.ReconstructEmbeddedBoundaryStencil(
          h_grid_eb_, h_grid_eb_gradients_, states, gradient_x, gradient_y,
          gradient_z, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_,
                                     h_grid_eb_gradients_, dt, dx[d], dir);
      Store(boundary_fluxes, boundary_flux_, cell);
      boundary_massflows(cell, d) = euler::Density(equation_, boundary_flux_);
    }

    h_grid_reconstruction_.ReconstructSinglyShieldedStencil(
        h_grid_singly_shielded_, h_grid_singly_shielded_gradients_, h_grid_eb_,
        h_grid_eb_gradients_, states, gradient_x, gradient_y, gradient_z, geom,
        cell, dt, dx, dir);
    FluxMethod::ComputeNumericFlux(
        singly_shielded_flux_, h_grid_singly_shielded_,
        h_grid_singly_shielded_gradients_, dt, dx[d], dir);

    if (betaL < betaR) {
      if (Contains(Box<0>(shielded_left_fluxes), faceR)) {
        Store(shielded_left_fluxes, singly_shielded_flux_, faceR);
        if (betaUsR > 0.0) {
          h_grid_reconstruction_.ReconstructRegularStencil(
              h_grid_regular_, h_grid_regular_gradients_, states, gradient_x,
              gradient_y, gradient_z, geom, faceR, dt, dx, dir);
          FluxMethod::ComputeNumericFlux(regular_flux_, h_grid_regular_,
                                         h_grid_regular_gradients_, dt, dx[d],
                                         dir);
          Store(regular_fluxes, regular_flux_, faceR);
        }
      }
    } else if (betaR < betaL) {
      if (Contains(Box<0>(shielded_right_fluxes), faceL)) {
        Store(shielded_right_fluxes, singly_shielded_flux_, faceL);
        if (betaUsL > 0.0) {
          h_grid_reconstruction_.ReconstructRegularStencil(
              h_grid_regular_, h_grid_regular_gradients_, states, gradient_x,
              gradient_y, gradient_z, geom, faceL, dt, dx, dir);
          FluxMethod::ComputeNumericFlux(regular_flux_, h_grid_regular_,
                                         h_grid_regular_gradients_, dt, dx[d],
                                         dir);
          Store(regular_fluxes, regular_flux_, faceL);
        }
      }
    }
  });

  // Do the convex combiination of regular and singly shielded fluxes
  ForEachComponent(
      [&](DataView<double> stabilised, DataView<double> shielded_left,
          DataView<double> shielded_right, DataView<const double> regular,
          DataView<const double> boundary) {
        MyStab_ComputeStableFluxComponents(stabilised, shielded_left,
                                           shielded_right, regular, boundary,
                                           geom, dt, dx[d], dir);
      },
      stabilised_fluxes, shielded_left_fluxes, shielded_right_fluxes,
      regular_fluxes, boundary_fluxes);
}

} // namespace fub

#endif
