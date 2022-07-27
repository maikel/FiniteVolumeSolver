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

template <int Rank> struct LinearOptimizationLimiter {
  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const;
};

extern template struct LinearOptimizationLimiter<2>;

template <int Rank> struct NoMdLimiter {
  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const;
};
extern template struct NoMdLimiter<2>;

template <int Rank> struct UpwindMdLimiter {
  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const;
};

extern template struct UpwindMdLimiter<2>;

template <int Rank> struct AnyLimiterBase {
  virtual ~AnyLimiterBase() = default;

  virtual std::unique_ptr<AnyLimiterBase<Rank>> Clone() const = 0;

  virtual void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const = 0;
};

template <int Rank> class AnyLimiter {
public:
  /// @{
  /// \name Constructors

  /// \brief This constructs a method that does nothing on invocation.
  AnyLimiter() = default;

  /// \brief Stores any object which satisfies the tagging method concept.
  template <typename T,
            typename = std::enable_if_t<!decays_to<T, AnyLimiter>()>>
  AnyLimiter(T&& tag);

  /// \brief Copies the implementation.
  AnyLimiter(const AnyLimiter& other);

  /// \brief Copies the implementation.
  AnyLimiter& operator=(const AnyLimiter& other);

  /// \brief Moves the `other` object without allocating and leaves an empty
  /// method.
  AnyLimiter(AnyLimiter&& other) noexcept = default;

  /// \brief Moves the `other` object without allocating and leaves an empty
  /// method.
  AnyLimiter& operator=(AnyLimiter&& other) noexcept = default;
  /// @}

  /// @{
  /// \name Member functions

  /// \brief Mask cells that need further refinement in a regridding procedure.
  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const;
  /// @}

private:
  std::unique_ptr<AnyLimiterBase<Rank>> limiter_;
};

template <int Rank, typename T>
struct LimiterWrapper_ : public AnyLimiterBase<Rank> {

  LimiterWrapper_(const T& limiter) : limiter_{limiter} {}
  LimiterWrapper_(T&& limiter) : limiter_{std::move(limiter)} {}

  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const override {
    limiter_.LimitGradientsAtIndex(grad_u, u, geom, index, dx);
  }

  std::unique_ptr<AnyLimiterBase<Rank>> Clone() const override {
    return std::make_unique<LimiterWrapper_<Rank, T>>(limiter_);
  };

  T limiter_;
};

template <int Rank>
AnyLimiter<Rank>::AnyLimiter(const AnyLimiter& other)
    : limiter_(other.limiter_ ? other.limiter_->Clone() : nullptr) {}

template <int Rank>
AnyLimiter<Rank>& AnyLimiter<Rank>::operator=(const AnyLimiter& other) {
  AnyLimiter tmp(other);
  return *this = std::move(tmp);
}

template <int Rank>
template <typename T, typename>
AnyLimiter<Rank>::AnyLimiter(T&& tag)
    : limiter_{std::make_unique<LimiterWrapper_<Rank, remove_cvref_t<T>>>(
          std::forward<T>(tag))} {}

template <int Rank>
void AnyLimiter<Rank>::LimitGradientsAtIndex(
    const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
    StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
    const Index<Rank>& index, const Coordinates<Rank>& dx) const {
  if (limiter_) {
    return limiter_->LimitGradientsAtIndex(grad_u, u, geom, index, dx);
  }
}

template <int Rank> class BasicHGridReconstruction {
public:
  explicit BasicHGridReconstruction(AnyLimiter<Rank> limiter)
      : limiter_(std::move(limiter)) {}

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
  AnyLimiter<Rank> limiter_;
};

extern template struct BasicHGridReconstruction<2>;

template <typename Equation, typename GradientMethod>
struct HGridReconstruction : BasicHGridReconstruction<Equation::Rank()> {
  using Base = BasicHGridReconstruction<Equation::Rank()>;

  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using Gradient = typename GradientMethod::Gradient;

  static constexpr int Rank = Equation::Rank();

  HGridReconstruction(const Equation& equation, AnyLimiter<Rank> limiter,
                      GradientMethod gradient_method)
      : Base(std::move(limiter)),
        equation_(equation), gradient_method_{std::move(gradient_method)} {}

  void CompleteFromGradient(Complete& dest, const Gradient& source) {
    if constexpr (std::is_same_v<Gradient, Conservative>) {
      CompleteFromCons(equation_, dest, source);
    } else {
      CompleteFromPrim(equation_, dest, source);
    }
  }

  void GradientFromComplete(Gradient& dest, const Complete& src) {
    if constexpr (std::is_same_v<Gradient, Conservative>) {
      dest = AsCons(src);
    } else {
      PrimFromComplete(equation_, dest, src);
    }
  }

  bool IntegrateInteriorCellState(Complete& integral,
                                  Gradient& integral_gradient,
                                  const View<const Complete>& states,
                                  const View<const Gradient>& gradient_x,
                                  const View<const Gradient>& gradient_y,
                                  const View<const Gradient>& gradient_z,
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
      span<const Gradient, Rank> grads{gradient_.data(), Rank};
      ApplyGradient(gradient_dir_, grads, xLssR - xL);
      Load(state_, states, iL);
      GradientFromComplete(scratch_, state_);
      scratch_ += gradient_dir_;

      ForEachVariable(
          [&](auto&& y, auto&& dy, const auto& xL, const auto& dxL) {
            y = (1.0 - alpha) * (xL + 0.5 * alpha * dx[d] * dxL);
            dy = (1.0 - alpha) * dxL;
          },
          integral_scratch_, integral_gradient, scratch_, grads[d]);

      const Coordinates<Rank> xC = GetAbsoluteVolumeCentroid(geom, index, dx);
      const Coordinates<Rank> xCssR =
          xLssR +
          0.5 * (1.0 + alpha) * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
      Load(gradient_[0], gradient_x, index);
      Load(gradient_[1], gradient_y, index);
      Load(gradient_[2], gradient_z, index);
      ApplyGradient(gradient_dir_, grads, xCssR - xC);
      Load(state_, states, index);
      GradientFromComplete(scratch_, state_);
      scratch_ += gradient_dir_;

      ForEachVariable(
          [&](auto&& y, auto&& dy, const auto& xR, const auto& dxR) {
            y += alpha * xR;
            dy += alpha * dxR;
          },
          integral_scratch_, integral_gradient, scratch_, grads[d]);

      CompleteFromGradient(integral, integral_scratch_);

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
      span<const Gradient, Rank> grads{gradient_.data(), Rank};
      ApplyGradient(gradient_dir_, grads, xRssR - xR);
      Load(state_, states, iR);
      GradientFromComplete(scratch_, state_);
      scratch_ += gradient_dir_;

      ForEachVariable(
          [&](auto&& y, auto&& dy, const auto& xR, const auto& dxR) {
            y = (1.0 - alpha) * (xR - 0.5 * alpha * dx[d] * dxR);
            dy = (1.0 - alpha) * dxR;
          },
          integral_scratch_, integral_gradient, scratch_, grads[d]);

      const Coordinates<Rank> xC = GetAbsoluteVolumeCentroid(geom, index, dx);
      const Coordinates<Rank> xCssR =
          xRssR -
          0.5 * (1.0 + alpha) * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
      Load(gradient_[0], gradient_x, index);
      Load(gradient_[1], gradient_y, index);
      Load(gradient_[2], gradient_z, index);
      ApplyGradient(gradient_dir_, grads, xCssR - xC);
      Load(state_, states, index);
      GradientFromComplete(scratch_, state_);
      scratch_ += gradient_dir_;

      ForEachVariable(
          [&](auto&& y, auto&& dy, const auto& xL, const auto& dxL) {
            y += alpha * xL;
            dy += alpha * dxL;
          },
          integral_scratch_, integral_gradient, scratch_, grads[d]);

      CompleteFromGradient(integral, integral_scratch_);
      return true;
    }

    return false;
  }

  bool IntegrateCellState(Complete& integral, Gradient& integral_gradient,
                          const View<const Complete>& states,
                          const View<const Gradient>& gradient_x,
                          const View<const Gradient>& gradient_y,
                          const View<const Gradient>& gradient_z,
                          const CutCellData<Rank>& geom,
                          const HGridIntegrationPoints<Rank>& integration,
                          const Coordinates<Rank>& h, Direction dir) {
    const double total_volume = std::accumulate(integration.volume.begin(),
                                                integration.volume.end(), 0.0);
    if (total_volume == 0.0) {
      return false;
    }
    SetZero(integral_scratch_);
    SetZero(integral_gradient);
    std::array<Gradient, 3> total_grads{
        Gradient(equation_), Gradient(equation_), Gradient(equation_)};
    Coordinates<Rank> total_xM = Coordinates<Rank>::Zero();
    for (std::size_t i = 0;
         i < HGridIntegrationPoints<Rank>::kMaxSources && integration.volume[i];
         ++i) {
      FUB_ASSERT(integration.volume[i] > 0);
      const Index<Rank> index = integration.index[i];
      if (Contains(Box<0>(states), index) &&
          geom.volume_fractions(index) > 0.0) {
        Load(gradient_[0], gradient_x, index);
        Load(gradient_[1], gradient_y, index);
        Load(gradient_[2], gradient_z, index);
        span<const Gradient, Rank> grads{gradient_.data(), Rank};
        const Coordinates<Rank> xC = GetAbsoluteVolumeCentroid(geom, index, h);
        const Coordinates<Rank> xM = integration.xM[i];
        const Coordinates<Rank> dx = xM - xC;
        ApplyGradient(gradient_dir_, grads, dx);
        Load(state_, states, index);
        GradientFromComplete(scratch_, state_);
        scratch_ += gradient_dir_;
        const double volume = integration.volume[i];
        const double lambda = volume / total_volume;
        total_xM += lambda * xM;
        ForEachVariable(
            [lambda](auto&& u, auto&& gradx_u, auto&& grady_u, auto&& gradz_u,
                     auto&& u_0, auto&& gradx_u_0, auto&& grady_u_0,
                     auto&& gradz_u_0) {
              u += lambda * u_0;
              gradx_u += lambda * gradx_u_0;
              grady_u += lambda * grady_u_0;
              gradz_u += lambda * gradz_u_0;
            },
            integral_scratch_, total_grads[0], total_grads[1], total_grads[2],
            scratch_, gradient_[0], gradient_[1], gradient_[2]);
      } else {
        return false;
      }
    }
    const Coordinates<Rank> xN = GetBoundaryNormal(geom, integration.iB);
    const Coordinates<Rank> xB =
        GetAbsoluteBoundaryCentroid(geom, integration.iB, h);
    const Coordinates<Rank> e_d =
        Eigen::Matrix<double, Rank, 1>::Unit(int(dir));
    const Coordinates<Rank> e_r = e_d - 2.0 * xN.dot(e_d) * xN;
    const int sign = (xN[int(dir)] > 0) - (xN[int(dir)] < 0);
    const Coordinates<Rank> x0 = xB - sign * 0.5 * h[int(dir)] * e_r;
    const Coordinates<Rank> dx = x0 - total_xM;
    span<const Gradient, Rank> grads{total_grads.data(), Rank};
    ApplyGradient(gradient_dir_, grads, dx);
    integral_scratch_ += gradient_dir_;
    CompleteFromGradient(integral, integral_scratch_);
    ApplyGradient(integral_gradient, grads, e_r);
    return true;
  }

  void ComputeGradients(span<Gradient, 2> gradient,
                        span<const Complete, 4> states,
                        span<const Coordinates<Rank>, 4> x) {
    std::array<Gradient, 4> states_as_gradient_type;
    states_as_gradient_type.fill(Gradient(equation_));
    for (int i = 0; i < states.size(); ++i) {
      GradientFromComplete(states_as_gradient_type[i], states[size_t(i)]);
    }
    ForEachComponent(
        [&](double& grad_x, double& grad_y, double uM, double u1, double u2,
            double u3) {
          std::array<double, 2> grad{0.0, 0.0};
          const std::array<double, 4> quantities{uM, u1, u2, u3};
          Base::ComputeGradients(grad, quantities, x);
          grad_x = grad[0];
          grad_y = grad[1];
        },
        gradient[0], gradient[1], states_as_gradient_type[0],
        states_as_gradient_type[1], states_as_gradient_type[2],
        states_as_gradient_type[3]);
  }

  void ComputeGradients(span<Gradient, 2> gradient,
                        span<const Complete, 5> states,
                        span<const Coordinates<Rank>, 5> x) {
    std::array<Gradient, 4> states_as_gradient_type;
    states_as_gradient_type.fill(Gradient(equation_));
    for (int i = 0; i < states.size(); ++i) {
      GradientFromComplete(states_as_gradient_type[i], states[size_t(i)]);
    }
    ForEachComponent(
        [&](double& grad_x, double& grad_y, double uM, double u1, double u2,
            double u3, double u4) {
          std::array<double, 2> grad{0.0, 0.0};
          const std::array<double, 5> quantities{uM, u1, u2, u3, u4};
          Base::ComputeGradients(grad, quantities, x);
          grad_x = grad[0];
          grad_y = grad[1];
        },
        gradient[0], gradient[1], states_as_gradient_type[0],
        states_as_gradient_type[1], states_as_gradient_type[2],
        states_as_gradient_type[3], states_as_gradient_type[4]);
  }

  void ComputeGradients(const View<Gradient>& gradient_x,
                        const View<Gradient>& gradient_y,
                        const View<Gradient>& gradient_z,
                        const View<const Complete>& states,
                        const StridedDataView<const char, Rank>& flags,
                        const CutCellData<Rank>& geom,
                        const Coordinates<Rank>& dx) {
    Gradient zero{equation_};
    if constexpr (Rank == 2) {
      std::array<Gradient, 2> gradient{Gradient{equation_},
                                       Gradient{equation_}};
      std::array<Complete, 5> u;
      u.fill(Complete{equation_});
      const IndexBox<Rank> box = Box<0>(gradient_x);
      FUB_ASSERT(box == Box<0>(gradient_y));
      FUB_ASSERT(box == Box<0>(gradient_z));
      ForEachIndex(box, [&](int i, int j) {
        ////////////////////////////////////////////////
        // All regular case
        if (geom.volume_fractions(i, j) == 1.0 &&
            geom.volume_fractions(i + 1, j) == 1.0 &&
            geom.volume_fractions(i - 1, j) == 1.0 &&
            geom.volume_fractions(i, j + 1) == 1.0 &&
            geom.volume_fractions(i, j - 1) == 1.0) {
          Load(u[0], states, {i - 1, j});
          Load(u[1], states, {i, j});
          Load(u[2], states, {i + 1, j});
          gradient_method_.ComputeGradient(
              gradient[0], span<const Complete>{u}.template subspan<0, 3>(),
              dx[0], Direction::X);

          Load(u[0], states, {i, j - 1});
          Load(u[1], states, {i, j});
          Load(u[2], states, {i, j + 1});
          gradient_method_.ComputeGradient(
              gradient[1], span<const Complete>{u}.template subspan<0, 3>(),
              dx[1], Direction::Y);
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
            Load(u[0], states, {i, j});
            Load(u[1], states, neighbors[is[0]]);
            Load(u[2], states, neighbors[is[1]]);
            Load(u[3], states, neighbors[is[2]]);
            std::array<Coordinates<Rank>, 4> xM;
            xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
            xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
            xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
            xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
            ComputeGradients(
                gradient, fub::span<const Complete>(u).template subspan<0, 4>(),
                xM);
          } else if (betas[is[3]] > 0.0) {
            FUB_ASSERT(betas[is[0]] > 0.0 && betas[is[2]] > 0.0 &&
                       betas[is[1]] > 0.0);
            Load(u[0], states, {i, j});
            Load(u[1], states, neighbors[is[0]]);
            Load(u[2], states, neighbors[is[1]]);
            Load(u[3], states, neighbors[is[2]]);
            Load(u[4], states, neighbors[is[3]]);
            std::array<Coordinates<Rank>, 5> xM;
            xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
            xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
            xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
            xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
            xM[4] = GetAbsoluteVolumeCentroid(geom, neighbors[is[3]], dx);
            ComputeGradients(gradient, u, xM);
          } else if (betas[is[2]] > 0.0) {
            FUB_ASSERT(betas[is[0]] > 0.0 && betas[is[1]] > 0.0);
            Load(u[0], states, {i, j});
            Load(u[1], states, neighbors[is[0]]);
            Load(u[2], states, neighbors[is[1]]);
            Load(u[3], states, neighbors[is[2]]);
            std::array<Coordinates<Rank>, 4> xM;
            xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
            xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
            xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
            xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
            ComputeGradients(
                gradient, fub::span<const Complete>(u).template subspan<0, 4>(),
                xM);
          }
        }
        Store(gradient_x, gradient[0], {i, j});
        Store(gradient_y, gradient[1], {i, j});
        Store(gradient_z, zero, {i, j});
      });
      ForEachComponent(
          [&](const StridedDataView<double, 2>& u_x,
              const StridedDataView<double, 2>& u_y,
              const StridedDataView<const double, 2>& u) {
            std::array<StridedDataView<double, 2>, 2> grad_u{u_x.Subview(box),
                                                             u_y.Subview(box)};
            Base::LimitGradients(grad_u, u.Subview(box), flags, geom, dx);
          },
          gradient_x, gradient_y, states);
    }
  }

  template <typename S> void SetZero(S& cons) {
    ForEachComponent([](auto&& u) { u = 0.0; }, cons);
  }

  void ReconstructSinglyShieldedStencil(
      span<Complete, 2> h_grid_singly_shielded,
      span<Gradient, 2> h_grid_singly_shielded_gradients,
      span<const Complete, 2> h_grid_eb,
      span<const Gradient, 2> h_grid_eb_gradients,
      const View<const Complete>& states,
      const View<const Gradient>& gradient_x,
      const View<const Gradient>& gradient_y,
      const View<const Gradient>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, const Coordinates<Rank>& dx, Direction dir) {
    const int d = static_cast<std::size_t>(dir);
    const Index<Rank> face_L = cell;
    const Index<Rank> face_R = Shift(face_L, dir, 1);
    const double betaL = cutcell_data.face_fractions[d](face_L);
    const double betaR = cutcell_data.face_fractions[d](face_R);

    const Index<Rank + 1> cell_d = EmbedIndex<Rank>(cell, dir);

    std::array<View<const Gradient>, 3> gradients{gradient_x, gradient_y,
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
      span<const Gradient, Rank> grads{gradient_.data(), Rank};
      ApplyGradient(gradient_dir_, grads, xRssR - xR);
      Load(h_grid_singly_shielded[1], states, iR);
      GradientFromComplete(scratch_, h_grid_singly_shielded[1]);
      scratch_ += gradient_dir_;
      CompleteFromGradient(h_grid_singly_shielded[1], scratch_);

      Load(h_grid_singly_shielded_gradients[1], gradients[d], iR);

      Load(state_, states, cell);
      Load(gradient_[0], gradient_x, cell);
      Load(gradient_[1], gradient_y, cell);
      Load(gradient_[2], gradient_z, cell);
      const Coordinates<Rank> xC =
          GetAbsoluteVolumeCentroid(cutcell_data, cell, dx);
      const Coordinates<Rank> xCssR =
          xB + 0.5 * alpha * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
      ApplyGradient(gradient_dir_, grads, xCssR - xC);
      GradientFromComplete(scratch_, state_);
      scratch_ += gradient_dir_;

      GradientFromComplete(boundary_gradient_, h_grid_eb[0]);
      ForEachComponent(
          [&](double& Q_l, double& dQ_l, double Q_b, double dQ_b, double Q_i,
              double dQ_i) {
            Q_l =
                (1 - alpha) * (Q_b + 0.5 * alpha * dx[d] * dQ_b) + alpha * Q_i;
            dQ_l = (1 - alpha) * dQ_b + alpha * dQ_i;
          },
          integral_scratch_, h_grid_singly_shielded_gradients[0],
          boundary_gradient_, h_grid_eb_gradients[0], scratch_, gradient_[d]);
      CompleteFromGradient(h_grid_singly_shielded[0], integral_scratch_);
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
      span<const Gradient, Rank> grads{gradient_.data(), Rank};
      ApplyGradient(gradient_dir_, grads, xLssR - xL);
      Load(h_grid_singly_shielded[0], states, iL);
      GradientFromComplete(scratch_, h_grid_singly_shielded[0]);
      scratch_ += gradient_dir_;
      CompleteFromGradient(h_grid_singly_shielded[0], scratch_);

      Load(h_grid_singly_shielded_gradients[0], gradients[d], iL);

      Load(state_, states, cell);
      GradientFromComplete(scratch_, state_);
      Load(gradient_[0], gradient_x, cell);
      Load(gradient_[1], gradient_y, cell);
      Load(gradient_[2], gradient_z, cell);
      const Coordinates<Rank> xC =
          GetAbsoluteVolumeCentroid(cutcell_data, cell, dx);
      const Coordinates<Rank> xCssR =
          xB - 0.5 * alpha * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
      ApplyGradient(gradient_dir_, grads, xCssR - xC);
      scratch_ += gradient_dir_;

      GradientFromComplete(boundary_gradient_, h_grid_eb[1]);
      ForEachComponent(
          [&](double& Q_r, double& dQ_r, double Q_b, double dQ_b, double Q_i,
              double dQ_i) {
            Q_r = (1.0 - alpha) * (Q_b - 0.5 * alpha * dx[d] * dQ_b) +
                  alpha * Q_i;
            dQ_r = (1.0 - alpha) * dQ_b + alpha * dQ_i;
          },
          integral_scratch_, h_grid_singly_shielded_gradients[1],
          boundary_gradient_, h_grid_eb_gradients[1], scratch_, gradient_[d]);
      CompleteFromGradient(h_grid_singly_shielded[1], integral_scratch_);
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
      span<Gradient, 2> h_grid_embedded_boundary_slopes,
      double required_massflux, const View<const Complete>& states,
      const View<const Gradient>& gradient_x,
      const View<const Gradient>& gradient_y,
      const View<const Gradient>& gradient_z,
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

    GradientFromComplete(scratch_, reconstructed_boundary_state_);
    GradientFromComplete(integral_scratch_, boundary_state_);
    ForEachComponent(
        [&](double& gradient, double x_rec, double x) {
          gradient = (x_rec - x) * 0.5 / dx[dir_v];
        },
        boundary_gradient_, scratch_, integral_scratch_);

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
      span<Gradient, 2> h_grid_embedded_boundary_slopes,
      const View<const Complete>& states,
      const View<const Gradient>& gradient_x,
      const View<const Gradient>& gradient_y,
      const View<const Gradient>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, Eigen::Matrix<double, Rank, 1> dx, Direction dir) {
    HGridIntegrationPoints boundary_aux_data =
        GetHGridIntegrationPoints(cutcell_data, cell, dx, dir);
    const Coordinates<Rank> normal = GetBoundaryNormal(cutcell_data, cell);
    SetZero(boundary_state_);
    SetZero(boundary_gradient_);
    if (!IntegrateCellState(boundary_state_, boundary_gradient_, states,
                            gradient_x, gradient_y, gradient_z, cutcell_data,
                            boundary_aux_data, dx, dir)) {
      Load(boundary_state_, states, cell);
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
      CompleteFromCons(equation_, h_grid_embedded_boundary[0], boundary_state_);
      h_grid_embedded_boundary_slopes[0] = boundary_gradient_;
      CompleteFromCons(equation_, h_grid_embedded_boundary[1], interior_state_);
      h_grid_embedded_boundary_slopes[1] = interior_gradient_;

    } else if (betaR < betaL) {
      CompleteFromCons(equation_, h_grid_embedded_boundary[1], boundary_state_);
      h_grid_embedded_boundary_slopes[1] = boundary_gradient_;
      CompleteFromCons(equation_, h_grid_embedded_boundary[0], interior_state_);
      h_grid_embedded_boundary_slopes[0] = interior_gradient_;
    }
  }

  void ReconstructRegularStencil(span<Complete, 2> h_grid_regular,
                                 span<Gradient, 2> h_grid_regular_gradients,
                                 const View<const Complete>& states,
                                 const View<const Gradient>& gradient_x,
                                 const View<const Gradient>& gradient_y,
                                 const View<const Gradient>& gradient_z,
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
    Load(state_, states, iL);
    GradientFromComplete(scratch_, state_);
    Load(gradient_[0], gradient_x, iL);
    Load(gradient_[1], gradient_y, iL);
    Load(gradient_[2], gradient_z, iL);

    span<const Gradient, Rank> grads{gradient_.data(), Rank};

    const Coordinates<Rank> delta_xL = xL_us - xL;
    ApplyGradient(gradient_dir_, grads, delta_xL);
    scratch_ += gradient_dir_;
    CompleteFromGradient(h_grid_regular[0], scratch_);
    h_grid_regular_gradients[0] = gradient_[int(dir)];

    Load(state_, states, iR);
    GradientFromComplete(scratch_, state_);
    Load(gradient_[0], gradient_x, iR);
    Load(gradient_[1], gradient_y, iR);
    Load(gradient_[2], gradient_z, iR);

    const Coordinates<Rank> delta_xR = xR_us - xR;
    ApplyGradient(gradient_dir_, grads, delta_xR);
    scratch_ += gradient_dir_;
    CompleteFromGradient(h_grid_regular[1], scratch_);
    h_grid_regular_gradients[1] = gradient_[int(dir)];
  }

  Equation equation_;
  GradientMethod gradient_method_;
  std::array<Gradient, 3> gradient_{Gradient(equation_), Gradient(equation_),
                                    Gradient(equation_)};
  Complete state_{equation_};
  Gradient scratch_{equation_};
  Gradient integral_scratch_{equation_};
  Gradient gradient_dir_{equation_};
  Complete reconstructed_boundary_state_{equation_};
  Complete boundary_state_{equation_};
  Gradient boundary_gradient_{equation_};
  Complete reconstructed_interior_state_{equation_};
  Complete reflected_interior_state_{equation_};
  Complete interior_state_{equation_};
  Gradient interior_gradient_{equation_};
};

template <typename EquationT, typename FluxMethod>
class MyCutCellMethod : public FluxMethod {
public:
  using Equation = EquationT;

  // Typedefs
  using GradientMethod = typename FluxMethod::GradientMethod;
  using Gradient = typename FluxMethod::Gradient;
  using GradientArray = typename FluxMethod::GradientArray;
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
  MyCutCellMethod(const Equation& equation, AnyLimiter<Rank> limiter);
  MyCutCellMethod(const Equation& equation, const FluxMethod& flux_method,
                  AnyLimiter<Rank> limiter);

  using FluxMethod::ComputeStableDt;

  // Compute
  void PreAdvanceSplitStep(
      const View<Complete>& references, const View<Gradient>& gradient_x,
      const View<Gradient>& gradient_y, const View<Gradient>& gradient_z,
      const View<const Complete>& states, const CutCellData<Rank>& cutcell_data,
      const Coordinates<Rank>& dx, Duration dt, Direction dir, int split_step,
      int total_split_steps);

  /// \todo compute stable dt inside of cutcells, i.e. in the reflection with
  /// their boundary state.
  double ComputeStableDt(const View<const Complete>& states,
                         const CutCellData<Rank>& cutcell_data, double dx,
                         Direction dir);

  void ComputeRegularFluxes(const View<Conservative>& fluxes,
                            const View<const Complete>& states,
                            const View<const Gradient>& gradient_x,
                            const View<const Gradient>& gradient_y,
                            const View<const Gradient>& gradient_z,
                            const CutCellData<Rank>& cutcell_data, Duration dt,
                            double dx, Direction dir);

  void
  ComputeCutCellFluxes(const View<Conservative>& stabilised_fluxes,
                       const View<Conservative>& shielded_left_fluxes,
                       const View<Conservative>& shielded_right_fluxes,
                       const View<Conservative>& doubly_shielded_fluxes,
                       const View<Conservative>& regular_fluxes,
                       const View<Conservative>& boundary_fluxes,
                       const View<const Complete>& boundary_reference_states,
                       const View<const Gradient>& gradient_x,
                       const View<const Gradient>& gradient_y,
                       const View<const Gradient>& gradient_z,
                       const View<const Complete>& states,
                       const CutCellData<Rank>& geom, Duration dt,
                       const Eigen::Matrix<double, Rank, 1>& dx, Direction dir);

  void ComputeGradients(const View<Gradient>& gradient_x,
                        const View<Gradient>& gradient_y,
                        const View<Gradient>& gradient_z,
                        const View<const Complete>& states,
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
  void IntegrateInTime(
    const View<const Complete>& states,
    const Index<Rank>& index, const CutCellData<Rank>& geom,
    const Coordinates<Rank>& dx, Duration dt, Direction dir);

  void ComputeCutCellFlux(const View<const Complete>& states,
                          span<View<const Gradient>, 3> gradients,
                          const Index<Rank>& index,
                          const CutCellData<Rank>& geom,
                          const Coordinates<Rank>& dx, Duration dt,
                          Direction dir);

  void Integrate(Complete& state, const View<const Complete>& states,
               span<View<const Gradient>, 3> gradients,
               const Index<Rank>& index, const CutCellData<Rank>& geom,
               const Coordinates<Rank>& dx, Duration dt, Direction dir);

  void ComputeReferenceState(const View<Complete>& references,
                             const View<const Complete>& states,
                             span<View<const Gradient>, 3> gradients,
                             const Index<Rank>& index,
                             const CutCellData<Rank>& geom,
                             const Coordinates<Rank>& dx, Duration dt);

  void ReconstructOnBoundary(Complete& state,
                             const View<const Complete>& states,
                             span<View<const Gradient>, 3> gradients,
                             const Index<Rank>& index,
                             const CutCellData<Rank>& geom,
                             const Coordinates<Rank>& dx, Duration dt);

  Equation equation_;
  HGridReconstruction<Equation, GradientMethod> h_grid_reconstruction_;

  std::array<Complete, 2> h_grid_eb_{};
  std::array<Gradient, 2> h_grid_eb_gradients_{};
  std::array<Complete, 2> h_grid_singly_shielded_{};
  std::array<Gradient, 2> h_grid_singly_shielded_gradients_{};
  std::array<Complete, 2> h_grid_regular_{};
  std::array<Gradient, 2> h_grid_regular_gradients_{};

  Conservative boundary_flux_{equation_};
  Conservative singly_shielded_flux_{equation_};
  Conservative regular_flux_{equation_};
  Conservative regular_flux_left_{equation_};
  Conservative regular_flux_right_{equation_};
  Conservative stable_flux_{equation_};

  Complete state_{equation_};
  Complete reflected_{equation_};
  Complete solution_{equation_};

  std::array<CompleteArray, 4> stencil_array_{};
  std::array<GradientArray, 2> gradient_array_{};
  ConservativeArray numeric_flux_array_{equation_};
};

// IMPLEMENTATION

template <typename Equation, typename FluxMethod>
MyCutCellMethod<Equation, FluxMethod>::MyCutCellMethod(const Equation& eq,
                                                       AnyLimiter<Rank> limiter)
    : MyCutCellMethod(eq, FluxMethod(eq), std::move(limiter)) {}

template <typename Equation, typename FluxMethod>
MyCutCellMethod<Equation, FluxMethod>::MyCutCellMethod(
    const Equation& eq, const FluxMethod& flux_method, AnyLimiter<Rank> limiter)
    : FluxMethod(flux_method), equation_(eq),
      h_grid_reconstruction_(eq, std::move(limiter),
                             FluxMethod::GetGradientMethod()) {
  h_grid_eb_.fill(Complete(equation_));
  h_grid_eb_gradients_.fill(Gradient(equation_));
  h_grid_singly_shielded_.fill(Complete(equation_));
  h_grid_singly_shielded_gradients_.fill(Gradient(equation_));
  h_grid_regular_.fill(Complete(equation_));
  h_grid_regular_gradients_.fill(Gradient(equation_));
  stencil_array_.fill(CompleteArray(equation_));
}

template <typename Equation, typename FluxMethod>
void MyCutCellMethod<Equation, FluxMethod>::IntegrateInTime(
    [[maybe_unused]] const View<const Complete>& states,
    [[maybe_unused]] const Index<Rank>& index,
    [[maybe_unused]] const CutCellData<Rank>& geom,
    [[maybe_unused]] const Coordinates<Rank>& dx, 
    [[maybe_unused]] Duration dt, 
    [[maybe_unused]] Direction dir)
{
  // Conservative state{equation_};
  // PatchDataView<const  double, Rank> betaUs =
  //     geom.unshielded_fractions_rel[r].Subview(faces);
  // PatchDataView<const double, Rank> betaL =
  //     geom.shielded_left_fractions_rel[r].Subview(faces);
  // PatchDataView<const double, Rank> betaR =
  //     geom.shielded_right_fractions_rel[r].Subview(faces);

  // Index<Rank> fL = LeftTo(index);
  // Index<Rank> fR = RightTo(index);

  // const double alpha = 0;
  // const double dt_over_alpha_h = dt.count() / dx[d] / alpha;
  // const double betaUs = 0;
  // const double betaSS = 0;
  // next = prev + dt_over_h * (betaUs * (fL - fR) + betaSS * (fSS - fB));
  // ForEachVariable([&](auto&& next, auto&& prev, ) {});
}

template <typename Equation, typename FluxMethod>
void MyCutCellMethod<Equation, FluxMethod>::ComputeCutCellFlux(
    const View<const Complete>& states, span<View<const Gradient>, 3> gradients,
    const Index<Rank>& index, const CutCellData<Rank>& geom,
    const Coordinates<Rank>& dx, Duration dt, Direction dir) {

  const int d = static_cast<int>(dir);
  // SetZero(regular_flux_);
  // SetZero(boundary_flux_);
  // SetZero(singly_shielded_flux_);

  const PatchDataView<const double, Rank>& betas = geom.face_fractions[d];
  const PatchDataView<const double, Rank>& betaUs = geom.unshielded_fractions[d];

  Index<Rank> faceL = index;
  Index<Rank> faceR = Shift(faceL, dir, 1);
  const double betaL = betas(faceL);
  const double betaR = betas(faceR);
  const double betaUsL = betaUs(faceL);
  const double betaUsR = betaUs(faceR);

  if (betaUsR > 0.0) {
    h_grid_reconstruction_.ReconstructRegularStencil(
        h_grid_regular_, h_grid_regular_gradients_, states, gradients[0],
        gradients[1], gradients[2], geom, faceR, dt, dx, dir);
    FluxMethod::ComputeNumericFlux(regular_flux_right_, h_grid_regular_,
                                    h_grid_regular_gradients_, dt, dx[d], dir);
  }

  if (betaUsL > 0.0) {
    h_grid_reconstruction_.ReconstructRegularStencil(
        h_grid_regular_, h_grid_regular_gradients_, states, gradients[0],
        gradients[1], gradients[2], geom, faceL, dt, dx, dir);
    FluxMethod::ComputeNumericFlux(regular_flux_left_, h_grid_regular_,
                                    h_grid_regular_gradients_, dt, dx[d], dir);
  }

  if (betaL == betaR) {
    return;
  }

  h_grid_reconstruction_.ReconstructEmbeddedBoundaryStencil(
      h_grid_eb_, h_grid_eb_gradients_, states, gradients[0], gradients[1],
      gradients[2], geom, index, dt, dx, dir);
  FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_,
                                 h_grid_eb_gradients_, dt, dx[d], dir);

  h_grid_reconstruction_.ReconstructSinglyShieldedStencil(
      h_grid_singly_shielded_, h_grid_singly_shielded_gradients_, h_grid_eb_,
      h_grid_eb_gradients_, states, gradients[0], gradients[1], gradients[2],
      geom, index, dt, dx, dir);
  FluxMethod::ComputeNumericFlux(singly_shielded_flux_, h_grid_singly_shielded_,
                                 h_grid_singly_shielded_gradients_, dt, dx[d],
                                 dir);
}

template <typename Equation, typename FluxMethod>
void MyCutCellMethod<Equation, FluxMethod>::Integrate(
    Complete& state, const View<const Complete>& states,
    span<View<const Gradient>, 3> gradients, const Index<Rank>& index,
    const CutCellData<Rank>& geom, const Coordinates<Rank>& dx, Duration dt,
    Direction dir) {
  ComputeCutCellFlux(states, gradients, index, geom, dx, dt, dir);
  // ComputeStableFlux(stable_flux_, regular_flux_, boundary_flux_,
  // singly_shielded_flux_, geom, index, dx, dt, dir);
  IntegrateInTime(state, states, geom, index, dx, dt, dir);
}

template <typename Equation, typename FluxMethod>
void MyCutCellMethod<Equation, FluxMethod>::ReconstructOnBoundary(
    Complete& state, const View<const Complete>& states,
    span<View<const Gradient>, 3> gradients, const Index<Rank>& index,
    const CutCellData<Rank>& geom, const Coordinates<Rank>& dx, Duration dt) {
  Advance(state, states, gradients, index, geom, dx, dt, Direction::X);
  Advance(state, states, gradients, index, geom, dx, dt, Direction::Y);
  // InterpolateToBoundary(state, states, gradients, index, geom, dx);
}

template <typename Equation, typename FluxMethod>
void MyCutCellMethod<Equation, FluxMethod>::ComputeReferenceState(
    const View<Complete>& references, const View<const Complete>& states,
    span<View<const Gradient>, 3> gradients, const Index<Rank>& index,
    const CutCellData<Rank>& geom, const Coordinates<Rank>& dx, Duration dt) {
  FUB_ASSERT(IsCutCell(geom, index));
  ReconstructOnBoundary(state_, states, gradients, geom, dx, 0.5 * dt);
  const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(Direction::X);
  const Eigen::Matrix<double, Rank, 1> normal = GetBoundaryNormal(geom, index);
  Rotate(state_, state_, MakeRotation(normal, unit), equation_);
  Reflect(reflected_, state_, unit, equation_);
  // riemann_solver_.SolveRiemannProblem(solution_, reflected_, state_,
  //                                     Direction::X);
  Rotate(solution_, solution_, MakeRotation(unit, normal), equation_);
  Store(references, solution_, index);
}

template <typename Equation, typename FluxMethod>
void MyCutCellMethod<Equation, FluxMethod>::PreAdvanceSplitStep(
    const View<Complete>& references, const View<Gradient>& gradient_x,
    const View<Gradient>& gradient_y, const View<Gradient>& gradient_z,
    const View<const Complete>& states, const CutCellData<Rank>& geom,
    const Coordinates<Rank>& dx, Duration dt, Direction dir, int split_step,
    int total_split_steps) {
  ComputeGradients(gradient_x, gradient_y, gradient_z, states, geom, dx);
  std::array<View<const Gradient>, 3> gradients{gradient_x, gradient_y,
                                                gradient_z};
  IndexBox<Rank> box =
      Shrink(Shrink(Box<0>(states), Direction::X, 1), Direction::Y, 1);
  ForEachIndex(box, [&](auto... is) {
    Index<Rank> index{is...};
    if (IsCutCell(geom, index)) {
      ComputeReferenceState(references, states, gradients, index, geom,
                            dt);
    }
  });
}

/// \todo compute stable dt inside of cutcells, i.e. in the reflection with
/// their boundary state.
template <typename Equation, typename FluxMethod>
double MyCutCellMethod<Equation, FluxMethod>::ComputeStableDt(
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

template <typename Equation, typename FluxMethod>
void MyCutCellMethod<Equation, FluxMethod>::ComputeRegularFluxes(
    const View<Conservative>& fluxes, const View<const Complete>& states,
    const View<const Gradient>& gradient_x,
    const View<const Gradient>& gradient_y,
    const View<const Gradient>& gradient_z,
    const CutCellData<Rank>& cutcell_data, Duration dt, double dx,
    Direction dir) {
  IndexBox<Rank> fluxbox = Box<0>(fluxes);
  IndexBox<Rank> cellbox = Grow(fluxbox, dir, {1, 0});
  View<const Complete> base = Subview(states, cellbox);
  using ArrayView = PatchDataView<const double, Rank, layout_stride>;
  ArrayView volumes = cutcell_data.volume_fractions.Subview(cellbox);
  const int d = static_cast<int>(dir);
  std::array<View<const Gradient>, 3> grads{gradient_x, gradient_y, gradient_z};
  View<const Gradient> base_grads = Subview(grads[d], cellbox);
  ArrayView faces = cutcell_data.face_fractions[d].Subview(fluxbox);
  std::array<View<const Complete>, 2> stencil_views{};
  std::array<View<const Gradient>, 2> gradient_views{};
  std::array<ArrayView, 2> stencil_volumes{};
  for (std::size_t i = 0; i < 2; ++i) {
    stencil_views[i] = Shrink(base, dir,
                              {static_cast<std::ptrdiff_t>(i),
                               static_cast<std::ptrdiff_t>(2 - i) - 1});
    gradient_views[i] = Shrink(base_grads, dir,
                               {{static_cast<std::ptrdiff_t>(i),
                                 static_cast<std::ptrdiff_t>(2 - i) - 1}});
    stencil_volumes[i] =
        volumes.Subview(Shrink(cellbox, dir,
                               {static_cast<std::ptrdiff_t>(i),
                                static_cast<std::ptrdiff_t>(2 - i) - 1}));
  }
  std::tuple views =
      std::tuple_cat(std::tuple(fluxes, faces), AsTuple(stencil_volumes),
                     AsTuple(stencil_views), AsTuple(gradient_views));
  ForEachRow(views, [this, dt, dx, dir](const Row<Conservative>& fluxes,
                                        span<const double> faces,
                                        span<const double> volumeL,
                                        span<const double> volumeR,
                                        const Row<const Complete>& statesL,
                                        const Row<const Complete>& statesR,
                                        const Row<const Gradient>& gradsL,
                                        const Row<const Gradient>& gradsR) {
    ViewPointer fit = Begin(fluxes);
    ViewPointer fend = End(fluxes);
    std::array<span<const double>, 2> volumes{volumeL, volumeR};
    std::array states{Begin(statesL), Begin(statesR)};
    std::array grads{Begin(gradsL), Begin(gradsR)};
    std::array<Array1d, 2> alphas{Array1d::Zero(), Array1d::Zero()};
    Array1d betas = Array1d::Zero();
    int n = static_cast<int>(get<0>(fend) - get<0>(fit));
    while (n >= kDefaultChunkSize) {
      betas = Array1d::Map(faces.data());
      for (std::size_t i = 0; i < 2; ++i) {
        Load(stencil_array_[i], states[i]);
        Load(gradient_array_[i], grads[i]);
        alphas[i] = Array1d::Map(volumes[i].data());
      }
      FluxMethod::ComputeNumericFlux(
          numeric_flux_array_, betas,
          span{stencil_array_}.template subspan<0, 2>(), gradient_array_,
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
        Advance(grads[i], kDefaultChunkSize);
        volumes[i] = volumes[i].subspan(kDefaultChunkSize);
      }
      faces = faces.subspan(kDefaultChunkSize);
      n = static_cast<int>(get<0>(fend) - get<0>(fit));
    }
    std::copy_n(faces.data(), n, betas.data());
    std::fill_n(betas.data() + n, kDefaultChunkSize - n, 0.0);
    for (std::size_t i = 0; i < 2; ++i) {
      LoadN(stencil_array_[i], states[i], n);
      LoadN(gradient_array_[i], grads[i], n);
      std::copy_n(volumes[i].data(), n, alphas[i].data());
      std::fill_n(alphas[i].data() + n, kDefaultChunkSize - n, 0.0);
    }
    FluxMethod::ComputeNumericFlux(
        numeric_flux_array_, betas,
        span{stencil_array_}.template subspan<0, 2>(), gradient_array_, alphas,
        dt, dx, dir);
    StoreN(fit, numeric_flux_array_, n);
  });
}

// struct HGrids
// {
//   std::array<Complete, 2> States_X;
//   std::array<Complete, 2> States_Y;
//   std::array<Gradient, 2> Gradient_X;
//   std::array<Gradient, 2> Gradient_Y;
// };

// template <typename Equation>
// void AdvectiveFlux(Equation& eq, Conservative<Equation>& flux, const Complete<Equation>& riemannProblemSolution, Direction dir)
// {
//   const auto dir_v = static_cast<int>(dir);
//   eq.Flux(flux, riemannProblemSolution);
//   const double pressure = euler::Pressure(eq, riemannProblemSolution);
//   const double velocity = euler::Velocity(eq, riemannProblemSolution, dir_v);
//   flux.momentum[dir_v] -= pressure;
//   flux.energy -= velocity * pressure;
// }

// template <typename Equation>
// struct HGrids {
//   std::array<Conservative<Equation>, 2> flux;
//   std::array<Complete<Equation>, 2> riemann_problem_solution;
// };

// template <typename Equation, int Rank>
// void FixAdvectiveFluxesForConservation(
//   Equation& eq,
//   HGrids<Equation>& boundary,
//   HGrids<Equation>& singly_shielded,
//   const CutCellData<Rank>& geom,
//   const Index<Rank>& index)
// {
//   const double u = euler::Velocity(eb.riemann_problem_solution[0], 0);
//   const double v = euler::Velocity(eb.riemann_problem_solution[1], 1);
//   FUB_ASSERT(u * v >= 0);
//   Eigen::Vector2d xN = GetBoundaryNormal(geom, index);
//   FUB_ASSERT(xN[0] * xN[1] <= 0);
//   const int in = (u * xN[0] >= 0);
//   const int out = 1 - in;

//   Conservative fadv_in(equation_;);
//   Conservative fadv_out(equation_;);
//   Conservative fadv_out_new(equation_;);

//   AdvectiveFlux(eq, fadv_in, eb.riemann_problem_solution[in], Direction(in));
//   AdvectiveFlux(eq, fadv_out, eb.riemann_problem_solution[out], Direction(out));
//   Reflect(fadv_out_new, fadv_in, xN, eq);
  
//   eb.flux[out] = eb.flux[out] + fadv_out_new - fadv_out;
//   singly_shielded.flux[out] = singly_shielded.flux[out] + (1.0 - d) * (fadv_out_new - fadv_out)
// }

// template <typename Equation, int Rank>
// void ComputeFluxesForConservativeFixup(
//   const std::array<View<Conservative<Equation>>, 2>& stabilised_fluxes,
//   const std::array<View<Conservative<Equation>>, 2>& shielded_left_fluxes,
//   const std::array<View<Conservative<Equation>>, 2>& shielded_right_fluxes,
//   [[maybe_unused]] const std::array<View<Conservative<Equation>>, 2>& doubly_shielded_fluxes,
//   const std::array<View<Conservative<Equation>>, 2>& regular_fluxes,
//   const std::array<View<Conservative<Equation>>, 2>& boundary_fluxes,
//   const std::array<View<const Complete<Equation>>, 2>& riemann_problem_solutions,
//   const CutCellData<Rank>& geom, Duration dt, const Eigen::Matrix<double, Rank, 1>& dx)
// {
//   const IndexBox<Rank> box = Box<0>(boundary_fluxes);
//   HGrids<Equation> eb{};
//   HGrids<Equation> singly_shielded{};
//   ForEachIndex(box, [&](auto... is) {
//     Index<Rank> index{is...};
//     if (IsCutCell(index)) {
//       Load(eb.flux[0], boundary_fluxes[0], index);
//       Load(eb.flux[1], boundary_fluxes[1], index);
//       Load(eb.riemann_problem_solution[0], riemann_problem_solutions[0], index);
//       Load(eb.riemann_problem_solution[1], riemann_problem_solutions[1], index);
//       if (betaLx < betaRx) {
//         Load(singly_shielded.flux[0], shielded_left_fluxes(fRx))
//       } else {
//         Load(singly_shielded.flux[0], shielded_right_fluxes(fLx))
//       }

//       if (betaLy < betaRy) {
//         Load(singly_shielded.flux[1], shielded_left_fluxes(fRy))
//       } else {
//         Load(singly_shielded.flux[1], shielded_right_fluxes(fLy))
//       }
//       FixAdvectiveFluxesForConservation(eq, boundary, singly_shielded, geom, index);

//       Store(regular_fluxes[0], zero, fLx);
//       Store(regular_fluxes[0], zero, fRx);
//       Store(regular_fluxes[1], zero, fLy);
//       Store(regular_fluxes[1], zero, fRy);

//       Store(boundary_fluxes[0], eb.flux[0], index);
//       Store(boundary_fluxes[1], eb.flux[1], index);

//       Store(boundary_fluxes[0], eb.flux[0], index);
//       Store(boundary_fluxes[1], eb.flux[1], index);
//     } else {
//       Store(regular_fluxes[0], zero, fLx);
//       Store(regular_fluxes[0], zero, fRx);
//       Store(regular_fluxes[1], zero, fLy);
//       Store(regular_fluxes[1], zero, fRy);
//     }
//   });
// }

template <typename Equation, typename FluxMethod>
void MyCutCellMethod<Equation, FluxMethod>::ComputeCutCellFluxes(
    const View<Conservative>& stabilised_fluxes,
    const View<Conservative>& shielded_left_fluxes,
    const View<Conservative>& shielded_right_fluxes,
    const View<Conservative>& /* doubly_shielded_fluxes */,
    const View<Conservative>& regular_fluxes,
    const View<Conservative>& boundary_fluxes,
    [[maybe_unused]] const View<const Complete>& boundary_reference_states,
    const View<const Gradient>& gradient_x,
    const View<const Gradient>& gradient_y,
    const View<const Gradient>& gradient_z, const View<const Complete>& states,
    const CutCellData<Rank>& geom, Duration dt,
    const Eigen::Matrix<double, Rank, 1>& dx, Direction dir) {

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

    if (betaUsR > 0.0 && Contains(Box<0>(regular_fluxes), faceR)) {
      h_grid_reconstruction_.ReconstructRegularStencil(
          h_grid_regular_, h_grid_regular_gradients_, states, gradient_x,
          gradient_y, gradient_z, geom, faceR, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(regular_flux_, h_grid_regular_,
                                     h_grid_regular_gradients_, dt, dx[d], dir);
      Store(regular_fluxes, regular_flux_, faceR);
    }

    if (betaUsL > 0.0 && Contains(Box<0>(regular_fluxes), faceL)) {
      h_grid_reconstruction_.ReconstructRegularStencil(
          h_grid_regular_, h_grid_regular_gradients_, states, gradient_x,
          gradient_y, gradient_z, geom, faceL, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(regular_flux_, h_grid_regular_,
                                     h_grid_regular_gradients_, dt, dx[d], dir);
      Store(regular_fluxes, regular_flux_, faceL);
    }

    if (betaL == betaR) {
      return;
    }

    h_grid_reconstruction_.ReconstructEmbeddedBoundaryStencil(
        h_grid_eb_, h_grid_eb_gradients_, states, gradient_x, gradient_y,
        gradient_z, geom, cell, dt, dx, dir);
    FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_,
                                   h_grid_eb_gradients_, dt, dx[d], dir);
    Store(boundary_fluxes, boundary_flux_, cell);

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
      }
    } else if (betaR < betaL) {
      if (Contains(Box<0>(shielded_right_fluxes), faceL)) {
        Store(shielded_right_fluxes, singly_shielded_flux_, faceL);
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
