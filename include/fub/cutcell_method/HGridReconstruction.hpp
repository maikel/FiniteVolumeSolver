#ifndef FUB_CUTCELL_H_GRID_RECONSTRUCTION_HPP
#define FUB_CUTCELL_H_GRID_RECONSTRUCTION_HPP

#include "fub/CutCellData.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"
#include "fub/StateUtil.hpp"
#include "CGAL/Exact_predicates_exact_constructions_kernel.h"

namespace fub {
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

template <typename S> void SetZero(S& state) {
  ForEachComponent([](auto&& u) { u = 0.0; }, state);
}

template <typename Equation, typename GradientT> class HGridIntegration {
public:
  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using Gradient = GradientT;
  using FT = CGAL::Epeck_ft;

  static constexpr int Rank = Equation::Rank();
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  static constexpr int N = HGridIntegrationPoints<Rank>::kMaxSources;
  static constexpr std::size_t sN = static_cast<std::size_t>(N);

  explicit HGridIntegration(const Equation& eq);

  bool IntegrateInteriorCellState(Complete& integral,
                                  Gradient& integral_gradient,
                                  const View<const Complete>& states,
                                  const View<const Gradient>& gradient_x,
                                  const View<const Gradient>& gradient_y,
                                  const View<const Gradient>& gradient_z,
                                  const CutCellData<Rank>& geom,
                                  const Index<Rank>& index,
                                  const Coordinates<Rank>& dx, Direction dir);

  bool IntegrateCellState(
    Complete& integral, Gradient& integral_gradient,
    const std::array<Gradient, sN>& states, const std::array<Gradient, sN>& gradient_x,
    const std::array<Gradient, sN>& gradient_y,
    const std::array<Gradient, sN>& gradient_z, const CutCellData<Rank>& geom,
    const HGridIntegrationPoints<Rank>& integration, const Coordinates<Rank>& h,
    Direction dir);

  bool IntegrateCellState(Complete& integral, Gradient& integral_gradient,
                          const View<const Complete>& states,
                          const View<const Gradient>& gradient_x,
                          const View<const Gradient>& gradient_y,
                          const View<const Gradient>& gradient_z,
                          const CutCellData<Rank>& geom,
                          const HGridIntegrationPoints<Rank>& integration,
                          const Coordinates<Rank>& h, Direction dir);

private:
  Equation equation_;
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

template <typename Equation, typename GradientT> class HGridReconstruction {
public:
  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using Gradient = GradientT;

  static constexpr int Rank = Equation::Rank();

  explicit HGridReconstruction(const Equation& equation);

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
      Duration /*dt*/, const Coordinates<Rank>& dx, Direction dir);

  void ReconstructEmbeddedBoundaryStencil(
      span<Complete, 2> h_grid_embedded_boundary,
      span<Gradient, 2> h_grid_embedded_boundary_slopes,
      const View<const Complete>& states,
      const View<const Gradient>& gradient_x,
      const View<const Gradient>& gradient_y,
      const View<const Gradient>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, Eigen::Matrix<double, Rank, 1> dx, Direction dir);

  void ReconstructRegularStencil(span<Complete, 2> h_grid_regular,
                                 span<Gradient, 2> h_grid_regular_gradients,
                                 const View<const Complete>& states,
                                 const View<const Gradient>& gradient_x,
                                 const View<const Gradient>& gradient_y,
                                 const View<const Gradient>& gradient_z,
                                 const CutCellData<Rank>& geom,
                                 const Index<Rank>& face, Duration /*dt*/,
                                 Eigen::Matrix<double, Rank, 1> dx,
                                 Direction dir);

  HGridIntegration<Equation, GradientT>& GetHGridIntegration() noexcept {
    return hgrid_integration_;
  }

  const HGridIntegration<Equation, GradientT>& GetHGridIntegration() const noexcept {
    return hgrid_integration_;
  }

private:
  Equation equation_;
  HGridIntegration<Equation, GradientT> hgrid_integration_{equation_};

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

///////////////////////////////////////////////////////////////////////////////
//                                             HGridIntegration IMPLEMENTATION

template <typename Equation, typename GradientT>
HGridIntegration<Equation, GradientT>::HGridIntegration(const Equation& equation)
      : equation_(equation) {}

template <typename Equation, typename GradientT>
bool HGridIntegration<Equation, GradientT>::IntegrateInteriorCellState(
    Complete& integral, Gradient& integral_gradient,
    const View<const Complete>& states, const View<const Gradient>& gradient_x,
    const View<const Gradient>& gradient_y,
    const View<const Gradient>& gradient_z, const CutCellData<Rank>& geom,
    const Index<Rank>& index, const Coordinates<Rank>& dx, Direction dir) {
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
        GetAbsoluteSinglyShieldedFromRightVolumeCentroid(geom, fL, Side::Lower,
                                                         dir, dx);
    span<const Gradient, Rank> grads{gradient_.data(), Rank};
    ApplyGradient(gradient_dir_, grads, xLssR - xL);
    Load(state_, states, iL);
    StateFromComplete(equation_, scratch_, state_);
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
    StateFromComplete(equation_, scratch_, state_);
    scratch_ += gradient_dir_;

    ForEachVariable(
        [&](auto&& y, auto&& dy, const auto& xR, const auto& dxR) {
          y += alpha * xR;
          dy += alpha * dxR;
        },
        integral_scratch_, integral_gradient, scratch_, grads[d]);

    CompleteFromState(equation_, integral, integral_scratch_);

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
    StateFromComplete(equation_, scratch_, state_);
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
    StateFromComplete(equation_, scratch_, state_);
    scratch_ += gradient_dir_;

    ForEachVariable(
        [&](auto&& y, auto&& dy, const auto& xL, const auto& dxL) {
          y += alpha * xL;
          dy += alpha * dxL;
        },
        integral_scratch_, integral_gradient, scratch_, grads[d]);

    CompleteFromState(equation_, integral, integral_scratch_);
    return true;
  }

  return false;
}

template <typename Equation, typename GradientT>
bool HGridIntegration<Equation, GradientT>::IntegrateCellState(
    Complete& integral, Gradient& integral_gradient,
    const std::array<Gradient, sN>& states, const std::array<Gradient, sN>& gradient_x,
    const std::array<Gradient, sN>& gradient_y,
    const std::array<Gradient, sN>& gradient_z, const CutCellData<Rank>& geom,
    const HGridIntegrationPoints<Rank>& integration, const Coordinates<Rank>& h,
    Direction dir) {
  const FT total_volume = std::accumulate(integration.volume.begin(),
                                              integration.volume.end(), FT{0.0});
  FT total_x = 0.0;
  FT total_y = 0.0;
  for (std::size_t i = 0; i < sN && integration.volume[i]; ++i) {
      const FT volume = integration.volume[i];
      const FT lambda = volume / total_volume;
      const FT x = integration.xM[i][0];
      const FT y = integration.xM[i][1];
      total_x += lambda * x;
      total_y += lambda * y;
  }
  const int n_components = GetNumberOfComponents(equation_, integral_scratch_);
  for (int ncomp = 0; ncomp < n_components; ++ncomp) {
    FT v = 0.0;
    FT dvdx = 0.0;
    FT dvdy = 0.0;
    FT dvdz = 0.0;
    for (std::size_t i = 0; i < sN && integration.volume[i]; ++i) {
      const FT volume = integration.volume[i];
      const FT lambda = volume / total_volume;
      const FT u0 = GetComponent(equation_, states[i], ncomp);
      const FT du0dx = GetComponent(equation_, gradient_x[i], ncomp);
      const FT du0dy = GetComponent(equation_, gradient_y[i], ncomp);
      const FT du0dz = GetComponent(equation_, gradient_z[i], ncomp);
      v += lambda * u0;
      dvdx += lambda * du0dx;
      dvdy += lambda * du0dx;
      dvdz += lambda * du0dx;
    }
    const double u = v.convert_to<double>();
    const double dudx = dvdx.convert_to<double>();
    const double dudy = dvdy.convert_to<double>();
    const double dudz = dvdz.convert_to<double>();
    GetComponent(equation_, integral_scratch_, ncomp) = u;
    GetComponent(equation_, gradient_[0], ncomp) = dudx;
    GetComponent(equation_, gradient_[1], ncomp) = dudy;
    GetComponent(equation_, gradient_[2], ncomp) = dudz;
  }
  const Coordinates<Rank> total_xM{total_x.convert_to<double>(), total_y.convert_to<double>()};
  const Coordinates<Rank> xN = GetBoundaryNormal(geom, integration.iB);
  const Coordinates<Rank> xB =
      GetAbsoluteBoundaryCentroid(geom, integration.iB, h);
  const Coordinates<Rank> e_d = Eigen::Matrix<double, Rank, 1>::Unit(int(dir));
  const Coordinates<Rank> e_r = e_d - 2.0 * xN.dot(e_d) * xN;
  const int sign = (xN[int(dir)] > 0) - (xN[int(dir)] < 0);
  const Coordinates<Rank> x0 = xB - sign * 0.5 * h[int(dir)] * e_r;
  const Coordinates<Rank> dx = x0 - total_xM;
  span<const Gradient, Rank> grads{gradient_.data(), Rank};
  ApplyGradient(gradient_dir_, grads, dx);
  integral_scratch_ += gradient_dir_;
  CompleteFromState(equation_, integral, integral_scratch_);
  ApplyGradient(integral_gradient, grads, e_r);

  return true;
}

template <typename Equation, typename GradientT>
bool HGridIntegration<Equation, GradientT>::IntegrateCellState(
    Complete& integral, Gradient& integral_gradient,
    const View<const Complete>& states, const View<const Gradient>& gradient_x,
    const View<const Gradient>& gradient_y,
    const View<const Gradient>& gradient_z, const CutCellData<Rank>& geom,
    const HGridIntegrationPoints<Rank>& integration, const Coordinates<Rank>& h,
    Direction dir) {
  const double total_volume = std::accumulate(integration.volume.begin(),
                                              integration.volume.end(), double{0.0});
  if (total_volume == 0.0) {
    return false;
  }
  SetZero(integral_scratch_);
  SetZero(integral_gradient);
  std::array<Gradient, sN> u{};
  std::array<Gradient, sN> dudx{};
  std::array<Gradient, sN> dudy{};
  std::array<Gradient, sN> dudz{};
  u.fill(Gradient{equation_});
  dudx.fill(Gradient{equation_});
  dudy.fill(Gradient{equation_});
  dudz.fill(Gradient{equation_});
  for (std::size_t i = 0; i < sN && integration.volume[i]; ++i) {
    FUB_ASSERT(integration.volume[i] > 0);
    const Index<Rank> index = integration.index[i];
    if (Contains(Box<0>(states), index) && geom.volume_fractions(index) > 0.0) {
      Load(gradient_[0], gradient_x, index);
      Load(gradient_[1], gradient_y, index);
      Load(gradient_[2], gradient_z, index);
      span<const Gradient, Rank> grads{gradient_.data(), Rank};
      const Coordinates<Rank> xC = GetAbsoluteVolumeCentroid(geom, index, h);
      const Coordinates<Rank> xM = integration.xM[i];
      const Coordinates<Rank> dx = xM - xC;
      ApplyGradient(gradient_dir_, grads, dx);
      Load(state_, states, index);
      StateFromComplete(equation_, u[i], state_);
      u[i] += gradient_dir_;
      dudx[i] = gradient_[0];
      dudy[i] = gradient_[1];
      dudz[i] = gradient_[2];
    } else {
      return false;
    }
  }
  return IntegrateCellState(integral, integral_gradient, u, dudx, dudy, dudz, geom, integration, h, dir);
}

///////////////////////////////////////////////////////////////////////////////
//                                          HGridReconstruction IMPLEMENTATION

template <typename Equation, typename GradientT>
HGridReconstruction<Equation, GradientT>::HGridReconstruction(const Equation& equation)
      : equation_(equation) {}

template <typename Equation, typename GradientT>
void HGridReconstruction<Equation, GradientT>::ReconstructSinglyShieldedStencil(
    span<Complete, 2> h_grid_singly_shielded,
    span<Gradient, 2> h_grid_singly_shielded_gradients,
    span<const Complete, 2> h_grid_eb,
    span<const Gradient, 2> h_grid_eb_gradients,
    const View<const Complete>& states, const View<const Gradient>& gradient_x,
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
    StateFromComplete(equation_, scratch_, h_grid_singly_shielded[1]);
    scratch_ += gradient_dir_;
    CompleteFromState(equation_, h_grid_singly_shielded[1], scratch_);

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
    StateFromComplete(equation_, scratch_, state_);
    scratch_ += gradient_dir_;

    StateFromComplete(equation_, boundary_gradient_, h_grid_eb[0]);
    ForEachComponent(
        [&](double& Q_l, double& dQ_l, double Q_b, double dQ_b, double Q_i,
            double dQ_i) {
          Q_l = (1 - alpha) * (Q_b + 0.5 * alpha * dx[d] * dQ_b) + alpha * Q_i;
          dQ_l = (1 - alpha) * dQ_b + alpha * dQ_i;
        },
        integral_scratch_, h_grid_singly_shielded_gradients[0],
        boundary_gradient_, h_grid_eb_gradients[0], scratch_, gradient_[d]);
    CompleteFromState(equation_, h_grid_singly_shielded[0], integral_scratch_);
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
    StateFromComplete(equation_, scratch_, h_grid_singly_shielded[0]);
    scratch_ += gradient_dir_;
    CompleteFromState(equation_, h_grid_singly_shielded[0], scratch_);

    Load(h_grid_singly_shielded_gradients[0], gradients[d], iL);

    Load(state_, states, cell);
    StateFromComplete(equation_, scratch_, state_);
    Load(gradient_[0], gradient_x, cell);
    Load(gradient_[1], gradient_y, cell);
    Load(gradient_[2], gradient_z, cell);
    const Coordinates<Rank> xC =
        GetAbsoluteVolumeCentroid(cutcell_data, cell, dx);
    const Coordinates<Rank> xCssR =
        xB - 0.5 * alpha * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
    ApplyGradient(gradient_dir_, grads, xCssR - xC);
    scratch_ += gradient_dir_;

    StateFromComplete(equation_, boundary_gradient_, h_grid_eb[1]);
    ForEachComponent(
        [&](double& Q_r, double& dQ_r, double Q_b, double dQ_b, double Q_i,
            double dQ_i) {
          Q_r =
              (1.0 - alpha) * (Q_b - 0.5 * alpha * dx[d] * dQ_b) + alpha * Q_i;
          dQ_r = (1.0 - alpha) * dQ_b + alpha * dQ_i;
        },
        integral_scratch_, h_grid_singly_shielded_gradients[1],
        boundary_gradient_, h_grid_eb_gradients[1], scratch_, gradient_[d]);
    CompleteFromState(equation_, h_grid_singly_shielded[1], integral_scratch_);
  }
}

template <typename Equation, typename GradientT>
void HGridReconstruction<Equation, GradientT>::
    ReconstructEmbeddedBoundaryStencil(
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
  if (!hgrid_integration_.IntegrateCellState(boundary_state_, boundary_gradient_, states,
                          gradient_x, gradient_y, gradient_z, cutcell_data,
                          boundary_aux_data, dx, dir)) {
    Load(boundary_state_, states, cell);
    SetZero(boundary_gradient_);
  }
  SetZero(interior_state_);
  SetZero(interior_gradient_);
  if (!hgrid_integration_.IntegrateInteriorCellState(interior_state_, interior_gradient_, states,
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

template <typename Equation, typename GradientT>
void HGridReconstruction<Equation, GradientT>::ReconstructRegularStencil(
    span<Complete, 2> h_grid_regular,
    span<Gradient, 2> h_grid_regular_gradients,
    const View<const Complete>& states, const View<const Gradient>& gradient_x,
    const View<const Gradient>& gradient_y,
    const View<const Gradient>& gradient_z, const CutCellData<Rank>& geom,
    const Index<Rank>& face, Duration /*dt*/, Eigen::Matrix<double, Rank, 1> dx,
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
  StateFromComplete(equation_, scratch_, state_);
  Load(gradient_[0], gradient_x, iL);
  Load(gradient_[1], gradient_y, iL);
  Load(gradient_[2], gradient_z, iL);

  span<const Gradient, Rank> grads{gradient_.data(), Rank};

  const Coordinates<Rank> delta_xL = xL_us - xL;
  ApplyGradient(gradient_dir_, grads, delta_xL);
  scratch_ += gradient_dir_;
  CompleteFromState(equation_, h_grid_regular[0], scratch_);
  h_grid_regular_gradients[0] = gradient_[int(dir)];

  Load(state_, states, iR);
  StateFromComplete(equation_, scratch_, state_);
  Load(gradient_[0], gradient_x, iR);
  Load(gradient_[1], gradient_y, iR);
  Load(gradient_[2], gradient_z, iR);

  const Coordinates<Rank> delta_xR = xR_us - xR;
  ApplyGradient(gradient_dir_, grads, delta_xR);
  scratch_ += gradient_dir_;
  CompleteFromState(equation_, h_grid_regular[1], scratch_);
  h_grid_regular_gradients[1] = gradient_[int(dir)];
}

} // namespace fub

#endif