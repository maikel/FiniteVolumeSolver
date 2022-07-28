#ifndef FUB_CUTCELL_METHOD_H_GRID_REC_2_HPP
#define FUB_CUTCELL_METHOD_H_GRID_REC_2_HPP

#include "fub/cutcell_method/HGridReconstruction.hpp"

namespace fub {

template <typename Equation, typename GradientT>
class HGridReconstruction2 : private HGridReconstruction<Equation, GradientT> {
public:
  using Base = HGridReconstruction<Equation, GradientT>;
  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using Gradient = GradientT;

  static constexpr int Rank = Equation::Rank();

  explicit HGridReconstruction2(const Equation& eq) : Base(eq), equation_{eq} {}

  using Base::ReconstructRegularStencil;
  using Base::ReconstructSinglyShieldedStencil;

  void ReconstructEmbeddedBoundaryStencil(
      span<Complete, 2> h_grid_embedded_boundary,
      span<Gradient, 2> h_grid_embedded_boundary_slopes,
      const View<const Complete>& states,
      const View<const Gradient>& gradient_x,
      const View<const Gradient>& gradient_y,
      const View<const Gradient>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, Eigen::Matrix<double, Rank, 1> dx, Direction dir);

private:
  Equation equation_;
  std::array<Gradient, 3> gradient_{Gradient(equation_), Gradient(equation_),
                                    Gradient(equation_)};

  Gradient boundary_scratch_{equation_};
  Gradient scratch_{equation_};
  Complete boundary_state_{equation_};
  Gradient boundary_gradient_{equation_};
  Complete interior_state_{equation_};
  Gradient interior_gradient_{equation_};
};

template <typename Equation, typename GradientT>
void HGridReconstruction2<Equation, GradientT>::
    ReconstructEmbeddedBoundaryStencil(
        span<Complete, 2> h_grid_embedded_boundary,
        span<Gradient, 2> h_grid_embedded_boundary_slopes,
        const View<const Complete>& states,
        const View<const Gradient>& gradient_x,
        const View<const Gradient>& gradient_y,
        const View<const Gradient>& gradient_z,
        const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
        Duration /*dt*/, Eigen::Matrix<double, Rank, 1> dx, Direction dir) {
  const auto d = static_cast<std::size_t>(dir);
  SetZero(interior_state_);
  SetZero(interior_gradient_);
  if (!Base::GetHGridIntegration().IntegrateInteriorCellState(
          interior_state_, interior_gradient_, states, gradient_x, gradient_y,
          gradient_z, cutcell_data, cell, dx, dir)) {
    Load(interior_state_, states, cell);
    SetZero(interior_gradient_);
  }

  SetZero(boundary_state_);
  SetZero(boundary_gradient_);
  HGridIntegrationPoints boundary_aux_data =
      GetHGridIntegrationPoints(cutcell_data, cell, dx, dir);
  if (!Base::GetHGridIntegration().IntegrateCellState(
          boundary_state_, boundary_gradient_, states, gradient_x, gradient_y,
          gradient_z, cutcell_data, boundary_aux_data, dx, dir)) {
    boundary_state_ = interior_state_;
    SetZero(boundary_gradient_);
  }
  const Coordinates<Rank> normal = GetBoundaryNormal(cutcell_data, cell);
  Reflect(boundary_gradient_, boundary_gradient_, normal, equation_);
  // boundary_gradient_ = interior_gradient_;
  StateFromComplete(equation_, scratch_, interior_state_);
  const Index<Rank> face_L = cell;
  const Index<Rank> face_R = Shift(face_L, dir, 1);
  const double betaL = cutcell_data.face_fractions[d](face_L);
  const double betaR = cutcell_data.face_fractions[d](face_R);

  const int sign = betaL < betaR ? -1 : 1;

  ForEachComponent([&](double& v, double u,
                       double du_dx) { v = u + sign * 0.5 * dx[d] * du_dx; },
                   boundary_scratch_, scratch_, interior_gradient_);

  Reflect(boundary_scratch_, boundary_scratch_, normal, equation_);

  ForEachComponent([&](double& v, double u,
                       double du_dx) { v = u + sign * 0.5 * dx[d] * du_dx; },
                   scratch_, boundary_scratch_, boundary_gradient_);

  CompleteFromState(equation_, boundary_state_, scratch_);

  if (betaL < betaR) {
    h_grid_embedded_boundary[0] = boundary_state_;
    h_grid_embedded_boundary_slopes[0] = boundary_gradient_;
    h_grid_embedded_boundary[1] = interior_state_;
    h_grid_embedded_boundary_slopes[1] = interior_gradient_;

  } else if (betaR < betaL) {
    h_grid_embedded_boundary[1] = boundary_state_;
    h_grid_embedded_boundary_slopes[1] = boundary_gradient_;
    h_grid_embedded_boundary[0] = interior_state_;
    h_grid_embedded_boundary_slopes[0] = interior_gradient_;
  }
}

} // namespace fub

#endif