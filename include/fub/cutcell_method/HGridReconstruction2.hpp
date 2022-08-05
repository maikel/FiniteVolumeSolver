#ifndef FUB_CUTCELL_METHOD_H_GRID_REC_2_HPP
#define FUB_CUTCELL_METHOD_H_GRID_REC_2_HPP

#include "fub/NewtonIteration.hpp"
#include "fub/cutcell_method/HGridReconstruction.hpp"

namespace fub {

template <typename Equation, typename GradientT, typename ReconstructionMethod>
class HGridReconstruction2 : private HGridReconstruction<Equation, GradientT> {
public:
  using Base = HGridReconstruction<Equation, GradientT>;
  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using Gradient = GradientT;

  static constexpr int Rank = Equation::Rank();

  explicit HGridReconstruction2(const Equation& eq, ReconstructionMethod rec)
      : Base(eq), equation_{eq}, rec_{std::move(rec)} {}

  using Base::ReconstructRegularStencil;

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
      Duration dt, const Coordinates<Rank>& dx, Direction dir);

  void ReconstructEmbeddedBoundaryStencil(
      span<Complete, 2> h_grid_embedded_boundary,
      span<Gradient, 2> h_grid_embedded_boundary_slopes,
      const View<const Complete>& states,
      const View<const Gradient>& gradient_x,
      const View<const Gradient>& gradient_y,
      const View<const Gradient>& gradient_z,
      const View<const Complete>& boundary_reference_states,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration dt, Eigen::Matrix<double, Rank, 1> dx, Direction dir);

  double FindPressure(double u_star, const Complete& right, Direction dir);

private:
  Equation equation_;
  ReconstructionMethod rec_;
  Gradient gradient_dir_{equation_};
  std::array<Gradient, 3> gradient_{Gradient(equation_), Gradient(equation_),
                                    Gradient(equation_)};

  Gradient integral_scratch_{equation_};
  Gradient boundary_scratch_{equation_};
  Gradient scratch_{equation_};
  Complete reference_state_{equation_};
  Complete boundary_rec_{equation_};
  Complete boundary_state_{equation_};
  Gradient boundary_gradient_{equation_};
  Complete interior_rec_{equation_};
  Complete interior_state_{equation_};
  Gradient interior_gradient_{equation_};
  double singly_shielded_factor_{0.0};
};

template <typename Equation, typename GradientT, typename ReconstructionMethod>
double
HGridReconstruction2<Equation, GradientT, ReconstructionMethod>::FindPressure(
    double u_star, const Complete& right, Direction dir) {
  const double g = equation_.gamma;
  const double gp1 = g + 1.0;
  const double gm1 = g - 1.0;
  const double exponent = 0.5 * gm1 / g;
  const double d_exponent = -0.5 * gp1 / g;
  const double rhoR = right.density;
  const double pR = right.pressure;
  const double aR = right.speed_of_sound;
  const double A_r = 2.0 / (gp1 * rhoR);
  const double B_r = gm1 / gp1 * pR;
  auto f_right = [=](double p) {
    if (p > pR) {
      return (p - pR) * std::sqrt(A_r / (p + B_r));
    } else {
      return 2 * aR / gm1 * (std::pow(p / pR, exponent) - 1.0);
    }
  };
  auto df_right = [=](double p) {
    if (p > pR) {
      return std::sqrt(A_r / (B_r + p)) * (1.0 - 0.5 * (p - pR) / (B_r + p));
    } else {
      return std::pow(p / pR, d_exponent) / (rhoR * aR);
    }
  };
  auto velocity = [&](const Complete& state) {
    return state.momentum[static_cast<int>(dir)] / state.density;
  };
  const double uR = velocity(right);
  const double du = uR - u_star;
  auto f = [=](double h) { return f_right(h) + du; };
  auto df = [=](double h) { return df_right(h); };
  const double pM = NewtonIteration(f, df, pR);
  return pM;
}

template <typename Equation, typename GradientT, typename ReconstructionMethod>
void HGridReconstruction2<Equation, GradientT, ReconstructionMethod>::
    ReconstructEmbeddedBoundaryStencil(
        span<Complete, 2> h_grid_embedded_boundary,
        span<Gradient, 2> h_grid_embedded_boundary_slopes,
        const View<const Complete>& states,
        const View<const Gradient>& gradient_x,
        const View<const Gradient>& gradient_y,
        const View<const Gradient>& gradient_z,
        const View<const Complete>& boundary_reference_states,
        const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
        Duration dt, Eigen::Matrix<double, Rank, 1> dx, Direction dir) {
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
  Reflect(boundary_state_, boundary_state_, normal, equation_);
  Reflect(boundary_gradient_, boundary_gradient_, normal, equation_);

  // // TODO Here we need to require a reference velocity (u,v)
  Load(reference_state_, boundary_reference_states, cell);

  const Index<Rank> face_L = cell;
  const Index<Rank> face_R = Shift(face_L, dir, 1);
  const double betaL = cutcell_data.face_fractions[d](face_L);
  const double betaR = cutcell_data.face_fractions[d](face_R);
  // const double dBeta = betaR - betaL;
  // const int sign = (dBeta > 0) - (dBeta < 0);
  singly_shielded_factor_ = 0.0;
  if (reference_state_.density > 0) {
    const double u_star = euler::Velocity(equation_, reference_state_, d);
    // mass flows out of the cut cell
    // take data from the boundary cell
    if (u_star * normal[d] >= 0) {
      const Side upper = u_star > 0 ? Side::Upper : Side::Lower;
      const Side lower = u_star > 0 ? Side::Lower : Side::Upper;
      rec_.Reconstruct(boundary_rec_, boundary_state_, boundary_gradient_, dt, dx[d], dir, upper);
      rec_.Reconstruct(interior_rec_, interior_state_, interior_gradient_, dt, dx[d], dir, lower);
      StateFromComplete(equation_, scratch_, boundary_rec_);
      const int u_star_sign = (u_star > 0) - (u_star < 0);
      interior_rec_.momentum[d] *= u_star_sign;
      const double p_star = FindPressure(std::abs(u_star), interior_rec_, dir);
      const double pL = boundary_rec_.pressure;
      const double rhoL_star = reference_state_.density;
      RequireMassflow_SolveExactRiemannProblem RequireMassflow{};
      RequireMassflow(equation_, scratch_, rhoL_star, std::abs(u_star), p_star, pL, dir);
      scratch_.velocity[d] *= u_star_sign;
      CompleteFromState(equation_, boundary_rec_, scratch_);

      rec_.ReconstructInvers(boundary_state_, boundary_rec_, boundary_gradient_, dt, dx[d], dir, upper);
    }
    // mass flows into the cut cell
    // take data from the interior cell
    else {
      const Side upper = u_star > 0 ? Side::Upper : Side::Lower;
      const Side lower = u_star > 0 ? Side::Lower : Side::Upper;
      rec_.Reconstruct(interior_rec_, interior_state_, interior_gradient_, dt, dx[d], dir, upper);
      rec_.Reconstruct(boundary_rec_, boundary_state_, boundary_gradient_, dt, dx[d], dir, lower);
      StateFromComplete(equation_, scratch_, interior_rec_);
      const int u_star_sign = (u_star > 0) - (u_star < 0);
      boundary_rec_.momentum[d] *= u_star_sign;
      const double p_star = FindPressure(std::abs(u_star), boundary_rec_, dir);
      const double pL = interior_rec_.pressure;
      const double rhoL_star = reference_state_.density;
      RequireMassflow_SolveExactRiemannProblem RequireMassflow{};
      RequireMassflow(equation_, scratch_, rhoL_star, u_star, p_star, pL, dir);
      scratch_.velocity[d] *= u_star_sign;
      singly_shielded_factor_ = rhoL_star / interior_rec_.density;
      CompleteFromState(equation_, interior_rec_, scratch_);

      rec_.ReconstructInvers(interior_state_, interior_rec_, interior_gradient_, dt, dx[d], dir, upper);
    }
  }

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

template <typename Equation, typename GradientT, typename ReconstructionMethod>
void HGridReconstruction2<Equation, GradientT, ReconstructionMethod>::
    ReconstructSinglyShieldedStencil(
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
    StateFromComplete(equation_, scratch_, h_grid_singly_shielded[1]);
    scratch_ += gradient_dir_;
    CompleteFromState(equation_, h_grid_singly_shielded[1], scratch_);

    Load(h_grid_singly_shielded_gradients[1], gradients[d], iR);

    if (singly_shielded_factor_) {
      StateFromComplete(equation_, boundary_scratch_, h_grid_eb[1]);
      // boundary_scratch_ = scratch_;
      // boundary_scratch_.density *= singly_shielded_factor_;
      ForEachComponent(
          [&](double& Q_r, double& dQ_r, double Q_b, double dQ_b, double Q_i,
              double dQ_i) {
            Q_r = (1.0 - alpha) * (Q_b + 0.5 * alpha * dx[d] * dQ_b) +
                  alpha * (Q_i + 0.5 * (1.0 - alpha) * dx[d] * dQ_i);
            dQ_r = (1.0 - alpha) * dQ_b + alpha * dQ_i;
          },
          integral_scratch_, h_grid_singly_shielded_gradients[1],
          boundary_scratch_, gradient_[d], scratch_, gradient_[d]);
      CompleteFromState(equation_, h_grid_singly_shielded[1],
                        integral_scratch_);
    }

    Load(boundary_state_, states, cell);
    Load(gradient_[0], gradient_x, cell);
    Load(gradient_[1], gradient_y, cell);
    Load(gradient_[2], gradient_z, cell);
    const Coordinates<Rank> xC =
        GetAbsoluteVolumeCentroid(cutcell_data, cell, dx);
    const Coordinates<Rank> xCssR =
        xB + 0.5 * alpha * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
    ApplyGradient(gradient_dir_, grads, xCssR - xC);
    StateFromComplete(equation_, scratch_, boundary_state_);
    scratch_ += gradient_dir_;

    StateFromComplete(equation_, boundary_scratch_, h_grid_eb[0]);
    ForEachComponent(
        [&](double& Q_l, double& dQ_l, double Q_b, double dQ_b, double Q_i,
            double dQ_i) {
          Q_l = (1 - alpha) * (Q_b + 0.5 * alpha * dx[d] * dQ_b) + alpha * Q_i;
          dQ_l = (1 - alpha) * dQ_b + alpha * dQ_i;
        },
        integral_scratch_, h_grid_singly_shielded_gradients[0],
        boundary_scratch_, h_grid_eb_gradients[0], scratch_, gradient_[d]);
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

    if (singly_shielded_factor_) {
      StateFromComplete(equation_, boundary_scratch_, h_grid_eb[0]);
      // boundary_scratch_ = scratch_;
      // boundary_scratch_.density *= singly_shielded_factor_;
      ForEachComponent(
          [&](double& Q_r, double& dQ_r, double Q_b, double dQ_b, double Q_i,
              double dQ_i) {
            Q_r = (1.0 - alpha) * (Q_b - 0.5 * alpha * dx[d] * dQ_b) +
                  alpha * (Q_i - 0.5 * (1.0 - alpha) * dx[d] * dQ_i);
            dQ_r = (1.0 - alpha) * dQ_b + alpha * dQ_i;
          },
          integral_scratch_, h_grid_singly_shielded_gradients[0],
          boundary_scratch_, gradient_[d], scratch_, gradient_[d]);
      CompleteFromState(equation_, h_grid_singly_shielded[0],
                        integral_scratch_);
    }

    Load(boundary_state_, states, cell);
    StateFromComplete(equation_, scratch_, boundary_state_);
    Load(gradient_[0], gradient_x, cell);
    Load(gradient_[1], gradient_y, cell);
    Load(gradient_[2], gradient_z, cell);
    const Coordinates<Rank> xC =
        GetAbsoluteVolumeCentroid(cutcell_data, cell, dx);
    const Coordinates<Rank> xCssR =
        xB - 0.5 * alpha * dx[d] * Eigen::Matrix<double, Rank, 1>::Unit(d);
    ApplyGradient(gradient_dir_, grads, xCssR - xC);
    scratch_ += gradient_dir_;

    StateFromComplete(equation_, boundary_scratch_, h_grid_eb[1]);
    ForEachComponent(
        [&](double& Q_r, double& dQ_r, double Q_b, double dQ_b, double Q_i,
            double dQ_i) {
          Q_r =
              (1.0 - alpha) * (Q_b - 0.5 * alpha * dx[d] * dQ_b) + alpha * Q_i;
          dQ_r = (1.0 - alpha) * dQ_b + alpha * dQ_i;
        },
        integral_scratch_, h_grid_singly_shielded_gradients[1],
        boundary_scratch_, h_grid_eb_gradients[1], scratch_, gradient_[d]);
    CompleteFromState(equation_, h_grid_singly_shielded[1], integral_scratch_);
  }
}

} // namespace fub

#endif