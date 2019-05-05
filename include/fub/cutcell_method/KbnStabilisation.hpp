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

#ifndef FUB_CUTCELL_METHOD_KBN_STABILISATION_HPP
#define FUB_CUTCELL_METHOD_KBN_STABILISATION_HPP

#include "fub/CutCellData.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"
#include "fub/core/span.hpp"

#include <algorithm>

namespace fub {

template <typename FM,
          typename RiemannSolver = ExactRiemannSolver<typename FM::Equation>>
class KbnCutCellMethod : public FM {
public:
  using Equation = typename FM::Equation;
  using Conservative = ::fub::Conservative<Equation>;
  using Complete = ::fub::Complete<Equation>;
  static constexpr int Rank = Equation::Rank();
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  static constexpr int StencilWidth = FM::GetStencilWidth();
  static_assert(StencilWidth > 0);
  static constexpr std::size_t StencilSize =
      static_cast<std::size_t>(2 * StencilWidth);

  KbnCutCellMethod(const FM& fm) : FM(fm) {
    stencil_.fill(Complete(FM::GetEquation()));
  }
  KbnCutCellMethod(const FM& fm, const RiemannSolver& rs)
      : FM(fm), riemann_solver_(rs) {
    stencil_.fill(Complete(FM::GetEquation()));
  }

  /// This function computes a reference state for each cut-cell.
  ///
  /// These states will be used to compute the embedded boundary fluxes and
  /// need to be fixed for a whole split cycle.
  ///
  /// \param[out] references the destination view of references states which
  ///                        will be filled by this method.
  ///
  /// \param[in] states
  void PreAdvanceHierarchy(const View<Complete>& references,
                           const View<const Complete>& states,
                           const CutCellData<Rank>& cutcell_data) {
    const Equation& equation = FM::GetEquation();
    const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(Direction::X);
    ForEachIndex(Box<0>(states), [&](auto... is) {
      std::array<std::ptrdiff_t, sRank> cell{is...};
      if (cutcell_data.flags(cell).isSingleValued()) {
        Load(state_, states, cell);
        const Eigen::Matrix<double, Rank, 1> normal =
            GetBoundaryNormal(cutcell_data, cell);
        Rotate(state_, state_, MakeRotation(normal, unit), equation);
        Reflect(reflected_, state_, unit, equation);
        riemann_solver_.SolveRiemannProblem(solution_, reflected_, state_,
                                            Direction::X);
        Rotate(solution_, solution_, MakeRotation(unit, normal), equation);
        Store(references, solution_, cell);
      }
    });
  }

  void
  ComputeBoundaryFlux(Conservative& flux, Complete& state,
                      const Complete& reference_state,
                      const Eigen::Matrix<double, Rank, 1>& boundary_normal,
                      Direction dir, Duration /* dt */, double /* dx */) {
    const Equation& equation = FM::GetEquation();
    const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(Direction::X);

    // Rotate states such that the boundary is left and the state is right
    // Reflect state in the split direction
    Rotate(state, state, MakeRotation(boundary_normal, unit), equation);
    Reflect(reflected_, state, unit, equation);
    riemann_solver_.SolveRiemannProblem(solution_, reflected_, state,
                                        Direction::X);
    Rotate(solution_, solution_, MakeRotation(unit, boundary_normal), equation);

    const int d = static_cast<int>(dir);
    const double u_advective =
        reference_state.momentum[d] / reference_state.density;
    const double u_solution = solution_.momentum[d] / solution_.density;

    equation.Flux(flux, reference_state, dir);
    flux.momentum[d] += solution_.pressure - reference_state.pressure;
    flux.energy += u_solution * solution_.pressure -
                   u_advective * reference_state.pressure;
  }

  void ComputeBoundaryFluxes(const View<Conservative>& boundary_fluxes,
                             const View<const Complete>& states,
                             const View<const Complete>& reference_states,
                             const CutCellData<Rank>& cutcell_data,
                             Direction dir, Duration dt, double dx) {
    FUB_ASSERT(Extents<0>(boundary_fluxes) == Extents<0>(states));
    const auto& flags = cutcell_data.flags;
    ForEachIndex(Box<0>(reference_states), [&](auto... is) {
      if (flags(is...).isSingleValued()) {
        // Get the state and the boundary normal in this cell.
        std::array<std::ptrdiff_t, sRank> cell{is...};
        Load(state_, states, cell);
        Load(reference_state_, reference_states, cell);
        const Eigen::Matrix<double, Rank, 1> normal =
            GetBoundaryNormal(cutcell_data, cell);

        this->ComputeBoundaryFlux(boundary_flux_left_, state_, reference_state_,
                                  normal, dir, dt, dx);

        // Store the result in our array
        Store(boundary_fluxes, boundary_flux_left_, cell);
      }
    });
  }

  void ReflectCoveredStates(
      const std::array<int, StencilSize>& is_covered,
      const std::array<std::ptrdiff_t, sRank> /* leftmost_cell */,
      const CutCellData<Rank>& /* cutcell_data */, Direction dir) {
    const Equation& equation = FM::GetEquation();
    if constexpr (StencilWidth == 1) {
      if (is_covered[0] == 2) {
        const Eigen::Matrix<double, Rank, 1> normal = UnitVector<Rank>(dir);
        Reflect(stencil_[0], stencil_[1], normal, equation);
      } else if (is_covered[1] == 2) {
        const Eigen::Matrix<double, Rank, 1> normal = -UnitVector<Rank>(dir);
        Reflect(stencil_[0], stencil_[1], normal, equation);
      }
    } else if constexpr (StencilWidth == 2) {
      // Do cases for | X | R | X | X |
      if (is_covered[1] == 0) {
        // Check for
        // | B | R | X | X | or
        // | R | R | X | X |
        if (is_covered[0] == 2) {
          const Eigen::Matrix<double, Rank, 1> normal = UnitVector<Rank>(dir);
          Reflect(stencil_[0], stencil_[1], normal, equation);
        }
        // Check for
        // | X | R | B | (R|SV) | or
        // | X | R | B | B |      or
        // | X | R | SV | B |
        // | X | R| R | B |
        if (is_covered[2] == 2) {
          const Eigen::Matrix<double, Rank, 1> normal = -UnitVector<Rank>(dir);
          Reflect(stencil_[2], stencil_[1], normal, equation);
          if (is_covered[3] == 2) {
            Reflect(stencil_[3], stencil_[0], normal, equation);
          }
        } else if (is_covered[2] == 1 && is_covered[3] == 2) {
          //          const Eigen::Matrix<double, Rank, 1> normal =
          //          GetBoundaryNormal(cutcell_data, Shift(leftmost_cell, dir,
          //          2)); Reflect(stencil_[3], stencil_[2], normal, equation);
          stencil_[3] = stencil_[2];
        } else if (is_covered[2] == 0 && is_covered[3] == 2) {
          const Eigen::Matrix<double, Rank, 1> normal = -UnitVector<Rank>(dir);
          Reflect(stencil_[3], stencil_[2], normal, equation);
        }
      } else {
        // Do cases for | X | X | R | X |
        FUB_ASSERT(is_covered[2] == 0);
        // Check for
        // | X | X | R | B | or
        // | X | X | R | (R|SV) |
        if (is_covered[3] == 2) {
          const Eigen::Matrix<double, Rank, 1> normal = -UnitVector<Rank>(dir);
          Reflect(stencil_[3], stencil_[2], normal, equation);
        }

        if (is_covered[1] == 2) {
          const Eigen::Matrix<double, Rank, 1> normal = UnitVector<Rank>(dir);
          Reflect(stencil_[1], stencil_[2], normal, equation);
          if (is_covered[0] == 2) {
            Reflect(stencil_[0], stencil_[3], normal, equation);
          }
        } else if (is_covered[0] == 2 && is_covered[1] == 1) {
          //          const Eigen::Matrix<double, Rank, 1> normal =
          //          GetBoundaryNormal(cutcell_data, Shift(leftmost_cell, dir,
          //          1)); Reflect(stencil_[0], stencil_[1], normal, equation);
          stencil_[0] = stencil_[1];
        }
        //        else if (is_covered[1] == 0 && is_covered[0] == 1) {
        //          const Eigen::Matrix<double, Rank, 1> normal =
        //          UnitVector<Rank>(dir); Reflect(stencil_[0], stencil_[1],
        //          normal, equation);
        //        }
      }
    } else {
      FUB_ASSERT(false);
      //      static_assert(false, "ReflectCoveredStates is not implemented for
      //      this stencil width");
    }
  }

  using FM::ComputeStableDt;

  /// \todo compute stable dt inside of cutcells, i.e. in the reflection with
  /// their boundary state.
  double ComputeStableDt(const View<const Complete>& states,
                         const CutCellData<Rank>& cutcell_data, Direction dir,
                         double dx) {
    double min_dt = std::numeric_limits<double>::infinity();
    auto&& flags = cutcell_data.flags;
    ForEachIndex(Shrink(Box<0>(states), dir, {0, StencilSize}),
                 [&](const auto... is) {
                   using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
                   const Index cell0{is...};
                   std::array<int, StencilSize> is_covered{};
                   for (std::size_t i = 0; i < stencil_.size(); ++i) {
                     const Index cell = Shift(cell0, dir, static_cast<int>(i));
                     Load(stencil_[i], states, cell);
                     if (flags(cell).isRegular()) {
                       is_covered[i] = 0;
                     } else if (flags(cell).isSingleValued()) {
                       is_covered[i] = 1;
                     } else {
                       is_covered[i] = 2;
                     }
                   }
                   const std::size_t iL = StencilWidth - 1;
                   const std::size_t iR = StencilWidth;
                   if (is_covered[iL] == 0 || is_covered[iR] == 0) {
                     ReflectCoveredStates(is_covered, cell0, cutcell_data, dir);
                     double dt = FM::ComputeStableDt(stencil_, dx, dir);
                     min_dt = std::min(dt, min_dt);
                   }
                   if (is_covered[iL] == 1 && is_covered[iR] == 1) {
                     for (std::size_t i = 0; i < iL; ++i) {
                       stencil_[i] = stencil_[iL];
                     }
                     for (std::size_t i = iR + 1; i < StencilSize; ++i) {
                       stencil_[i] = stencil_[iR];
                     }
                     double dt = FM::ComputeStableDt(stencil_, dx, dir);
                     min_dt = std::min(dt, min_dt);
                   }
                 });
    return min_dt;
  }

  void ComputeCutCellFluxes(const View<Conservative>& stabilised_fluxes,
                            const View<Conservative>& regular_fluxes,
                            const View<Conservative>& shielded_left_fluxes,
                            const View<Conservative>& shielded_right_fluxes,
                            const View<Conservative>& doubly_shielded_fluxes,
                            const View<const Conservative>& boundary_fluxes,
                            const View<const Complete>& states,
                            const CutCellData<Rank>& cutcell_data,
                            Direction dir, Duration dt, double dx) {
    const std::size_t dir_v = static_cast<std::size_t>(dir);
    using Index = std::array<std::ptrdiff_t, sRank>;
    const auto& flags = cutcell_data.flags;
    ForEachIndex(Box<0>(regular_fluxes), [&](auto... is) {
      const Index face{is...};
      const Index leftmost_cell = Shift(face, dir, -StencilWidth);
      std::array<int, StencilSize> is_covered{};
      for (std::size_t i = 0; i < stencil_.size(); ++i) {
        const Index cell = Shift(leftmost_cell, dir, static_cast<int>(i));
        Load(stencil_[i], states, cell);
        if (flags(cell).isRegular()) {
          is_covered[i] = 0;
        } else if (flags(cell).isSingleValued()) {
          is_covered[i] = 1;
        } else {
          is_covered[i] = 2;
        }
      }
      const std::size_t iL = StencilWidth - 1;
      const std::size_t iR = StencilWidth;
      if (is_covered[iL] == 0 || is_covered[iR] == 0) {
        ReflectCoveredStates(is_covered, leftmost_cell, cutcell_data, dir);
        FM::ComputeNumericFlux(regular_flux_, stencil_, dt, dx, dir);
        Store(regular_fluxes, regular_flux_, face);
      } else if (is_covered[iL] == 1 && is_covered[iR] == 1) {
        for (std::size_t i = 0; i < iL; ++i) {
          if (is_covered[i] == 2) {
            stencil_[i] = stencil_[iL];
          }
        }
        for (std::size_t i = iR + 1; i < StencilSize; ++i) {
          if (is_covered[i] == 2) {
            stencil_[i] = stencil_[iR];
          }
        }
        FM::ComputeNumericFlux(regular_flux_, stencil_, dt, dx, dir);
        Store(regular_fluxes, regular_flux_, face);
      }
    });

    auto boundary_centeroid = [dir_v, &cutcell_data](const Index& cell) {
      const double center = std::apply(
          [&](auto... is) {
            return cutcell_data.boundary_centeroids(is..., dir_v);
          },
          cell);
      return std::clamp(center, -0.5, +0.5);
    };

    ForEachIndex(Box<0>(stabilised_fluxes), [&](auto... is) {
      const Index face = Index{is...};
      const double beta_unshielded = cutcell_data.unshielded_fractions(face);
      const double beta_left = cutcell_data.shielded_left_fractions(face);
      const double beta_right = cutcell_data.shielded_right_fractions(face);
      const double beta_ds = cutcell_data.doubly_shielded_fractions(face);
      FUB_ASSERT(beta_unshielded >= 0 && beta_left >= 0 && beta_right >= 0 &&
                 beta_ds >= 0);

      const Index cell_right{is...};
      const Index cell_left = Shift(cell_right, dir, -1);
      if (flags(cell_left).isSingleValued() ||
          flags(cell_right).isSingleValued()) {
        // The unshielded face fraction does not need any stabilisation.
        // Take the regular flux here!
        Load(regular_flux_, regular_fluxes, face);

        // Compute stabilisation for shielded face fraction from left
        if (beta_left > 0.0) {
          FUB_ASSERT(flags(cell_left).isSingleValued());
          Load(boundary_flux_left_, boundary_fluxes, cell_left);
          const double center = boundary_centeroid(cell_left);
          // d is the average distance to the boundary
          const double d = std::clamp(0.5 - center, 0.0, 1.0);
          ForEachComponent<Conservative>(
              [=](double& f_shielded, double f_regular, double f_B) {
                f_shielded = d * f_regular + (1.0 - d) * f_B;
              },
              shielded_left_flux_, regular_flux_, boundary_flux_left_);
        } else {
          ForEachComponent<Conservative>([](double& x) { x = 0.0; },
                                         shielded_left_flux_);
        }

        // Compute stabilisation for shielded face fraction from right
        if (beta_right > 0.0) {
          FUB_ASSERT(flags(cell_right).isSingleValued());
          Load(boundary_flux_right_, boundary_fluxes, cell_right);
          const double center = boundary_centeroid(cell_right);
          const double d = std::clamp(0.5 + center, 0.0, 1.0);
          ForEachComponent<Conservative>(
              [d](double& f_shielded, double f_regular, double f_B) {
                f_shielded = (d * f_regular + (1.0 - d) * f_B);
              },
              shielded_right_flux_, regular_flux_, boundary_flux_right_);
        } else {
          ForEachComponent<Conservative>([](double& x) { x = 0.0; },
                                         shielded_right_flux_);
        }

        if (beta_ds > 0.0) {
          FUB_ASSERT(flags(cell_left).isSingleValued() &&
                     flags(cell_right).isSingleValued());
          Load(boundary_flux_left_, boundary_fluxes, cell_left);
          Load(boundary_flux_right_, boundary_fluxes, cell_right);
          Load(stencil_[0], states, cell_left);
          Load(stencil_[1], states, cell_right);
          const double dL = 0.5 - boundary_centeroid(cell_left);
          const double dR = 0.5 + boundary_centeroid(cell_right);
          FUB_ASSERT(0.0 < dL && dL < 1.0);
          FUB_ASSERT(0.0 < dR && dR < 1.0);
          const double one_over_ds = 1.0 / (dL + dR);
          const double dLdR_over_ds = dL * dR * one_over_ds;
          const double dL_over_ds = dL * one_over_ds;
          const double dR_over_ds = dR * one_over_ds;
          ForEachComponent<Conservative>(
              [=](double& f_ds, double uL, double uR, double fL_B,
                  double fR_B) {
                f_ds = (dLdR_over_ds * (uL - uR) + dL_over_ds * fR_B +
                        dR_over_ds * fL_B);
              },
              doubly_shielded_flux_, stencil_[0], stencil_[1],
              boundary_flux_left_, boundary_flux_right_);
        } else {
          ForEachComponent<Conservative>([](double& x) { x = 0.0; },
                                         doubly_shielded_flux_);
        }

        // Now assemble the total stabilised flux as an area weighted sum of
        // unshielded, shielded and doubly shielded parts.
        const double beta_total =
            beta_unshielded + beta_left + beta_right + beta_ds;
        if (beta_total > 0.0) {
          const double beta_us_frac = beta_unshielded / beta_total;
          const double beta_left_frac = beta_left / beta_total;
          const double beta_right_frac = beta_right / beta_total;
          ForEachComponent<Conservative>(
              [=](double& stabilized_flux, double regular_flux,
                  double shielded_left_flux, double shielded_right_flux) {
                stabilized_flux = beta_us_frac * regular_flux +
                                  beta_left_frac * shielded_left_flux +
                                  beta_right_frac * shielded_right_flux;
              },
              cutcell_flux_, regular_flux_, shielded_left_flux_,
              shielded_right_flux_);
        } else {
          ForEachComponent<Conservative>([](double& x) { x = 0.0; },
                                         cutcell_flux_);
        }

        Store(shielded_left_fluxes, shielded_left_flux_, face);
        Store(shielded_right_fluxes, shielded_right_flux_, face);
        Store(doubly_shielded_fluxes, doubly_shielded_flux_, face);
        Store(stabilised_fluxes, cutcell_flux_, face);
      } else {
        Load(regular_flux_, regular_fluxes, face);
        ForEachComponent<Conservative>([](double& x) { x = 0.0; },
                                       cutcell_flux_);
        Store(shielded_left_fluxes, cutcell_flux_, face);
        Store(shielded_right_fluxes, cutcell_flux_, face);
        Store(doubly_shielded_fluxes, cutcell_flux_, face);
        Store(stabilised_fluxes, regular_flux_, face);
      }
    });
  }

  void UpdateConservatively(const View<Conservative>& next,
                            const View<const Conservative>& prev,
                            const View<const Conservative>& fluxes,
                            Direction dir, Duration dt, double dx) {
    regular_integrator_.UpdateConservatively(next, fluxes, prev, dir, dt, dx);
  }

  /// \brief This is a specialized version of a conservative time update
  /// tailored to this particular KBN stabilisation scheme.
  ///
  /// It is less general than HyperbolicSplitCutCellPatchIntegrator but should
  /// be more stable.
  ///
  /// \param[out] nexts The state array for the next time level.
  ///
  /// \param[in] prevs The state array for the current time level.
  ///
  /// \param[in] stabilised fluxes An array for stabilised fluxes, which are
  /// convex combinations of regular and shielded ones.
  void UpdateConservatively(
      const View<Conservative>& nexts, const View<const Conservative>& prevs,
      const View<const Conservative>& stabilised_fluxes,
      const View<const Conservative>& regular_fluxes,
      const View<const Conservative>& shielded_left_fluxes,
      const View<const Conservative>& shielded_right_fluxes,
      const View<const Conservative>& /* doubly_shielded_fluxes */,
      const View<const Conservative>& fluxes_boundary,
      const CutCellData<Rank>& cutcell_data, Direction dir, Duration dt,
      double dx) {
    using Index = std::array<std::ptrdiff_t, static_cast<std::size_t>(Rank)>;
    const double dt_over_dx = dt.count() / dx;
    FUB_ASSERT(Extents<0>(nexts) == Extents<0>(prevs));
    auto volume = cutcell_data.volume_fractions;
    auto face = cutcell_data.face_fractions;
    auto flags = cutcell_data.flags;
    auto center = cutcell_data.boundary_centeroids;
    const std::ptrdiff_t dir_v = static_cast<std::ptrdiff_t>(dir);
    ForEachIndex(Shrink(Box<0>(regular_fluxes), dir, {0, 1}), [&](auto... is) {
      const Index cell{is...};
      double alpha = volume(cell);
      // if the cell is regular we take stabilised fluxes and do the standard
      // forward euler step.
      if (flags(cell).isRegular()) {
        const Index left_face = cell;
        const Index right_face = Shift(cell, dir, 1);
        Load(prev_, prevs, cell);
        Load(flux_left_, stabilised_fluxes, left_face);
        Load(flux_right_, stabilised_fluxes, right_face);
        ForEachComponent<Conservative>(
            [=](double& next, double prev, double f_left, double f_right) {
              next = prev + dt_over_dx * (f_left - f_right);
            },
            next_, prev_, flux_left_, flux_right_);
        FUB_ASSERT(next_.density > 0.0);
        FUB_ASSERT(next_.energy >
                   0.5 * next_.momentum.matrix().squaredNorm() / next_.density);
        Store(nexts, next_, cell);
      } else if (flags(cell).isSingleValued() && alpha > 0.0) {
        const Index left_face = cell;
        const Index right_face = Shift(cell, dir, 1);
        Load(prev_, prevs, cell);
        Load(flux_boundary_, fluxes_boundary, cell);
        const double beta_left = face(left_face);
        const double beta_right = face(right_face);
        if (beta_left == beta_right) {
          // if both sides have the same sized area fraction we do the standard
          // forward euler step.
          Load(flux_left_, stabilised_fluxes, left_face);
          Load(flux_right_, stabilised_fluxes, right_face);
          ForEachComponent<Conservative>(
              [=](double& next, double prev, double f_left, double f_right) {
                next = prev + dt_over_dx * (f_left - f_right);
              },
              next_, prev_, flux_left_, flux_right_);
        } else if (beta_left > beta_right) {
          Load(flux_left_, regular_fluxes, left_face);
          Load(flux_shielded_left_, shielded_left_fluxes, left_face);
          Load(flux_right_, stabilised_fluxes, right_face);

          const double beta_ss_left =
              cutcell_data.shielded_left_fractions(left_face);
          const double beta_ss_right =
              cutcell_data.shielded_right_fractions(left_face);
          const double beta_us = cutcell_data.unshielded_fractions(left_face);
          const double distance_to_boundary = 0.5 + center(is..., dir_v);
          FUB_ASSERT(distance_to_boundary > 0.0);
          const double beta_us_frac =
              beta_right == 0.0 ? 1.0
                                : std::clamp(beta_us / beta_right, 0.0, 1.0);
          const double beta_ssL_frac =
              beta_right == 0.0
                  ? 0.0
                  : std::clamp(beta_ss_left / beta_right, 0.0, 1.0);
          FUB_ASSERT(std::abs(1.0 - (beta_us_frac + beta_ssL_frac)) < 1e-6);

          const double betaR_over_alpha =
              std::clamp(beta_right / alpha, 0.0, 1.0);
          const double dBeta_over_alpha = std::clamp(
              (beta_ss_right * distance_to_boundary) / alpha, 0.0, 1.0);

          ForEachComponent<Conservative>(
              [=](double& next, double prev, double f_left, double f_ssL,
                  double f_right, double f_B) {
                const double dF =
                    beta_us_frac * f_left + beta_ssL_frac * f_ssL - f_right;
                const double dBoundary = f_left - f_B;
                next = prev + dt_over_dx * betaR_over_alpha * dF +
                       dt_over_dx * dBeta_over_alpha * dBoundary;
              },
              next_, prev_, flux_left_, flux_shielded_left_, flux_right_,
              flux_boundary_);

        } else if (beta_left < beta_right) {
          Load(flux_left_, stabilised_fluxes, left_face);
          Load(flux_shielded_right_, shielded_right_fluxes, right_face);
          Load(flux_right_, regular_fluxes, right_face);

          const double beta_ss_left =
              cutcell_data.shielded_left_fractions(right_face);
          const double beta_ss_right =
              cutcell_data.shielded_right_fractions(right_face);
          const double beta_us = cutcell_data.unshielded_fractions(right_face);
          const double distance_to_boundary = 0.5 - center(is..., dir_v);
          FUB_ASSERT(distance_to_boundary > 0.0);
          const double beta_us_frac =
              beta_left == 0.0 ? 1.0
                               : std::clamp(beta_us / beta_left, 0.0, 1.0);
          const double beta_ssR_frac =
              beta_left == 0.0
                  ? 0.0
                  : std::clamp(beta_ss_right / beta_left, 0.0, 1.0);
          FUB_ASSERT(std::abs(1.0 - (beta_us_frac + beta_ssR_frac)) < 1e-6);

          const double betaL_over_alpha =
              std::clamp(beta_left / alpha, 0.0, 1.0);
          const double dBeta_over_alpha = std::clamp(
              (beta_ss_left * distance_to_boundary) / alpha, 0.0, 1.0);

          ForEachComponent<Conservative>(
              [=](double& next, double prev, double f_left, double f_right,
                  double f_ssR, double f_B) {
                const double dF =
                    f_left - (beta_us_frac * f_right + beta_ssR_frac * f_ssR);
                const double dBoundary = f_B - f_right;
                next = prev + dt_over_dx * betaL_over_alpha * dF +
                       dt_over_dx * dBeta_over_alpha * dBoundary;
              },
              next_, prev_, flux_left_, flux_right_, flux_shielded_right_,
              flux_boundary_);
        }
        Store(nexts, next_, cell);
      }
    });
  }

private:
  std::array<Complete, StencilSize> stencil_{};
  Complete state_{FM::GetEquation()};
  Complete solution_{FM::GetEquation()};
  Complete reflected_{FM::GetEquation()};
  Complete reference_state_{FM::GetEquation()};
  Conservative cutcell_flux_{FM::GetEquation()};
  Conservative regular_flux_{FM::GetEquation()};
  Conservative boundary_flux_left_{FM::GetEquation()};
  Conservative boundary_flux_right_{FM::GetEquation()};
  Conservative shielded_right_flux_{FM::GetEquation()};
  Conservative shielded_left_flux_{FM::GetEquation()};
  Conservative doubly_shielded_flux_{FM::GetEquation()};

  HyperbolicSplitPatchIntegrator<Equation> regular_integrator_{
      FM::GetEquation()};
  RiemannSolver riemann_solver_{FM::GetEquation()};
  Conservative prev_{FM::GetEquation()};
  Conservative next_{FM::GetEquation()};
  Conservative flux_shielded_left_{FM::GetEquation()};
  Conservative flux_shielded_right_{FM::GetEquation()};
  Conservative flux_left_{FM::GetEquation()};
  Conservative flux_right_{FM::GetEquation()};
  Conservative flux_boundary_{FM::GetEquation()};
};

} // namespace fub

#endif
