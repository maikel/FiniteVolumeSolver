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

template <typename FM> class KbnCutCellMethod : public FM {
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

  KbnCutCellMethod(const FM& fm) : FM(fm) {}

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
        ExactRiemannSolver<Equation> riemann_solver{equation};
        riemann_solver.SolveRiemannProblem(solution_, reflected_, state_,
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
    ExactRiemannSolver<Equation> riemann_solver{equation};
    riemann_solver.SolveRiemannProblem(solution_, reflected_, state,
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
};

} // namespace fub

#endif
