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
  static constexpr int StencilWidth = FM::GetStencilWidth();

  KbnCutCellMethod(const FM& fm) : FM(fm) {}

  void
  ComputeBoundaryFlux(Conservative& flux, Complete& reflected, Complete& state,
                      const Eigen::Matrix<double, Rank, 1>& boundary_normal,
                      Direction dir, Duration dt, double dx) {
    const Equation& equation = FM::GetEquation();
    const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(dir);

    // Rotate states such that the boundary is left and the state is right
    //    Rotate(state, state, MakeRotation(boundary_normal, unit), equation);

    // Reflect state in the split direction
    Reflect(reflected, state, unit, equation);

    // Fill the stencil of our numeric flux method with reflected states
    {
      Eigen::Index d = static_cast<Eigen::Index>(dir);
      const int half_n = static_cast<int>(stencil_.size()) / 2;
      const auto first = stencil_.begin();
      const auto last = stencil_.end();
      const auto half = std::next(first, half_n);
      if (boundary_normal[d] > 0.0) {
        std::fill(first, half, reflected);
        std::fill(half, last, state);
      } else {
        std::fill(first, half, state);
        std::fill(half, last, reflected);
      }
    }

    // Compute the numeric boundary flux with this stencil.
    // This will always yield an first order approximation
    FM::ComputeNumericFlux(flux, stencil_, dt, dx, dir);

    // Rotate the computed numeric flux into the normal direction.
    //    Rotate(flux, flux, MakeRotation(unit, boundary_normal), equation);
  }

  template <typename LB, typename LS>
  void ComputeBoundaryFluxes(const View<Conservative, LB>& boundary_fluxes,
                             const View<const Complete, LS>& states,
                             const CutCellData<Rank>& cutcell_data,
                             Direction dir, Duration dt, double dx) {
    FUB_ASSERT(Extents<0>(boundary_fluxes) == Extents<0>(states));
    const auto& flags = cutcell_data.flags;
    ForEachIndex(Box<0>(states), [&](auto... is) {
      if (flags(is...).isSingleValued()) {
        // Get the state and the boundary normal in this cell.
        std::array<std::ptrdiff_t, Rank> cell{is...};
        Load(state_, states, cell);
        const Eigen::Matrix<double, Rank, 1> normal =
            GetBoundaryNormal(cutcell_data, cell);
        // TODO: assertion: almost equal to 1.0
        FUB_ASSERT(normal.squaredNorm() > 0.0);

        this->ComputeBoundaryFlux(boundary_flux_left_, reflected_, state_,
                                  normal, dir, dt, dx);

        // Store the result in our array
        Store(boundary_fluxes, boundary_flux_left_, cell);
      }
    });
  }

  bool ReflectCoveredStates(const std::array<int, 2 * StencilWidth>& is_covered,
                            Direction dir) {
    const Equation& equation = FM::GetEquation();
    constexpr int mid_L = StencilWidth - 1;
    constexpr int mid_R = StencilWidth;
    std::size_t j = static_cast<std::size_t>(mid_R);
    FUB_ASSERT(mid_L >= 0);
    int i = mid_L;
    Eigen::Matrix<double, Rank, 1> normal{};
    normal[static_cast<int>(dir)] = 1.0;
    while (i >= 0) {
      std::size_t is = static_cast<std::size_t>(i);
      std::size_t js = static_cast<std::size_t>(j);
      if (is_covered[is]) {
        if (is_covered[js]) {
          stencil_[is] = stencil_[is + 1];
        } else {
          Reflect(stencil_[is], stencil_[js], normal, equation);
        }
      }
      --i;
      ++j;
    }
    j = mid_L;
    i = mid_R;
    while (i < 2 * StencilWidth) {
      FUB_ASSERT(j >= 0);
      std::size_t is = static_cast<std::size_t>(i);
      std::size_t js = static_cast<std::size_t>(j);
      if (is_covered[is]) {
        if (is_covered[js]) {
          stencil_[is] = stencil_[is - 1];
        } else {
          Reflect(stencil_[is], stencil_[js], normal, equation);
        }
      }
      ++i;
      --j;
    }
    return true;
  }

  using FM::ComputeStableDt;

  /// \todo compute stable dt inside of cutcells, i.e. in the reflection with
  /// their boundary state.
  template <typename Layout>
  double ComputeStableDt(const View<const Complete, Layout>& states,
                         const CutCellData<Rank>& cutcell_data, Direction dir,
                         double dx) {
    double min_dt = std::numeric_limits<double>::infinity();
    auto&& flags = cutcell_data.flags;
    ForEachIndex(Shrink(Box<0>(states), dir, {0, 2 * StencilWidth}),
                 [&](const auto... is) {
                   using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
                   const Index cell0{is...};
                   std::array<int, 2 * StencilWidth> is_covered{};
                   for (std::size_t i = 0; i < stencil_.size(); ++i) {
                     const Index cell = Shift(cell0, dir, static_cast<int>(i));
                     Load(stencil_[i], states, cell);
                     is_covered[i] = !(flags(cell).isSingleValued() ||
                                       flags(cell).isRegular());
                   }
                   const std::size_t iL = StencilWidth - 1;
                   const std::size_t iR = StencilWidth;
                   if (!is_covered[iL] || !is_covered[iR]) {
                     ReflectCoveredStates(is_covered, dir);
                     double dt = FM::ComputeStableDt(stencil_, dx, dir);
                     min_dt = std::min(dt, min_dt);
                   }
                 });
    return min_dt;
  }

  template <typename L1, typename L2, typename L3, typename L4, typename L5,
            typename L6, typename L7>
  void
  ComputeCutCellFluxes(const View<Conservative, L1>& stabilised_fluxes,
                       const View<Conservative, L2>& regular_fluxes,
                       const View<Conservative, L3>& shielded_left_fluxes,
                       const View<Conservative, L4>& shielded_right_fluxes,
                       const View<Conservative, L5>& doubly_shielded_fluxes,
                       const View<const Conservative, L6>& boundary_fluxes,
                       const View<const Complete, L7>& states,
                       const CutCellData<Rank>& cutcell_data, Direction dir,
                       Duration dt, double dx) {
    const std::size_t dir_v = static_cast<std::size_t>(dir);
    using Index = std::array<std::ptrdiff_t, Rank>;
    const auto& flags = cutcell_data.flags;
    ForEachIndex(Box<0>(regular_fluxes), [&](auto... is) {
      const Index face{is...};
      const Index leftmost_cell = Shift(face, dir, -StencilWidth);
      std::array<int, 2 * StencilWidth> is_covered{};
      for (std::size_t i = 0; i < stencil_.size(); ++i) {
        const Index cell = Shift(leftmost_cell, dir, static_cast<int>(i));
        Load(stencil_[i], states, cell);
        is_covered[i] =
            !(flags(cell).isSingleValued() || flags(cell).isRegular());
      }
      const std::size_t iL = StencilWidth - 1;
      const std::size_t iR = StencilWidth;
      if (!is_covered[iL] || !is_covered[iR]) {
        ReflectCoveredStates(is_covered, dir);
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
  std::array<Complete, 2 * StencilWidth> stencil_{};
  Complete state_{FM::GetEquation()};
  Complete reflected_{FM::GetEquation()};
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
