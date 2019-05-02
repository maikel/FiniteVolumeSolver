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

#ifndef FUB_HYPERBOLIC_SPLIT_CUTCELL_PATCH_INTEGRATOR_HPP
#define FUB_HYPERBOLIC_SPLIT_CUTCELL_PATCH_INTEGRATOR_HPP

#include "fub/CutCellData.hpp"
#include "fub/HyperbolicSplitPatchIntegrator.hpp"

namespace fub {

template <typename Eq> class HyperbolicSplitCutCellPatchIntegrator {
public:
  static constexpr int Rank = Eq::Rank();
  using Equation = Eq;
  using Conservative = ::fub::Conservative<Eq>;

  HyperbolicSplitCutCellPatchIntegrator(const Eq& equation)
      : regular_integrator_(equation) {}

  const Eq& GetEquation() const noexcept {
    return regular_integrator_.GetEquation();
  }

  void UpdateConservatively(const View<Conservative>& next,
                            const View<const Conservative>& prev,
                            const View<const Conservative>& fluxes,
                            Direction dir, Duration dt, double dx) {
    regular_integrator_.UpdateConservatively(next, fluxes, prev, dir, dt, dx);
  }

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
        Store(nexts, next_, cell);
      } else if (flags(cell).isSingleValued() && alpha > 0.0) {
        const Index left_face = cell;
        const Index right_face = Shift(cell, dir, 1);
        Load(prev_, prevs, cell);
        Load(flux_boundary_, fluxes_boundary, cell);
        const double beta_left = face(left_face);
        const double beta_right = face(right_face);
        if (beta_left == beta_right) {
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
  HyperbolicSplitPatchIntegrator<Eq> regular_integrator_;
  Conservative prev_{GetEquation()};
  Conservative next_{GetEquation()};
  Conservative flux_shielded_left_{GetEquation()};
  Conservative flux_shielded_right_{GetEquation()};
  Conservative flux_left_{GetEquation()};
  Conservative flux_right_{GetEquation()};
  Conservative flux_boundary_{GetEquation()};
};

} // namespace fub

#endif
