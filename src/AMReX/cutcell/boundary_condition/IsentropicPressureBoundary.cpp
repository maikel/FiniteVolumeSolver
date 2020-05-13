// Copyright (c) 2018 Maikel Nadolski
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

#include "fub/AMReX/cutcell/boundary_condition/IsentropicPressureBoundary.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AnyBoundaryCondition.hpp"
#include "fub/ForEach.hpp"

#include "fub/ext/Log.hpp"

namespace fub::amrex::cutcell {
IsentropicPressureBoundaryOptions::IsentropicPressureBoundaryOptions(
    const ProgramOptions& options) {
  channel_name = GetOptionOr(options, "channel_name", channel_name);
  coarse_inner_box = GetOptionOr(options, "coarse_inner_box", coarse_inner_box);
  outer_pressure = GetOptionOr(options, "outer_pressure", outer_pressure);
  direction = GetOptionOr(options, "direction", direction);
  side = GetOptionOr(options, "side", side);
}

void IsentropicPressureBoundaryOptions::Print(SeverityLogger& log) const {
  BOOST_LOG(log) << fmt::format(" - channel_name = '{}'", channel_name);
  std::array<int, AMREX_SPACEDIM> lower, upper;
  std::copy_n(coarse_inner_box.smallEnd().getVect(), AMREX_SPACEDIM,
              lower.data());
  std::copy_n(coarse_inner_box.bigEnd().getVect(), AMREX_SPACEDIM,
              upper.data());
  BOOST_LOG(log) << fmt::format(" - coarse_inner_box = {{{{{}}}, {{{}}}}} [-]",
                                lower, upper);
  BOOST_LOG(log) << fmt::format(" - outer_pressure = {} [Pa]", outer_pressure);
  BOOST_LOG(log) << fmt::format(" - direction = {} [-]",
                                static_cast<int>(direction));
  BOOST_LOG(log) << fmt::format(" - side = {} [-]", side);
}

namespace {
int Sign(double x) { return (x > 0) - (x < 0); }

void IsentropicExpansionWithoutDissipation_(
    IdealGasMix<AMREX_SPACEDIM>& eq,
    Complete<IdealGasMix<AMREX_SPACEDIM>>& dest,
    const Complete<IdealGasMix<AMREX_SPACEDIM>>& src, double dest_pressure,
    double efficiency = 1.0) {
  Array<double, AMREX_SPACEDIM, 1> old_velocity = src.momentum / src.density;
  eq.SetReactorStateFromComplete(src);
  eq.CompleteFromReactor(dest);
  // Ekin = 0 here!
  const double h_before =
      dest.energy / dest.density + dest.pressure / dest.density;
  eq.GetReactor().SetPressureIsentropic(dest_pressure);
  eq.CompleteFromReactor(dest);
  const double h_after =
      dest.energy / dest.density + dest.pressure / dest.density;
  const double enthalpyDifference = h_before - h_after;
  const double u_border =
      Sign(enthalpyDifference) *
      std::sqrt(efficiency * std::abs(enthalpyDifference) * 2 +
                old_velocity[0] * old_velocity[0]);
  dest.momentum[0] = dest.density * u_border;
  dest.momentum[1] = dest.density * old_velocity[1];
  dest.momentum[2] = dest.density * old_velocity[2];
  dest.energy += 0.5 * dest.momentum.matrix().squaredNorm() / dest.density;
}

} // namespace

IsentropicPressureBoundary::IsentropicPressureBoundary(
    const IdealGasMix<AMREX_SPACEDIM>& eq,
    const IsentropicPressureBoundaryOptions& options)
    : equation_(eq), options_{options} {}

void IsentropicPressureBoundary::FillBoundary(::amrex::MultiFab& mf,
                                              const GriddingAlgorithm& grid,
                                              int level) {
  Complete<IdealGasMix<AMREX_SPACEDIM>> state(equation_);
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  FillBoundary(mf, alphas, geom, state);
}

void IsentropicPressureBoundary::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
    const ::amrex::Geometry& geom,
    const Complete<IdealGasMix<AMREX_SPACEDIM>>& /*state*/) {
  const int ngrow = mf.nGrow(int(options_.direction));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }
  Complete<IdealGasMix<AMREX_SPACEDIM>> state_(equation_);
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    for (const ::amrex::Box& boundary : boundaries) {
      ::amrex::Box shifted = ::amrex::shift(boundary, int(options_.direction),
                                            GetSign(options_.side) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      if (!box_to_fill.isEmpty()) {
        const int dir_v = static_cast<int>(options_.direction);
        const int x0 = options_.side ? box_to_fill.smallEnd(dir_v) - 1
                                     : box_to_fill.bigEnd(dir_v) + 1;
        auto states = MakeView<Complete<IdealGasMix<AMREX_SPACEDIM>>>(
            fab, equation_, mfi.growntilebox());
        ForEachIndex(box_to_fill, [this, &alpha, &state_, &states, x0,
                                   dir_v](auto... is) {
          std::array<std::ptrdiff_t, AMREX_SPACEDIM> dest{int(is)...};
          std::array<std::ptrdiff_t, AMREX_SPACEDIM> src = dest;
          src[static_cast<std::size_t>(dir_v)] = x0;
          ::amrex::IntVect dest_iv{int(is)...};
          ::amrex::IntVect src_iv{
              AMREX_D_DECL(int(src[0]), int(src[1]), int(src[2]))};
          if (alpha(dest_iv) > 0.0 && alpha(src_iv)) {
            Load(state_, states, src);
            const double c = state_.speed_of_sound;
            const double c2 = c * c;
            const double u2 =
                (state_.momentum / state_.density).matrix().squaredNorm();
            if (u2 < c2) {
              equation_.GetReactor().SetDensity(state_.density);
              equation_.GetReactor().SetMassFractions(state_.species);
              equation_.GetReactor().SetTemperature(state_.temperature);
              equation_.GetReactor().SetPressure(options_.outer_pressure);
              equation_.CompleteFromReactor(state_);
              IsentropicExpansionWithoutDissipation_(equation_, state_, state_,
                                                     options_.outer_pressure);
            }
            Store(states, state_, dest);
          }
        });
      }
    }
  });
}

namespace perfect_gas {

Complete<PerfectGas<AMREX_SPACEDIM>> IsentropicExpansionWithoutDissipation_(
    const PerfectGas<AMREX_SPACEDIM>& equation,
    const Complete<PerfectGas<AMREX_SPACEDIM>>& state, double outer_pressure,
    double efficiency) {
  // p_1 / p_2 = (T_1 / T_2)^(gamma / (gamma - 1))
  FUB_ASSERT(state.density > 0.0 && state.pressure > 0.0);
  FUB_ASSERT(equation.Rspec > 0.0);
  const double gamma = equation.gamma;
  const double R = equation.Rspec;
  const double gamma_minus_1_inv = equation.gamma_minus_1_inv;
  const double temperature_old = state.pressure / state.density / R;
  const double pressure_new_over_old = outer_pressure / state.pressure;
  const double exponent = (gamma - 1.0) / gamma;
  const double temperature_new =
      std::pow(pressure_new_over_old, exponent) * temperature_old;
  const double density_new = outer_pressure / temperature_new / R;
  const double h_old =
      state.pressure * (1.0 + gamma_minus_1_inv) / state.density;
  const double h_new =
      (outer_pressure * (1.0 + gamma_minus_1_inv)) / density_new;
  const double enthalpyDifference = h_old - h_new;
  const Array<double, AMREX_SPACEDIM, 1> u_old = state.momentum / state.density;
  const double u_old_norm = u_old.matrix().squaredNorm();
  const double u_border =
      Sign(enthalpyDifference) *
      std::sqrt(efficiency * std::abs(enthalpyDifference) * 2 + u_old_norm);
  const Array<double, AMREX_SPACEDIM, 1> u_new = u_border / u_old_norm * u_old;
  return equation.CompleteFromPrim(density_new, u_new, outer_pressure);
}

IsentropicPressureExpansionBoundary::IsentropicPressureExpansionBoundary(
    const PerfectGas<AMREX_SPACEDIM>& eq,
    const IsentropicPressureBoundaryOptions& options)
    : equation_(eq), options_{options} {}

void IsentropicPressureExpansionBoundary::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& grid, int level) {
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  const int ngrow = mf.nGrow(int(options_.direction));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }
  const int dir_v = static_cast<int>(options_.direction);
  Complete<PerfectGas<AMREX_SPACEDIM>> state{equation_};
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    for (const ::amrex::Box& boundary : boundaries) {
      ::amrex::Box shifted =
          ::amrex::shift(boundary, dir_v, GetSign(options_.side) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      const int x0 = options_.side ? box_to_fill.smallEnd(dir_v) - 1
                                   : box_to_fill.bigEnd(dir_v) + 1;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<PerfectGas<AMREX_SPACEDIM>>>(
            fab, equation_, mfi.growntilebox());
        ForEachIndex(box_to_fill, [this, x0, dir_v, &alpha, &state,
                                   &states](auto... is) {
          std::array<std::ptrdiff_t, AMREX_SPACEDIM> dest{int(is)...};
          std::array<std::ptrdiff_t, AMREX_SPACEDIM> src = dest;
          src[static_cast<std::size_t>(dir_v)] = x0;
          Load(state, states, src);
          Complete<PerfectGas<AMREX_SPACEDIM>> state =
              IsentropicExpansionWithoutDissipation_(
                  equation_, state, options_.outer_pressure, 1.0);
          ::amrex::IntVect iv{int(is)...};
          if (alpha(iv) > 0.0) {
            Store(states, state, dest);
          }
        });
      }
    }
  });
}

} // namespace perfect_gas

} // namespace fub::amrex::cutcell
