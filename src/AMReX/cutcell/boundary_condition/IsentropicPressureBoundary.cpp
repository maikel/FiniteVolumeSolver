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
  std::copy_n(coarse_inner_box.smallEnd().getVect(), AMREX_SPACEDIM, lower.data());
  std::copy_n(coarse_inner_box.bigEnd().getVect(), AMREX_SPACEDIM, upper.data());
  BOOST_LOG(log) << fmt::format(" - coarse_inner_box = {{{{{}}}, {{{}}}}} [-]", lower, upper);
  BOOST_LOG(log) << fmt::format(" - outer_pressure = {} [Pa]", outer_pressure);
  BOOST_LOG(log) << fmt::format(" - direction = {} [-]", static_cast<int>(direction));
  BOOST_LOG(log) << fmt::format(" - side = {} [-]", side);
}

namespace {
inline int GetSign(int side) { return (side == 0) - (side == 1); }

int Sign(double x) { return (x > 0) - (x < 0); }

void IsentropicExpansionWithoutDissipation(
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
                old_velocity.matrix().squaredNorm());
  dest.momentum[0] = dest.density * u_border;
  dest.energy += 0.5 * dest.density * u_border * u_border;
}

} // namespace

IsentropicPressureBoundary::IsentropicPressureBoundary(
    const IdealGasMix<AMREX_SPACEDIM>& eq,
    const IsentropicPressureBoundaryOptions& options)
    : equation_(eq), options_{options} {}

namespace {
double TotalVolume(const PatchHierarchy& hier, int level,
                   const ::amrex::Box& box) {
  double local_volume(0.0);
  const ::amrex::EBFArrayBoxFactory& factory = *hier.GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory.getVolFrac();
  ForEachFab(alphas, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    ForEachIndex(box & mfi.tilebox(), [&](auto... is) {
      ::amrex::IntVect index{static_cast<int>(is)...};
      const double frac = alpha(index);
      local_volume += frac;
    });
  });
  double global_volume = 0.0;
  MPI_Allreduce(&local_volume, &global_volume, 1, MPI_DOUBLE, MPI_SUM,
                ::amrex::ParallelDescriptor::Communicator());
  return global_volume;
}

void AverageState(Complete<IdealGasMix<AMREX_SPACEDIM>>& state,
                  const PatchHierarchy& hier, int level,
                  const ::amrex::Box& box) {
  const double total_volume = TotalVolume(hier, level, box);
  const int ncomp = hier.GetDataDescription().n_state_components;
  std::vector<double> state_buffer(static_cast<std::size_t>(ncomp));
  const ::amrex::EBFArrayBoxFactory& factory = *hier.GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory.getVolFrac();
  const ::amrex::MultiFab& datas = hier.GetPatchLevel(level).data;
  ForEachFab(alphas, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    const ::amrex::FArrayBox& data = datas[mfi];
    IndexBox<AMREX_SPACEDIM> section =
        AsIndexBox<AMREX_SPACEDIM>(box & mfi.tilebox());
    for (int comp = 0; comp < ncomp; ++comp) {
      ForEachIndex(section, [&](auto... is) {
        ::amrex::IntVect index{static_cast<int>(is)...};
        const double frac = alpha(index);
        state_buffer[comp] += frac / total_volume * data(index, comp);
      });
    }
  });
  std::vector<double> global_state_buffer(static_cast<std::size_t>(ncomp));
  MPI_Allreduce(state_buffer.data(), global_state_buffer.data(), ncomp,
                MPI_DOUBLE, MPI_SUM,
                ::amrex::ParallelDescriptor::Communicator());
  int comp = 0;
  ForEachComponent(
      [&comp, &global_state_buffer](auto&& var) {
        var = global_state_buffer[comp];
        comp += 1;
      },
      state);
}

template <typename GriddingAlgorithm>
int FindLevel(const ::amrex::Geometry& geom,
              const GriddingAlgorithm& gridding) {
  for (int level = 0; level < gridding.GetPatchHierarchy().GetNumberOfLevels();
       ++level) {
    if (geom.Domain() ==
        gridding.GetPatchHierarchy().GetGeometry(level).Domain()) {
      return level;
    }
  }
  return -1;
}

} // namespace

void IsentropicPressureBoundary::FillBoundary(::amrex::MultiFab& mf,
                                              const ::amrex::Geometry& geom,
                                              Duration t,
                                              const GriddingAlgorithm& grid) {
  Complete<IdealGasMix<AMREX_SPACEDIM>> state(equation_);
  int level = FindLevel(geom, grid);
  ::amrex::Box refined_inner_box = options_.coarse_inner_box;
  for (int l = 1; l <= level; ++l) {
    refined_inner_box.refine(
        grid.GetPatchHierarchy().GetRatioToCoarserLevel(l));
  }
  AverageState(state, grid.GetPatchHierarchy(), level, refined_inner_box);
  equation_.CompleteFromCons(state, state);

  boost::log::sources::severity_channel_logger<
      boost::log::trivial::severity_level>
      log(boost::log::keywords::channel = options_.channel_name,
          boost::log::keywords::severity = boost::log::trivial::debug);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", t.count());
  double rho = state.density;
  double u = state.momentum[0] / rho;
  double p = state.pressure;
  BOOST_LOG(log) << fmt::format("Average inner state: {} kg/m3, {} m/s, {} Pa",
                                rho, u, p);

  equation_.GetReactor().SetDensity(state.density);
  equation_.GetReactor().SetMassFractions(state.species);
  equation_.GetReactor().SetTemperature(state.temperature);
  equation_.GetReactor().SetPressure(options_.outer_pressure);
  equation_.CompleteFromReactor(state);
  IsentropicExpansionWithoutDissipation(equation_, state, state, p);
  if (options_.side == 1) {
    state.momentum[0] = -state.momentum[0];
  }
  rho = state.density;
  u = state.momentum[0] / rho;
  p = state.pressure;
  BOOST_LOG(log) << fmt::format("Outer State: {} kg/m3, {} m/s, {} Pa", rho, u,
                                p);
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  FillBoundary(mf, alphas, geom, state);
}

void IsentropicPressureBoundary::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
    const ::amrex::Geometry& geom,
    const Complete<IdealGasMix<AMREX_SPACEDIM>>& state) {
  const int ngrow = mf.nGrow(int(options_.direction));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }
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
        auto states = MakeView<Complete<IdealGasMix<AMREX_SPACEDIM>>>(
            fab, equation_, mfi.growntilebox());
        ForEachIndex(box_to_fill, [&alpha, &state, &states](auto... is) {
          std::array<std::ptrdiff_t, AMREX_SPACEDIM> dest{int(is)...};
          ::amrex::IntVect iv{int(is)...};
          if (alpha(iv) > 0.0) {
            Store(states, state, dest);
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
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom, Duration /* dt */,
    const GriddingAlgorithm& grid) {
  int level = FindLevel(geom, grid);
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
        ForEachIndex(box_to_fill, [this, x0, dir_v, &alpha, &state, &states](auto... is) {
          std::array<std::ptrdiff_t, AMREX_SPACEDIM> dest{int(is)...};
          std::array<std::ptrdiff_t, AMREX_SPACEDIM> src = dest;
          src[dir_v] = x0;
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
