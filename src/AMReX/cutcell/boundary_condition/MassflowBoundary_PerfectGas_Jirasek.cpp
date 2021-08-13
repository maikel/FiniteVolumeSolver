// Copyright (c) 2021 Christian Zenker
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

#include "fub/AMReX/cutcell/boundary_condition/MassflowBoundary_PerfectGas_Jirasek.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/cutcell/AverageState.hpp"
#include "fub/AnyBoundaryCondition.hpp"
#include "fub/ForEach.hpp"

namespace fub::amrex::cutcell {
MassflowBoundary_PerfectGas_JirasekOptions::
    MassflowBoundary_PerfectGas_JirasekOptions(const ProgramOptions& options) {
  FUB_GET_OPTION_VAR(options, channel_name);
  FUB_GET_OPTION_VAR(options, coarse_inner_box);
  FUB_GET_OPTION_VAR(options, required_massflow);
  FUB_GET_OPTION_VAR(options, ramp_time);
  FUB_GET_OPTION_VAR(options, surface_area);
  FUB_GET_OPTION_VAR(options, dir);
  FUB_GET_OPTION_VAR(options, side);
}

void MassflowBoundary_PerfectGas_JirasekOptions::Print(SeverityLogger& log) const {
  BOOST_LOG(log) << fmt::format(" - required_massflow = {} [kg/s]",
                                required_massflow);
  BOOST_LOG(log) << fmt::format(" - surface_area = {} [m2]", surface_area);
  BOOST_LOG(log) << fmt::format(" - direction = {} [-]", static_cast<int>(dir));
  BOOST_LOG(log) << fmt::format(" - side = {} [-]", side);
  BOOST_LOG(log) << fmt::format(" - ramp_time = {} [s]", ramp_time.count());
  std::array<int, AMREX_SPACEDIM> lower, upper;
  std::copy_n(coarse_inner_box.smallEnd().getVect(), AMREX_SPACEDIM,
              lower.data());
  std::copy_n(coarse_inner_box.bigEnd().getVect(), AMREX_SPACEDIM,
              upper.data());
  BOOST_LOG(log) << fmt::format(" - coarse_inner_box = {{{{{}}}, {{{}}}}} [-]",
                                lower, upper);
}

MassflowBoundary_PerfectGas_Jirasek::MassflowBoundary_PerfectGas_Jirasek(
    const PerfectGas<AMREX_SPACEDIM>& eq,
    const MassflowBoundary_PerfectGas_JirasekOptions& options)
    : equation_{eq}, options_{options} {}

void MassflowBoundary_PerfectGas_Jirasek::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& grid, int level) {
  boost::log::sources::severity_channel_logger<
      boost::log::trivial::severity_level>
      log(boost::log::keywords::channel = options_.channel_name,
          boost::log::keywords::severity = boost::log::trivial::debug);
  const Duration current_time = grid.GetTimePoint();
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", current_time.count());

  Complete<PerfectGas<AMREX_SPACEDIM>> state(equation_);
  AverageState(state, grid.GetPatchHierarchy(), 0, options_.coarse_inner_box);
  equation_.CompleteFromCons(state, state);
  double rho = state.density;
  double u = state.momentum[int(options_.dir)] / rho;
  double p = state.pressure;

  const double gamma = equation_.gamma;
  const double c = state.speed_of_sound;
  const double Ma = u / c;
  BOOST_LOG(log) << fmt::format(
      "Average inner state: {} kg/m3, {} m/s, {} Pa, Ma = {}", rho, u, p, Ma);

  const double gammaMinus = gamma - 1.0;
  const double gammaPlus = gamma + 1.0;
  const double rGammaPlus = 1.0 / gammaPlus;
  const double c_critical =
      std::sqrt(c * c + 0.5 * gammaMinus * u * u) * std::sqrt(2 * rGammaPlus);
  
  static const double eps = std::numeric_limits<double>::epsilon();
  const double ramp_time = options_.ramp_time.count();
  const double ramp_factor = current_time.count() / (ramp_time+eps);
  const double ramp_coeff = (current_time.count() < ramp_time) ? ramp_factor : 1.0;
  BOOST_LOG(log) << fmt::format(
      "MassFlow: time {}, offset {}, ramp_coeff {}", current_time.count(), ramp_time, ramp_coeff);
  const double u_n = ramp_coeff * options_.required_massflow / rho / options_.surface_area;
  const double lambda = u / c_critical;
  const double lambda_n = u_n / c_critical;
  const double gammaQuot = gammaMinus * rGammaPlus;
  const double p0_n =
      p * std::pow(1. - gammaQuot * lambda_n * lambda_n, -gamma / gammaMinus);
  const double p_n =
      p0_n * std::pow(1. - gammaQuot * lambda * lambda, gamma / gammaMinus);

  Eigen::Array<double, AMREX_SPACEDIM, 1> velocity =
      Eigen::Array<double, AMREX_SPACEDIM, 1>::Zero();
  velocity[int(options_.dir)] = u_n;

  Complete<PerfectGas<AMREX_SPACEDIM>> new_state =
      equation_.CompleteFromPrim(state.density, velocity, p_n);

  rho = new_state.density;
  const double Ma_n = u_n / new_state.speed_of_sound;
  p = new_state.pressure;
  BOOST_LOG(log) << fmt::format("Outer State: {} kg/m3, {} m/s, {} Pa, Ma = {}",
                                rho, u_n, p, Ma_n);
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  FillBoundary(mf, alphas, geom, new_state);
}

void MassflowBoundary_PerfectGas_Jirasek::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
    const ::amrex::Geometry& geom,
    const Complete<PerfectGas<AMREX_SPACEDIM>>& state) {
  const int ngrow = mf.nGrow(int(options_.dir));
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
      ::amrex::Box shifted = ::amrex::shift(boundary, int(options_.dir),
                                            GetSign(options_.side) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<PerfectGas<AMREX_SPACEDIM>>>(
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

} // namespace fub::amrex::cutcell
