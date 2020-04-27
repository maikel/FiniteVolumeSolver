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

#include "fub/AMReX/boundary_condition/MassflowBoundary.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AnyBoundaryCondition.hpp"
#include "fub/ForEach.hpp"

#include "fub/AMReX/AverageState.hpp"

#include "fub/ext/Log.hpp"

namespace fub::amrex {

MassflowBoundaryOptions::MassflowBoundaryOptions(
    const ProgramOptions& options) {
  channel_name = GetOptionOr(options, "channel_name", channel_name);
  coarse_inner_box = GetOptionOr(options, "coarse_inner_box", coarse_inner_box);
  required_massflow =
      GetOptionOr(options, "required_massflow", required_massflow);
  surface_area = GetOptionOr(options, "surface_area", surface_area);
  dir = GetOptionOr(options, "direction", dir);
  side = GetOptionOr(options, "side", side);
}

void MassflowBoundaryOptions::Print(SeverityLogger& log) const {
  BOOST_LOG(log) << fmt::format(" - required_massflow = {} [kg/s]",
                                required_massflow);
  BOOST_LOG(log) << fmt::format(" - surface_area = {} [m2]", surface_area);
  BOOST_LOG(log) << fmt::format(" - direction = {} [-]", static_cast<int>(dir));
  BOOST_LOG(log) << fmt::format(" - side = {} [-]", side);
  std::array<int, AMREX_SPACEDIM> lower, upper;
  std::copy_n(coarse_inner_box.smallEnd().getVect(), AMREX_SPACEDIM,
              lower.data());
  std::copy_n(coarse_inner_box.bigEnd().getVect(), AMREX_SPACEDIM,
              upper.data());
  BOOST_LOG(log) << fmt::format(" - coarse_inner_box = {{{{{}}}, {{{}}}}} [-]",
                                lower, upper);
}

template <int Rank>
MassflowBoundary<Rank>::MassflowBoundary(const IdealGasMix<Rank>& eq,
                                         const MassflowBoundaryOptions& options)
    : equation_{eq}, options_{options} {}

template <int Rank>
void MassflowBoundary<Rank>::FillBoundary(::amrex::MultiFab& mf,
                                          const GriddingAlgorithm& grid,
                                          int level, Direction dir) {
  if (dir == options_.dir) {
    FillBoundary(mf, grid, level);
  }
}

template <int Rank>
void MassflowBoundary<Rank>::ComputeBoundaryState(
    Complete<IdealGasMix<Rank>>& boundary,
    const Complete<IdealGasMix<Rank>>& state) {
  const double rho = state.density;
  const double u = state.momentum[int(options_.dir)] / rho;
  const double p = state.pressure;
  const double gamma = state.gamma;
  const double c = state.speed_of_sound;
  const double gammaMinus = gamma - 1.0;
  const double gammaPlus = gamma + 1.0;
  const double rGammaPlus = 1.0 / gammaPlus;
  const double c_critical =
      std::sqrt(c * c + 0.5 * gammaMinus * u * u) * std::sqrt(2 * rGammaPlus);
  const double u_n = options_.required_massflow / rho / options_.surface_area;
  const double lambda = u / c_critical;
  const double lambda_n = u_n / c_critical;
  const double gammaQuot = gammaMinus * rGammaPlus;
  const double p0_n =
      p * std::pow(1. - gammaQuot * lambda_n * lambda_n, -gamma / gammaMinus);
  const double p_n =
      p0_n * std::pow(1. - gammaQuot * lambda * lambda, gamma / gammaMinus);

  equation_.GetReactor().SetDensity(state.density);
  equation_.GetReactor().SetMassFractions(state.species);
  equation_.GetReactor().SetTemperature(state.temperature);
  equation_.GetReactor().SetPressureIsentropic(p_n);
  Eigen::Array<double, Rank, 1> velocity =
      Eigen::Array<double, Rank, 1>::Zero();
  velocity[int(options_.dir)] = u_n;
  equation_.CompleteFromReactor(boundary, velocity);
}

template <int Rank>
void MassflowBoundary<Rank>::FillBoundary(::amrex::MultiFab& mf,
                                          const GriddingAlgorithm& grid,
                                          int level) {
  boost::log::sources::severity_channel_logger<
      boost::log::trivial::severity_level>
      log(boost::log::keywords::channel = options_.channel_name,
          boost::log::keywords::severity = boost::log::trivial::debug);
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  const Duration t = grid.GetPatchHierarchy().GetTimePoint(level);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", t.count());

  Complete<IdealGasMix<Rank>> state(equation_);
  AverageState(state, mf, geom, options_.coarse_inner_box);
  equation_.CompleteFromCons(state, state);
  double rho = state.density;
  double u = state.momentum[static_cast<int>(options_.dir)] / rho;
  double p = state.pressure;
  BOOST_LOG(log) << fmt::format(
      "Average inner state: rho = {} [kg/m3], u = {} [m/s], p = {} [Pa]", rho,
      u, p);

  ComputeBoundaryState(state, state);

  rho = state.density;
  u = state.momentum[static_cast<int>(options_.dir)] / rho;
  p = state.pressure;
  const double c = state.speed_of_sound;
  const double Ma = u / c;
  BOOST_LOG(log) << fmt::format("Outer state: rho = {} [kg/m3], u = {} [m/s], "
                                "p = {} [Pa], Ma = {} [-]",
                                rho, u, p, Ma);
  FillBoundary(mf, geom, state);
}

template <int Rank>
void MassflowBoundary<Rank>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
    const Complete<IdealGasMix<Rank>>& state) {
  const int ngrow = mf.nGrow(int(options_.dir));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    for (const ::amrex::Box& boundary : boundaries) {
      ::amrex::Box shifted = ::amrex::shift(boundary, int(options_.dir),
                                            GetSign(options_.side) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<IdealGasMix<Rank>>>(fab, equation_,
                                                            mfi.growntilebox());
        ForEachIndex(AsIndexBox<Rank>(box_to_fill),
                     [&state, &states](auto... is) {
                       Index<Rank> dest{is...};
                       Store(states, state, dest);
                     });
      }
    }
  });
}

template class MassflowBoundary<1>;
template class MassflowBoundary<2>;
template class MassflowBoundary<3>;

} // namespace fub::amrex
