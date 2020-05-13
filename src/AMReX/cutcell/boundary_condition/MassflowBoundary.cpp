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

#include "fub/AMReX/cutcell/boundary_condition/MassflowBoundary.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AnyBoundaryCondition.hpp"
#include "fub/ForEach.hpp"

namespace fub::amrex::cutcell {
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

MassflowBoundary::MassflowBoundary(const IdealGasMix<AMREX_SPACEDIM>& eq,
                                   const MassflowBoundaryOptions& options)
    : equation_{eq}, options_{options} {}

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
        const auto c = static_cast<std::size_t>(comp);
        state_buffer[c] += frac / total_volume * data(index, comp);
      });
    }
  });
  std::vector<double> global_state_buffer(static_cast<std::size_t>(ncomp));
  MPI_Allreduce(state_buffer.data(), global_state_buffer.data(), ncomp,
                MPI_DOUBLE, MPI_SUM,
                ::amrex::ParallelDescriptor::Communicator());
  CopyFromBuffer(state, global_state_buffer);
}

} // namespace

void MassflowBoundary::FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid, int level) {
  boost::log::sources::severity_channel_logger<
      boost::log::trivial::severity_level>
      log(boost::log::keywords::channel = options_.channel_name,
          boost::log::keywords::severity = boost::log::trivial::debug);
  const Duration t = grid.GetTimePoint();
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", t.count());

  Complete<IdealGasMix<AMREX_SPACEDIM>> state(equation_);
  AverageState(state, grid.GetPatchHierarchy(), 0, options_.coarse_inner_box);
  equation_.CompleteFromCons(state, state);
  double rho = state.density;
  double u = state.momentum[int(options_.dir)] / rho;
  double p = state.pressure;
  BOOST_LOG(log) << fmt::format("Average inner state: {} kg/m3, {} m/s, {} Pa",
                                rho, u, p);

  equation_.GetReactor().SetDensity(state.density);
  equation_.GetReactor().SetMassFractions(state.species);
  equation_.GetReactor().SetTemperature(state.temperature);
  //  const double rho = state.density
  const double gamma = state.gamma;
  const double c = state.speed_of_sound;
  //  const double u = state.momentum[int(options_.dir)] / rho;
  const double Ma = u / c;
  const double gammaMinus = gamma - 1.0;
  const double gammaPlus = gamma + 1.0;
  const double rGammaPlus = 1.0 / gammaPlus;
  const double c_critical =
      std::sqrt(c * c + 0.5 * gammaMinus * u * u) * std::sqrt(2 * rGammaPlus);
  const double u_n = options_.required_massflow / rho / options_.surface_area;
  const double lambda = u / c_critical;
  const double lambda_n = u_n / c_critical;
  const double gammaQuot = gammaMinus * rGammaPlus;
  //  const double p = state.pressure;
  const double p0_n =
      p * std::pow(1. - gammaQuot * lambda_n * lambda_n, -gamma / gammaMinus);
  const double p_n =
      p0_n * std::pow(1. - gammaQuot * lambda * lambda, gamma / gammaMinus);

  equation_.GetReactor().SetPressureIsentropic(p_n);
  Eigen::Array<double, AMREX_SPACEDIM, 1> velocity =
      Eigen::Array<double, AMREX_SPACEDIM, 1>::Zero();
  velocity[int(options_.dir)] = u_n;

  equation_.CompleteFromReactor(state, velocity);

  rho = state.density;
  u = u_n;
  p = state.pressure;
  BOOST_LOG(log) << fmt::format("Outer State: {} kg/m3, {} m/s, {} Pa, Ma = {}",
                                rho, u, p, Ma);
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  FillBoundary(mf, alphas, geom, state);
}

void MassflowBoundary::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
    const ::amrex::Geometry& geom,
    const Complete<IdealGasMix<AMREX_SPACEDIM>>& state) {
  const int ngrow = mf.nGrow(int(options_.dir));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  //  Complete<IdealGasMix<AMREX_SPACEDIM>> state{equation_};
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

} // namespace fub::amrex::cutcell
