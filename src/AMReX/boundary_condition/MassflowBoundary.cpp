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
#include "fub/ForEach.hpp"

namespace fub::amrex {
namespace {
inline int GetSign(int side) { return (side == 0) - (side == 1); }
} // namespace
MassflowBoundary::MassflowBoundary(const IdealGasMix<AMREX_SPACEDIM>& eq,
                                   const ::amrex::Box& coarse_inner_box,
                                   double required_massflow,
                                   double surface_area, Direction dir, int side)
    : MassflowBoundary("MassflowBoundary", eq, coarse_inner_box,
                       required_massflow, surface_area, dir, side) {}

MassflowBoundary::MassflowBoundary(const std::string& name,
                                   const IdealGasMix<AMREX_SPACEDIM>& eq,
                                   const ::amrex::Box& coarse_inner_box,
                                   double required_massflow,
                                   double surface_area, Direction dir, int side)
    : log_(boost::log::keywords::channel = name), equation_{eq},
      coarse_inner_box_{coarse_inner_box},
      required_massflow_{required_massflow},
      surface_area_{surface_area}, dir_{dir}, side_{side} {}

namespace {
double TotalVolume(const PatchHierarchy& hier, int level,
                   const ::amrex::Box& box) {
  double local_volume(0.0);
  const ::amrex::Geometry& geom = hier.GetGeometry(level);
  const double cell_volume = goem.CellVolume();
  return box.numcells() * cell_volume;
}

template <int Rank>
void AverageState(Complete<IdealGasMix<Rank>>& state,
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

} // namespace

void MassflowBoundary::FillBoundary(::amrex::MultiFab& mf,
                                    const GriddingAlgorithm& grid, int level,
                                    Direction dir) {
  if (dir == dir_) {
    FillBoundary(mf, grid, level);
  }
}

void MassflowBoundary::FillBoundary(::amrex::MultiFab& mf,
                                    const GriddingAlgorithm& grid, int level) {
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  Complete<IdealGasMix<AMREX_SPACEDIM>> state(equation_);
  AverageState(state, grid.GetPatchHierarchy(), level, coarse_inner_box_);
  equation_.CompleteFromCons(state, state);
  double rho = state.density;
  double u = state.momentum[int(dir_)] / rho;
  double p = state.pressure;
  BOOST_LOG(log_) << fmt::format(
      "Average inner state: rho = {} [kg/m3], u = {} [m/s], p = {} [Pa]", rho,
      u, p);

  equation_.GetReactor().SetDensity(state.density);
  equation_.GetReactor().SetMassFractions(state.species);
  equation_.GetReactor().SetTemperature(state.temperature);
  const double gamma = state.gamma;
  const double c = state.speed_of_sound;
  const double Ma = u / c;
  const double gammaMinus = gamma - 1.0;
  const double gammaPlus = gamma + 1.0;
  const double rGammaPlus = 1.0 / gammaPlus;
  const double c_critical =
      std::sqrt(c * c + 0.5 * gammaMinus * u * u) * std::sqrt(2 * rGammaPlus);
  const double u_n = required_massflow_ / rho / surface_area_;
  const double lambda = u / c_critical;
  const double lambda_n = u_n / c_critical;
  const double gammaQuot = gammaMinus * rGammaPlus;
  const double p0_n =
      p * std::pow(1. - gammaQuot * lambda_n * lambda_n, -gamma / gammaMinus);
  const double p_n =
      p0_n * std::pow(1. - gammaQuot * lambda * lambda, gamma / gammaMinus);

  equation_.GetReactor().SetPressureIsentropic(p_n);
  Eigen::Array<double, AMREX_SPACEDIM, 1> velocity =
      Eigen::Array<double, AMREX_SPACEDIM, 1>::Zero();
  velocity[int(dir_)] = u_n;

  equation_.CompleteFromReactor(state, velocity);

  rho = state.density;
  u = u_n;
  p = state.pressure;
  BOOST_LOG(log_) << boost::log::add_value("Time", t.count())
                  << boost::log::add_value("Level", level)
                  << fmt::format("outer state: rho = {} [kg/m3], u = {} [m/s], "
                                 "p = {} [Pa], Ma = {} [-]",
                                 rho, u, p, Ma);
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  FillBoundary(mf, alphas, geom, state);
}

void MassflowBoundary::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
    const Complete<IdealGasMix<AMREX_SPACEDIM>>& state) {
  const int ngrow = mf.nGrow(int(dir_));
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
      ::amrex::Box shifted =
          ::amrex::shift(boundary, int(dir_), GetSign(side_) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<IdealGasMix<AMREX_SPACEDIM>>>(
            fab, equation_, mfi.growntilebox());
        ForEachIndex(box_to_fill, [&alpha, &state, &states](auto... is) {
          std::array<std::ptrdiff_t, AMREX_SPACEDIM> dest{int(is)...};
          Store(states, state, dest);
        });
      }
    }
  });
}

} // namespace fub::amrex
