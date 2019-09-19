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

namespace fub::amrex::cutcell {
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
  const double h_before = dest.energy / dest.density + dest.pressure / dest.density;
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
    const std::string& name, const IdealGasMix<AMREX_SPACEDIM>& eq,
    const ::amrex::Box& coarse_inner_box, double outer_pressure, Direction dir,
    int side)
    : log_(boost::log::keywords::channel = name),
      time_attr_{0.0}, equation_{eq}, coarse_inner_box_{coarse_inner_box},
      outer_pressure_{outer_pressure}, dir_{dir}, side_{side} {
  log_.add_attribute("Time", time_attr_);
}

namespace {
// double CellVolume(const ::amrex::Geometry& geom) {
//  return AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
//}

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
  AverageState(state, grid.GetPatchHierarchy(), 0, coarse_inner_box_);
  equation_.CompleteFromCons(state, state);
  time_attr_.set(t.count());
  double rho = state.density;
  double u = state.momentum[0] / rho;
  double p = state.pressure;
  BOOST_LOG(log_) << fmt::format("Average inner state: {} kg/m3, {} m/s, {} Pa",
                                 rho, u, p);

  equation_.GetReactor().SetMassFractions(state.species);
  equation_.GetReactor().SetTemperature(300.0);
  equation_.GetReactor().SetPressure(outer_pressure_);
  equation_.CompleteFromReactor(state);
  IsentropicExpansionWithoutDissipation(equation_, state, state, p);
  rho = state.density;
  u = state.momentum[0] / rho;
  p = state.pressure;
  BOOST_LOG(log_) << fmt::format("Outer State: {} kg/m3, {} m/s, {} Pa", rho, u,
                                 p);
  int level = FindLevel(geom, grid);
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  FillBoundary(mf, alphas, geom, state);
}

void IsentropicPressureBoundary::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
    const ::amrex::Geometry& geom,
    const Complete<IdealGasMix<AMREX_SPACEDIM>>& state) {
  const int ngrow = mf.nGrow(int(dir_));
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
