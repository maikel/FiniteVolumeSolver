// Copyright (c) 2020 Maikel Nadolski
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

#ifndef FUB_AMREX_CUTCELL_PRESSURE_OUTFLOW_BOUNDARY_HPP
#define FUB_AMREX_CUTCELL_PRESSURE_OUTFLOW_BOUNDARY_HPP

#include "fub/AMReX/cutcell/boundary_condition/PressureOutflowBoundary.hpp"

namespace fub::amrex::cutcell {
namespace {
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

PressureOutflowOptions::PressureOutflowOptions(
    const PerfectGas<AMREX_SPACEDIM>& eq, const ProgramOptions& options)
    : equation_(eq), options_(options) {}

void PressureOutlfowBoundary::FillBoundary(::amrex::MultiFab& mf,
                                           const ::amrex::Geometry& geom,
                                           Duration dt,
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
  const double pb = options_.outer_pressure;
  const double kappa = equation_.gamma;
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
        ForEachIndex(box_to_fill, [&](auto... is) {
          std::array<std::ptrdiff_t, AMREX_SPACEDIM> dest{int(is)...};
          ::amrex::IntVect iv{int(is)...};
          if (alpha(iv) > 0.0) {
            std::array<std::ptrdiff_t, AMREX_SPACEDIM> src = dest;
            src[dir_v] = x0;
            Load(state, states, src);
            const double c = state.speed_of_sound;
            const double c2 = c * c;
            auto velocity = equation_.Velocity(state);
            const double double velocity_norm2 = velocity.squaredNorm();
            if (velocity_norm2 <= c2) {
              Store(states, state, data);
            } else {
              const double rho_new = kappa * pb / c2;
              state = equation_.CompleteFromPrim(rho_new, velocity, pb);
            }
            Store(states, state, dest);
          }
        });
      }
    }
  });
}

} // namespace fub::amrex::cutcell

#endif // !FUB_AMREX_CUTCELL_PRESSURE_OUTFLOW_BOUNDARY_HPP
