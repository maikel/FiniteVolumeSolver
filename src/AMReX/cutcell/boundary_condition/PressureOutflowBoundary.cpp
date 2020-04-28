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

#include "fub/AMReX/cutcell/boundary_condition/PressureOutflowBoundary.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

namespace fub::amrex::cutcell {

PressureOutflowOptions::PressureOutflowOptions(const ProgramOptions& options) {
  outer_pressure = GetOptionOr(options, "outer_pressure", outer_pressure);
  direction = GetOptionOr(options, "direction", direction);
  side = GetOptionOr(options, "side", side);
}

void PressureOutflowOptions::Print(SeverityLogger& log) {
  BOOST_LOG(log) << " - outer_pressure = " << outer_pressure;
  BOOST_LOG(log) << " - direction = " << int(direction);
  BOOST_LOG(log) << " - side = " << side;
}

PressureOutflowBoundary::PressureOutflowBoundary(
    const PerfectGas<AMREX_SPACEDIM>& eq, const PressureOutflowOptions& options)
    : equation_(eq), options_(options) {}

void PressureOutflowBoundary::FillBoundary(::amrex::MultiFab& mf,
                                           const GriddingAlgorithm& grid,
                                           int level) {
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
  Complete<PerfectGas<AMREX_SPACEDIM>> zeros{equation_};
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
          std::array<std::ptrdiff_t, AMREX_SPACEDIM> src = dest;
          const auto d = static_cast<std::size_t>(dir_v);
          src[d] = x0;
          ::amrex::IntVect src_iv{
              AMREX_D_DECL(int(src[0]), int(src[1]), int(src[2]))};
          if (alpha(iv) > 0.0 && alpha(src_iv) > 0.0) {
            Load(state, states, src);
            const double c = state.speed_of_sound;
            FUB_ASSERT(c > 0.0);
            const double c2 = c * c;
            auto velocity = equation_.Velocity(state);
            const double velocity_norm2 = velocity.matrix().squaredNorm();
            if (velocity_norm2 <= c2) {
              const double rho_new = kappa * pb / c2;
              state = equation_.CompleteFromPrim(rho_new, velocity, pb);
            }
            Store(states, state, dest);
          } else {
            Store(states, zeros, dest);
          }
        });
      }
    }
  });
}

} // namespace fub::amrex::cutcell
