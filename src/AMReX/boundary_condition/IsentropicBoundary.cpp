#include "fub/AMReX/boundary_condition/IsentropicBoundary.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/ForEach.hpp"

namespace fub::amrex {

namespace {
std::array<std::ptrdiff_t, 1>
MapToSrc(const std::array<std::ptrdiff_t, 1>& dest,
         const ::amrex::Geometry& geom, int side, Direction dir) {
  const int boundary = (side == 0) * geom.Domain().smallEnd(int(dir)) +
                       (side == 1) * geom.Domain().bigEnd(int(dir));
  //    const int distance = dest[int(dir)] - boundary;
  //  const int sign = int((distance > 0) - (distance < 0));
  std::array<std::ptrdiff_t, 1> src{dest};
  //  src[int(dir)] -= 2 * distance - sign;
  src[std::size_t(dir)] = boundary;
  return src;
}

} // namespace

IsentropicBoundary::IsentropicBoundary(const IdealGasMix<1>& eq,
                                       double outer_pressure, Direction dir,
                                       int side)
    : equation_{eq}, outer_pressure_{outer_pressure}, dir_{dir}, side_{side} {}

void IsentropicBoundary::FillBoundary(::amrex::MultiFab& mf,
                                      const ::amrex::Geometry& geom, Duration,
                                      const GriddingAlgorithm&) {
  FillBoundary(mf, geom);
}

void IsentropicBoundary::FillBoundary(::amrex::MultiFab& mf,
                                      const ::amrex::Geometry& geom) {
  const int ngrow = mf.nGrow(int(dir_));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  Complete<IdealGasMix<1>> state{equation_};
  if (boundaries.isEmpty()) {
    return;
  }
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    for (const ::amrex::Box& boundary : boundaries) {
      ::amrex::Box shifted =
          ::amrex::shift(boundary, int(dir_), GetSign(side_) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<IdealGasMix<1>>>(fab, equation_,
                                                         mfi.growntilebox());
        ForEachIndex(box_to_fill, [this, &geom, &state,
                                   &states](std::ptrdiff_t i, auto...) {
          std::array<std::ptrdiff_t, 1> dest{i};
          std::array<std::ptrdiff_t, 1> src = MapToSrc(dest, geom, side_, dir_);
          Load(state, states, src);
          equation_.SetReactorStateFromComplete(state);
          const double rhoE_prev = state.energy;
          equation_.GetReactor().SetPressureIsentropic(outer_pressure_);
          equation_.CompleteFromReactor(state, state.momentum / state.density);
          const double rhoE_after = state.energy;
          const double rhoE_difference = rhoE_after - rhoE_prev;
          const double velocity =
              std::sqrt(std::abs(2.0 * rhoE_difference / state.density));
          state.momentum += state.density * velocity;
          state.energy += rhoE_difference;
          Store(states, state, dest);
        });
      }
    }
  });
}
} // namespace fub::amrex
