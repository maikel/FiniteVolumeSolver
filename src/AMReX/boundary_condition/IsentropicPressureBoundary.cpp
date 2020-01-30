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

#include "fub/AMReX/boundary_condition/IsentropicPressureBoundary.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/ForEach.hpp"

#include <cmath>

namespace fub::amrex {

namespace {
std::array<std::ptrdiff_t, 1>
MapToSrc(const std::array<std::ptrdiff_t, 1>& dest,
         const ::amrex::Geometry& geom, int side, Direction dir) {
  const int boundary = (side == 0) * geom.Domain().smallEnd(int(dir)) +
                       (side == 1) * geom.Domain().bigEnd(int(dir));
  const int distance = dest[int(dir)] - boundary;
  const int sign = int((distance > 0) - (distance < 0));
  std::array<std::ptrdiff_t, 1> src{dest};
  src[int(dir)] -= 2 * distance - sign;
  //  src[std::size_t(dir)] = boundary;
  return src;
}

int Sign(double x) { return (x > 0) - (x < 0); }

void IsentropicExpansionWithoutDissipation(IdealGasMix<1>& eq,
                                           Complete<IdealGasMix<1>>& dest,
                                           const Complete<IdealGasMix<1>>& src,
                                           double dest_pressure,
                                           double efficiency = 1.0) {
  double old_velocity = src.momentum[0] / src.density;
  eq.SetReactorStateFromComplete(src);
  eq.CompleteFromReactor(dest);
  const double h_before =
      dest.energy / dest.density + dest.pressure / dest.density;
  eq.GetReactor().SetPressureIsentropic(dest_pressure);
  eq.CompleteFromReactor(dest);
  const double h_after =
      dest.energy / dest.density + dest.pressure / dest.density;
  const double enthalpyDifference = h_before - h_after;
  const double u_border = [&] {
    if (Sign(old_velocity) == Sign(enthalpyDifference)) {
      return Sign(enthalpyDifference) *
             std::sqrt(efficiency * std::abs(enthalpyDifference) * 2 +
                       old_velocity * old_velocity);
    } else {
      double inner = efficiency * std::abs(enthalpyDifference) * 2 -
                     old_velocity * old_velocity;
      return inner > 0.0 ? Sign(enthalpyDifference) * std::sqrt(std::abs(inner))
                         : Sign(inner) * std::sqrt(std::abs(inner));
    }
  }();
  dest.momentum[0] = dest.density * u_border;
  dest.energy += 0.5 * u_border * dest.momentum[0];
}

} // namespace

IsentropicPressureBoundary::IsentropicPressureBoundary(const IdealGasMix<1>& eq,
                                                       double outer_pressure,
                                                       Direction dir, int side)
    : equation_{eq}, outer_pressure_{outer_pressure}, dir_{dir}, side_{side} {}

void IsentropicPressureBoundary::FillBoundary(::amrex::MultiFab& mf,
                                              const ::amrex::Geometry& geom,
                                              Duration,
                                              const GriddingAlgorithm&) {
  FillBoundary(mf, geom);
}

void IsentropicPressureBoundary::FillBoundary(::amrex::MultiFab& mf,
                                              const ::amrex::Geometry& geom,
                                              Duration,
                                              const GriddingAlgorithm&,
                                              Direction dir) {
  if (dir == dir_) {
    FillBoundary(mf, geom);
  }
}

void IsentropicPressureBoundary::FillBoundary(::amrex::MultiFab& mf,
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
          IsentropicExpansionWithoutDissipation(equation_, state, state,
                                                outer_pressure_);
          Store(states, state, dest);
        });
      }
    }
  });
}
} // namespace fub::amrex
