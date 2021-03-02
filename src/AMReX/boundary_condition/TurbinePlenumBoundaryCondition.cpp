// Copyright (c) 2021 Maikel Nadolski
// Copyright (c) 2021 Christian Zenker
// Copyright (c) 2021 Rupert Klein
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

#include "fub/AMReX/boundary_condition/TurbinePlenumBoundaryCondition.hpp"
#include "fub/ForEach.hpp"
#include "fub/AMReX/ForEachFab.hpp"

namespace fub::amrex {

namespace {
template <typename I>
static I MapToSrc(I& dest, const ::amrex::Geometry& geom, int side,
                  Direction dir) noexcept {
  const int boundary = (side == 0) * geom.Domain().smallEnd(int(dir)) +
                       (side == 1) * geom.Domain().bigEnd(int(dir));
  I src{dest};
  src[int(dir)] = boundary;
  return src;
}
} // namespace

void TurbinePlenumBoundaryCondition::FillBoundary(::amrex::MultiFab& mf,
                                                  const GriddingAlgorithm& grid,
                                                  int level)
{
  /* open b.c's with pressure relaxation to a given value */
  // Fetch Turbine Plenum state
  const double pflush = control_state_->turbine.pressure;
  const double Tflush = control_state_->turbine.temperature;
  const double rhoflush = pflush / Tflush;// * equation_.ooRspec;

  const int ngrow = mf.nGrow(0);
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }
  double Msqinv = 1.0/equation_.Msq;
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    Complete<PerfectGasMix<1>> state(equation_);
    for (const ::amrex::Box& boundary : boundaries) {
      ::amrex::Box shifted = ::amrex::shift(boundary, 0, GetSign(1) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<PerfectGasMix<1>>>(fab, equation_,
                                                           mfi.growntilebox());
        ForEachIndex(
            AsIndexBox<1>(box_to_fill), [&](std::ptrdiff_t i) {
              Index<1> dest{i};
              Index<1> image = MapToSrc(dest, geom, 1, Direction::X);
              Load(state, states, image);
              double rho = state.density;
              double u = state.momentum[0] / rho;
              double ek = 0.5 * rho * (u * u);
              double p = state.pressure;
              double T = p / rho;// * equation_.ooRspec;
              double pratio = pflush / p;
              double rhoratio = pow(pratio, equation_.gamma_inv);
              double Tratio = pratio / rhoratio;
              double uout{};
              double rhoout{};
              if (pflush > p) {
                if (u <= 0.0) {
                  /* flow from plenum into tube */
                  // std::sqrt(2.0 * eq.gamma_over_gamma_minus_one *
                  // std::max(0.0, Tpv - Tin));
                  uout = -std::sqrt(
                      std::max(0.0, 2.0 * equation_.gamma_over_gamma_minus_one *
                                        (Tratio - 1.0) * T)) * Msqinv;
                  rhoout = rhoflush / rhoratio;
                  state.density = rhoout;
                  state.momentum[0] = uout * rhoout;
                  state.energy = p * equation_.gamma_minus_one_inv +
                                 0.5 * rhoout * (uout * uout);
                  for (int nsp = 0; nsp < state.species.size(); nsp++) {
                    state.species[nsp] = state.species[nsp] * rhoout / rho;
                  }
                } else {
                  uout = std::sqrt(std::max(
                      0.0, u * u - 2.0 * equation_.gamma_over_gamma_minus_one *
                                       (Tratio - 1.0) * T) * Msqinv);
                  rhoout = rho * rhoratio;
                  state.density = rhoout;
                  state.momentum[0] = uout * rhoout;
                  state.energy = pflush * equation_.gamma_minus_one_inv +
                                 0.5 * rhoout * (uout * uout);
                  for (int nsp = 0; nsp < state.species.size(); nsp++) {
                    state.species[nsp] = state.species[nsp] * rhoratio;
                  }
                }
              } else {
                if (u <= 0.0) {
                  state.density = state.density * rhoratio;
                  rhoout = state.density;
                  uout = 0.0;
                  state.momentum[0] = 0.0;
                  state.energy = pflush * equation_.gamma_minus_one_inv +
                                 (ek - 0.5 * rho * u * u) * rhoratio;
                  for (int nsp = 0; nsp < state.species.size(); nsp++) {
                    state.species[nsp] = state.species[nsp] * rhoratio;
                  }
                } else {
                  state.density = state.density * rhoratio;
                  state.momentum[0] = state.momentum[0] * rhoratio;
                  rhoout = state.density;
                  uout = state.momentum[0] / rhoout;
                  state.energy =
                      pflush * equation_.gamma_minus_one_inv + ek * rhoratio;
                  for (int nsp = 0; nsp < state.species.size(); nsp++) {
                    state.species[nsp] = state.species[nsp] * rhoratio;
                  }
                }
              }
              equation_.CompleteFromCons(state, state);
              Store(states, state, dest);
            });
      }
    }
  });
}

} // namespace fub::amrex