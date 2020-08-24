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

#include "fub/AMReX/cutcell/AxiSymmetricSourceTerm.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/equations/IdealGasMix.hpp"

namespace fub::amrex {

AxiSymmetricSourceTerm::AxiSymmetricSourceTerm(const IdealGasMix<2>& eq)
    : equation_(eq) {}

Duration AxiSymmetricSourceTerm::ComputeStableDt(int /* level */) {
  return Duration(std::numeric_limits<double>::max());
}

Result<void, TimeStepTooLarge>
AxiSymmetricSourceTerm::AdvanceLevel(cutcell::IntegratorContext& simulation_data, int level,
                              Duration time_step_size, const ::amrex::IntVect&) {
  ::amrex::MultiFab& data = simulation_data.GetScratch(level);
  const ::amrex::Geometry& geom = simulation_data.GetGeometry(level);
  auto&& cutcell_data = simulation_data.GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = cutcell_data.getVolFrac();
  const double dt = time_step_size.count();
  Complete<IdealGasMix<2>> state(equation_);
  ForEachFab(execution::seq, data, [&](const ::amrex::MFIter& mfi) {
    ::amrex::Box box = mfi.growntilebox();
    ::amrex::FArrayBox& fab = data[mfi];
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    auto states = MakeView<Complete<IdealGasMix<2>>>(fab, equation_, box);
    FlameMasterReactor& reactor = equation_.GetReactor();
    ForEachIndex(Box<0>(states), [&](auto i, auto j) {
      static constexpr int i_radial = 1;
      static constexpr int i_axial = 0;
      Index<2> ij{i, j};
      Load(state, states, ij);
      if (alpha({int(i), int(j)}) == 1.0) {
        equation_.SetReactorStateFromComplete(state);
        const double H = reactor.GetEnthalpy();
        const double r = geom.CellCenter(ij[i_radial], i_radial);
        if (r > 0.0) {
          const double rho = state.density;
          const double u = state.momentum[i_radial] / rho;
          const double v = state.momentum[i_axial] / rho;
          const double lambda = dt / r;
          state.density -= lambda * rho * u;
          FUB_ASSERT(state.density > 0.0);
          state.momentum[i_radial] -= lambda * rho * u * u;
          state.momentum[i_axial] -= lambda * rho * u * v;
          state.energy -= lambda * H * u;
          state.species -= lambda * rho * state.species * u;
          equation_.CompleteFromCons(state, state);
        }
      }
    });
  });
  return boost::outcome_v2::success();
}

} // namespace fub::amrex
