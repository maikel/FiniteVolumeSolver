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

#include <boost/numeric/odeint.hpp>

namespace fub::amrex {

AxiSymmetricSourceTerm::AxiSymmetricSourceTerm(const IdealGasMix<2>& eq)
    : equation_(eq) {}

Duration AxiSymmetricSourceTerm::ComputeStableDt(int /* level */) {
  return Duration(std::numeric_limits<double>::max());
}

Result<void, TimeStepTooLarge> AxiSymmetricSourceTerm::AdvanceLevel(
    cutcell::IntegratorContext& simulation_data, int level,
    Duration time_step_size, const ::amrex::IntVect&) {
  ::amrex::MultiFab& data = simulation_data.GetScratch(level);
  const ::amrex::Geometry& geom = simulation_data.GetGeometry(level);
  auto&& cutcell_data = simulation_data.GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = cutcell_data.getVolFrac();
  const double dt = time_step_size.count();
  Complete<IdealGasMix<2>> state(equation_.Get());
  Conservative<IdealGasMix<2>> cons(equation_.Get());
  ForEachFab(execution::seq, data, [&](const ::amrex::MFIter& mfi) {
    ::amrex::Box box = mfi.growntilebox();
    ::amrex::FArrayBox& fab = data[mfi];
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    IdealGasMix<2>& equation = equation_.Get();
    auto states = MakeView<Complete<IdealGasMix<2>>>(fab, equation, box);
    FlameMasterReactor& reactor = equation.GetReactor();
    const int ncons = simulation_data.GetFluxes(0, Direction::X).nComp();
    std::vector<double> state_buffer(ncons);
    ForEachIndex(Box<0>(states), [&](auto i, auto j) {
      static constexpr int i_radial = 1;
      static constexpr int i_axial = 0;
      Index<2> ij{i, j};
      if (alpha({AMREX_D_DECL(int(i), int(j), 0)}) > 0.0) {
        const double r = geom.CellCenter(ij[i_radial], i_radial);
        if (r > 0.0) {
          Load(cons, AsCons(states), ij);
          CopyToBuffer(span{state_buffer}, cons);
          const double rho = cons.density;
          const double u = cons.momentum[i_radial] / rho;
          const double step_size = std::min(dt, r / u);
          boost::numeric::odeint::runge_kutta4<std::vector<double>> stepper;
          boost::numeric::odeint::integrate_const(
              stepper,
              [dt, r, &cons, &state, &equation,
               &reactor](const std::vector<double>& U,
                         std::vector<double>& dUdt, double) {
                CopyFromBuffer(cons, U);
                equation.CompleteFromCons(state, cons);
                equation.SetReactorStateFromComplete(state);
                const double rho = state.density;
                const double rhoH = state.energy + state.pressure;
                const double u = state.momentum[i_radial] / rho;
                const double v = state.momentum[i_axial] / rho;
                dUdt[0] = -rho * u / r;
                dUdt[1 + i_radial] = -rho *  u * u / r;
                dUdt[1 + i_axial] = -rho * u * v / r;
                dUdt[3] = -rhoH * u / r;
                for (std::size_t s = 4; s < dUdt.size(); ++s) {
                  dUdt[s] = -rho * U[s] * u / r;
                }
              },
              state_buffer, 0.0, dt, step_size);
          CopyFromBuffer(cons, state_buffer);
          if (cons.density <= 0.0) {
            throw std::runtime_error(
                "ODE integration for the axi-symmetric source term failed.");
          }
          equation.CompleteFromCons(state, cons);
          FUB_ASSERT(state.pressure > 0.0);
          FUB_ASSERT(state.speed_of_sound > 0.0);
          Store(states, state, ij);
        }
      }
    });
  });
  return boost::outcome_v2::success();
}

} // namespace fub::amrex
