// Copyright (c) 2020 Maikel Nadolski
// Copyright (c) 2020 Christian Zenker
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

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/ForEach.hpp"
#include "fub/equations/perfect_gas_mix/IgnitionDelayKinetics.hpp"

namespace fub::perfect_gas_mix {

template <int Rank>
IgnitionDelayKinetics<Rank>::IgnitionDelayKinetics(
    const PerfectGasMix<Rank>& eq)
    : equation_{eq}, state_{Complete<PerfectGasMix<Rank>>(eq)}, kinetic_state_{KineticState<PerfectGasMix<Rank>>(eq)} {}

template <int Rank> Duration IgnitionDelayKinetics<Rank>::ComputeStableDt(int) {
  return Duration(std::numeric_limits<double>::max());
}

template <int Rank>
Result<void, TimeStepTooLarge> IgnitionDelayKinetics<Rank>::AdvanceLevel(
    amrex::IntegratorContext& simulation_data, int level, Duration dt,
    const ::amrex::IntVect& ngrow) {
  Timer advance_timer = simulation_data.GetCounterRegistry()->get_timer(
      "IgnitionDelayKinetics::AdvanceLevel");
  ::amrex::MultiFab& data = simulation_data.GetScratch(level);
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
  for (::amrex::MFIter mfi(data, ::amrex::IntVect(8)); mfi.isValid(); ++mfi) {
    using Complete = ::fub::Complete<PerfectGasMix<Rank>>;
    View<Complete> states = amrex::MakeView<Complete>(data[mfi], *equation_,
                                                      mfi.growntilebox(ngrow));
    PerfectGasMix<Rank>& eq = *equation_;
    Complete& state = *state_;
    auto& kinetic_state = *kinetic_state_;
    ForEachIndex(Box<0>(states), [&](auto... is) {
      std::array<std::ptrdiff_t, sRank> index{is...};
      Load(state, states, index);

      FUB_ASSERT(eq.n_species >= 2);
      euler::KineticStateFromComplete(eq, kinetic_state, state);
      Array<double, -1, 1>& X = kinetic_state.mole_fractions;
      const double T = kinetic_state.temperature;
      FUB_ASSERT(dt > Duration(0.0));
      const double lambda =
          -std::log(options_.Yign / options_.Yinit / options_.tau);
      const double rate =
          lambda +
          std::pow(T, options_.R_B) *
              std::exp(options_.R_C * (1.0 - options_.R_E / T)) +
          options_.R2_A * std::exp(options_.R2_C * (1.0 - options_.R2_E / T));
      const double Fdiff = (1.0 - std::exp(-rate * dt.count())) * X[0];
      X[0] = X[0] - Fdiff;
      X[1] = std::clamp(X[1] + Fdiff, 0.0, 1.0);

      const double rate2 = std::exp(options_.C * (1.0 - options_.E / T));
      const int activator = X[0] <= options_.Yign * (1.0 - X[2]);
      
      const double coeff = activator * (1.0 - std::exp(-rate2 * dt.count()));
      FUB_ASSERT(0.0 <= coeff && coeff <= 1.0);
      const double FRdiff1 = coeff * X[1];
      FUB_ASSERT(FRdiff1 <= X[1]);
      X[1] = X[1] - FRdiff1;
      FUB_ASSERT(X[1] >= 0.0);
      X[2] = std::clamp(X[2] + FRdiff1, 0.0, 1.0);
      
      const double FRdiff0 = coeff * X[0];
      X[0] = X[0] - FRdiff0;
      FUB_ASSERT(X[0] >= 0.0);
      X[2] = std::clamp(X[2] + FRdiff0, 0.0, 1.0);

      const double T_new = T + options_.Tdiff * (FRdiff0 + FRdiff1);
      kinetic_state.temperature = T_new;

      const Array<double, Rank, 1> velocity = euler::Velocity(eq, state);
      
      euler::CompleteFromKineticState(eq, state, kinetic_state, velocity);
      Store(states, state, index);
    });
  }

  if (level == 0) {
    simulation_data.FillGhostLayerSingleLevel(level);
  } else {
    simulation_data.FillGhostLayerTwoLevels(level, level -1);
  }
  return boost::outcome_v2::success();
}

template class IgnitionDelayKinetics<1>;
template class IgnitionDelayKinetics<2>;
template class IgnitionDelayKinetics<3>;

} // namespace fub::perfect_gas_mix
