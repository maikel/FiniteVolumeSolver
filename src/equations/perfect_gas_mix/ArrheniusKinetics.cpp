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

#include "fub/equations/perfect_gas_mix/ArrheniusKinetics.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/ForEach.hpp"
#include "fub/StateRow.hpp"

#include <Eigen/Dense>

namespace fub::perfect_gas_mix {

template <int Rank>
ArrheniusKinetics<Rank>::ArrheniusKinetics(
    const PerfectGasMix<Rank>& eq)
    : equation_{eq}, state_{Complete<PerfectGasMix<Rank>>(eq)},
      kinetic_state_{KineticState<PerfectGasMix<Rank>>(eq)} {}

template <int Rank> Duration ArrheniusKinetics<Rank>::ComputeStableDt(int) {
  return Duration(std::numeric_limits<double>::max());
}

/* namespace {
template <typename Equation>
void Advance_(Equation& eq, CompleteArray<Equation>& state,
              KineticStateArray<Equation>& kinetic_state, Duration dt,
              const ArrheniusKineticsOptions& options) {
  FUB_ASSERT(eq.n_species >= 2);
  euler::KineticStateFromComplete(eq, kinetic_state, state);
  ArrayXd& X = kinetic_state.mole_fractions;
  FUB_ASSERT(X.colwise().sum().isApproxToConstant(1.0));
  const Array1d T = kinetic_state.temperature;
  FUB_ASSERT(dt > Duration(0.0));
  const Array1d lambda =
      Array1d::Constant(-std::log(options.Yign / options.Yinit / options.tau));
  const Array1d rate =
      lambda +
      Eigen::pow(T, options.R_B) *
          Eigen::exp(options.R_C * (1.0 - options.R_E / T)) +
      options.R2_A * Eigen::exp(options.R2_C * (1.0 - options.R2_E / T));
  const Array1d Fdiff = (1.0 - Eigen::exp(-rate * dt.count())) * X.row(0);
  const Array1d Xrow0a = X.row(0) - Fdiff;
  const Array1d Xrow1a = (X.row(1) + Fdiff).max(0.0).min(1.0);
  X.row(0) = Xrow0a;
  X.row(1) = Xrow1a;
  const Array1d rate2 = Eigen::exp(options.C * (1.0 - options.E / T));
  const MaskArray activator = X.row(0) <= (options.Yign * (1.0 - X.row(2)));
  const Array1d coeff =
      activator.select((1.0 - Eigen::exp(-rate2 * dt.count())), 0.0);
  FUB_ASSERT((0.0 <= coeff && coeff <= 1).all());
  const Array1d FRdiff1 = coeff * X.row(1);
  const Array1d Xrow1b = X.row(1) - FRdiff1;
  const Array1d Xrow2b = (X.row(2) + FRdiff1).max(0.0).min(1.0);
  X.row(1) = Xrow1b;
  X.row(2) = Xrow2b;
  const Array1d FRdiff0 = coeff * X.row(0);
  const Array1d Xrow0c = X.row(0) - FRdiff0;
  const Array1d Xrow2c = (X.row(2) + FRdiff0).max(0.0).min(1.0);
  X.row(0) = Xrow0c;
  X.row(2) = Xrow2c;
  FUB_ASSERT(X.colwise().sum().isApproxToConstant(1.0));
  const Array1d T_new = T + options.Tdiff * (FRdiff0 + FRdiff1);
  kinetic_state.temperature = T_new;
  const Array<double, Equation::Rank()> velocity = euler::Velocity(eq, state);
  euler::CompleteFromKineticState(eq, state, kinetic_state, velocity);
}
} // namespace */

/* template <typename Equation> struct ArrheniusKinetics_Rows {
  Equation* equation_;
  CompleteArray<Equation>* state_;
  KineticStateArray<Equation>* kinetic_state_;
  const ArrheniusKineticsOptions* options_;
  Duration dt_;

  void operator()(const Row<Complete<Equation>>& row) const {
    ViewPointer in = AsConst(Begin(row));
    ViewPointer end = AsConst(End(row));
    ViewPointer<Complete<Equation>> out = Begin(row);
    Equation& eq = *equation_;
    CompleteArray<Equation>& complete = *state_;
    KineticStateArray<Equation>& kin = *kinetic_state_;
    int n = static_cast<int>(get<0>(end) - get<0>(in));
    while (n >= kDefaultChunkSize) {
      Load(complete, in);
      Advance_(eq, complete, kin, dt_, *options_);
      Store(out, complete);
      Advance(in, kDefaultChunkSize);
      Advance(out, kDefaultChunkSize);
      n = static_cast<int>(get<0>(end) - get<0>(in));
    }
    LoadN(complete, in, n);
    Advance_(eq, complete, kin, dt_, *options_);
    StoreN(out, complete, n);
  }
}; */

template <int Rank>
Result<void, TimeStepTooLarge> ArrheniusKinetics<Rank>::AdvanceLevel(
    amrex::IntegratorContext& simulation_data, int level, Duration dt,
    const ::amrex::IntVect& ngrow) {
  Timer advance_timer = simulation_data.GetCounterRegistry()->get_timer(
      "ArrheniusKinetics::AdvanceLevel");
  ::amrex::MultiFab& data = simulation_data.GetScratch(level);
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
  for (::amrex::MFIter mfi(data, ::amrex::IntVect(8)); mfi.isValid(); ++mfi) {
    using Complete = ::fub::Complete<PerfectGasMix<Rank>>;
    View<Complete> states = amrex::MakeView<Complete>(data[mfi], *equation_,
                                                      mfi.growntilebox(ngrow));
    PerfectGasMix<Rank>& eq = *equation_;
    auto& state = *state_;
    auto& kinetic_state = *kinetic_state_;
    /* ForEachRow(std::tuple{states},
               ArrheniusKinetics_Rows<PerfectGasMix<Rank>>{
                   &eq, &state, &kinetic_state, &options, dt}); */

      ForEachIndex(Box<0>(states), [&](auto... is) {
        std::array<std::ptrdiff_t, sRank> index{is...};
        Load(state, states, index);

        FUB_ASSERT(eq.n_species >= 1);
        euler::KineticStateFromComplete(eq, kinetic_state, state);
        Array<double, -1, 1>& X = kinetic_state.mole_fractions;
        FUB_ASSERT(X.colwise().sum().isApproxToConstant(1.0));
        const double T = kinetic_state.temperature;
        FUB_ASSERT(dt > Duration(0.0));
        
        const int activator = (T >= options.T_switch);
        const double rate = std::exp( options.EA * ( 1.0 -1.0 / T) );
        const double Xnew = X[0] * std::exp( - options.B * dt.count() * rate );
        
        const double dX = activator * (Xnew - X[0]);
        X[0] = X[0] + dX;
        X[1] = X[1] - dX;
        FUB_ASSERT(X.colwise().sum().isApproxToConstant(1.0));

        const Array<double, Rank, 1> velocity = euler::Velocity(eq, state);
        euler::CompleteFromKineticState(eq, state, kinetic_state, velocity);

        state.energy += - options.Q * state.density * dX;

        fub::CompleteFromCons(eq, state, state);

        Store(states, state, index);
      });
  }

  if (level == 0) {
    simulation_data.FillGhostLayerSingleLevel(level);
  } else {
    simulation_data.FillGhostLayerTwoLevels(level, level - 1);
  }
  return boost::outcome_v2::success();
}

template class ArrheniusKinetics<1>;
template class ArrheniusKinetics<2>;
template class ArrheniusKinetics<3>;

} // namespace fub::perfect_gas_mix
