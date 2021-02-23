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

ArrheniusKineticsOptions::ArrheniusKineticsOptions(
    const ProgramOptions& options) {
  Q = GetOptionOr(options, "Q", Q);
  EA = GetOptionOr(options, "EA", EA);
  B = GetOptionOr(options, "B", B);
  T_switch = GetOptionOr(options, "T_switch", T_switch);
}

template <int Rank>
ArrheniusKinetics<Rank>::ArrheniusKinetics(const PerfectGasMix<Rank>& eq,
                                           const ArrheniusKineticsOptions& opts)
    : options{opts}, equation_{eq}, state_{CompleteArray<PerfectGasMix<Rank>>(
                                        eq)},
      kinetic_state_{KineticStateArray<PerfectGasMix<Rank>>(eq)} {}

template <int Rank>
Duration
ArrheniusKinetics<Rank>::ComputeStableDt(const amrex::IntegratorContext&, int) {
  return Duration(std::numeric_limits<double>::max());
}

namespace {
template <typename Equation>
void Advance_(Equation& eq, CompleteArray<Equation>& state,
              KineticStateArray<Equation>& kinetic_state, Duration dt,
              const ArrheniusKineticsOptions& options) {
  FUB_ASSERT(eq.n_species >= 1);
  euler::KineticStateFromComplete(eq, kinetic_state, state);
  ArrayXd& X = kinetic_state.mole_fractions;
  FUB_ASSERT(X.colwise().sum().isApproxToConstant(1.0));
  const Array1d T = kinetic_state.temperature;
  FUB_ASSERT(dt > Duration(0.0));

  const MaskArray activator = (T >= options.T_switch);
  const Array1d rate = Eigen::exp(options.EA * (1.0 - 1.0 / T));
  const Array1d Xnew = X.row(0) * Eigen::exp(-options.B * dt.count() * rate);

  const Array1d dX = activator.select((Xnew - X.row(0)), 0.0);
  X.row(0) = X.row(0) + dX;
  X.row(1) = X.row(1) - dX;
  FUB_ASSERT(X.colwise().sum().isApproxToConstant(1.0));

  const Array<double, Equation::Rank()> velocity = euler::Velocity(eq, state);
  euler::CompleteFromKineticState(eq, state, kinetic_state, velocity);

  state.energy += -options.Q * state.density * dX;

  fub::CompleteFromCons(eq, state, state);
}
} // namespace

template <typename Equation> struct ArrheniusKinetics_Rows {
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
};

template <int Rank>
Result<void, TimeStepTooLarge>
ArrheniusKinetics<Rank>::AdvanceLevel(amrex::IntegratorContext& simulation_data,
                                      int level, Duration dt,
                                      const ::amrex::IntVect& ) {
  Timer advance_timer = simulation_data.GetCounterRegistry()->get_timer(
      "ArrheniusKinetics::AdvanceLevel");
  ::amrex::MultiFab& data = simulation_data.GetScratch(level);

  amrex::ForEachFab(execution::openmp, data, [&](const ::amrex::MFIter& mfi) {
    using Complete = ::fub::Complete<PerfectGasMix<Rank>>;
    View<Complete> states = amrex::MakeView<Complete>(data[mfi], *equation_,
                                                      mfi.growntilebox());
    PerfectGasMix<Rank>& eq = *equation_;
    auto& state = *state_;
    auto& kinetic_state = *kinetic_state_;
    ForEachRow(std::tuple{states},
               ArrheniusKinetics_Rows<PerfectGasMix<Rank>>{
                   &eq, &state, &kinetic_state, &options, dt});
  });

  for (Direction dir : {Direction::X, Direction::Y, Direction::Z}) {
    simulation_data.ApplyBoundaryCondition(level, dir);
  }
  return boost::outcome_v2::success();
}

template class ArrheniusKinetics<1>;
template class ArrheniusKinetics<2>;
template class ArrheniusKinetics<3>;

} // namespace fub::perfect_gas_mix