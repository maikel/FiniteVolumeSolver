// Copyright (c) 2018 Maikel Nadolski
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

#include "fub/equations/ideal_gas_mix/KineticSourceTerm.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/ForEach.hpp"

namespace fub::ideal_gas {

template <int Rank>
KineticSourceTerm<Rank>::KineticSourceTerm(const IdealGasMix<Rank>& eq)
    : equation_{eq}, state_{Complete<IdealGasMix<Rank>>(eq)} {}

template <int Rank> Duration KineticSourceTerm<Rank>::ComputeStableDt(const amrex::IntegratorContext&, int) {
  return Duration(std::numeric_limits<double>::max());
}

template <int Rank> Duration KineticSourceTerm<Rank>::ComputeStableDt(int) {
  return Duration(std::numeric_limits<double>::max());
}

template <int Rank>
Result<void, TimeStepTooLarge>
KineticSourceTerm<Rank>::AdvanceLevel(amrex::IntegratorContext& simulation_data,
                                      int level, Duration dt,
                                      const ::amrex::IntVect& ngrow) {
  Timer advance_timer = simulation_data.GetCounterRegistry()->get_timer(
      "KineticSourceTerm::AdvanceLevel");
  ::amrex::MultiFab& data = simulation_data.GetScratch(level);
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
  for (::amrex::MFIter mfi(data, ::amrex::IntVect(8)); mfi.isValid(); ++mfi) {
    using Complete = ::fub::Complete<IdealGasMix<Rank>>;
    View<Complete> states = amrex::MakeView<Complete>(data[mfi], *equation_,
                                                      mfi.growntilebox(ngrow));
    FlameMasterReactor& reactor = equation_->GetReactor();
    ForEachIndex(Box<0>(states), [&](auto... is) {
      std::array<std::ptrdiff_t, sRank> index{is...};
      Load(*state_, states, index);
      // equation_->SetReactorStateFromComplete(*state_);
      reactor.SetMassFractions(state_->species);
      reactor.SetTemperature(state_->temperature);
      reactor.SetDensity(state_->density);
      reactor.Advance(dt.count());
      Eigen::Matrix<double, Rank, 1> velocity =
          state_->momentum / state_->density;
      equation_->CompleteFromReactor(*state_, velocity);
      Store(states, *state_, index);
    });
  }
  if (level == 0) {
    simulation_data.FillGhostLayerSingleLevel(level);
  } else {
    simulation_data.FillGhostLayerTwoLevels(level, level - 1);
  }
  return boost::outcome_v2::success();
}

template class KineticSourceTerm<1>;
template class KineticSourceTerm<2>;
template class KineticSourceTerm<3>;

} // namespace fub::ideal_gas
