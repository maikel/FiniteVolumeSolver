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

namespace fub::ideal_gas {

template <int Rank>
void KineticSourceTerm<Rank>::ResetHierarchyConfiguration(
    std::shared_ptr<amrex::GriddingAlgorithm>&& gridding) {
  gridding_ = std::move(gridding);
}

template <int Rank>
void KineticSourceTerm<Rank>::ResetHierarchyConfiguration(
    const std::shared_ptr<amrex::GriddingAlgorithm>& gridding) {
  gridding_ = gridding;
}

template <int Rank>
KineticSourceTerm<Rank>::KineticSourceTerm(
    const IdealGasMix<Rank>& eq,
    std::shared_ptr<amrex::GriddingAlgorithm> gridding)
    : equation_{eq}, state_{Complete<IdealGasMix<Rank>>(eq)},
      gridding_{std::move(gridding)} {}

template <int Rank> Duration KineticSourceTerm<Rank>::ComputeStableDt(int) {
  return Duration(std::numeric_limits<double>::infinity());
}

template <int Rank>
Result<void, TimeStepTooLarge>
KineticSourceTerm<Rank>::AdvanceLevel(int level, Duration dt) {
  ::amrex::MultiFab& data =
      gridding_->GetPatchHierarchy().GetPatchLevel(level).data;
  fub::amrex::ForEachFab(
      execution::openmp, data, [&](const ::amrex::MFIter& mfi) {
        using Complete = ::fub::Complete<IdealGasMix<Rank>>;
        View<Complete> states =
            amrex::MakeView<Complete>(data[mfi], *equation_, mfi.tilebox());
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
      });
  return boost::outcome_v2::success();
}

template <int Rank> Duration KineticSourceTerm<Rank>::GetTimePoint() const {
  return GetPatchHierarchy().GetTimePoint();
}

template <int Rank> std::ptrdiff_t KineticSourceTerm<Rank>::GetCycles() const {
  return GetPatchHierarchy().GetCycles();
}

template <int Rank>
const amrex::PatchHierarchy&
KineticSourceTerm<Rank>::GetPatchHierarchy() const {
  return gridding_->GetPatchHierarchy();
}

template <int Rank>
amrex::PatchHierarchy& KineticSourceTerm<Rank>::GetPatchHierarchy() {
  return gridding_->GetPatchHierarchy();
}

template class KineticSourceTerm<1>;
template class KineticSourceTerm<2>;
template class KineticSourceTerm<3>;

} // namespace fub::ideal_gas
