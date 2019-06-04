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

#ifndef FUB_IDEAL_GAS_MIX_KINETIC_SOURCE_TERM_HPP
#define FUB_IDEAL_GAS_MIX_KINETIC_SOURCE_TERM_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/equations/IdealGasMix.hpp"
#include "fub/ext/outcome.hpp"

#include <functional>
#include <optional>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fub {
namespace ideal_gas {

template <int R> class KineticSourceTerm {
public:
  static constexpr int Rank = R;
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);

  KineticSourceTerm(const IdealGasMix<Rank>& eq,
                    std::shared_ptr<amrex::GriddingAlgorithm> gridding);

  void ResetHierarchyConfiguration(
      std::shared_ptr<amrex::GriddingAlgorithm>&& griding);

  void ResetHierarchyConfiguration(
      const std::shared_ptr<amrex::GriddingAlgorithm>& griding);

  Duration ComputeStableDt();
  Duration ComputeStableDt(amrex::PatchHandle patch);

  Result<void, TimeStepTooLarge> AdvanceHierarchy(Duration dt);
  Result<void, TimeStepTooLarge> AdvancePatch(amrex::PatchHandle patch,
                                              Duration dt);

  Duration GetTimePoint() const;

  std::ptrdiff_t GetCycles() const;

  amrex::PatchHierarchy& GetPatchHierarchy();
  const amrex::PatchHierarchy& GetPatchHierarchy() const;

private:
  IdealGasMix<Rank> equation_;
  std::shared_ptr<amrex::GriddingAlgorithm> gridding_;
  Complete<IdealGasMix<Rank>> state_{equation_};
};

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
    : equation_{eq}, gridding_{std::move(gridding)} {}

template <int Rank> Duration KineticSourceTerm<Rank>::ComputeStableDt() {
  const int nlevels = GetPatchHierarchy().GetNumberOfLevels();
  double dt(std::numeric_limits<double>::infinity());
  for (int level = 0; level < nlevels; ++level) {
    dt = std::min(
        dt, GetPatchHierarchy().Minimum(
                level, [solver = *this](amrex::PatchHandle patch) mutable {
                  return solver.ComputeStableDt(patch).count();
                }));
  }
  return Duration(dt);
}

template <int Rank>
Duration KineticSourceTerm<Rank>::ComputeStableDt(amrex::PatchHandle patch) {
  const amrex::PatchLevel& level =
      GetPatchHierarchy().GetPatchLevel(patch.level);
  const IndexBox<Rank> tilebox =
      amrex::AsIndexBox<Rank>(patch.iterator->tilebox());
  const View<const Complete<IdealGasMix<Rank>>> states = Subview(
      amrex::MakeView<BasicView<const Complete<IdealGasMix<Rank>>>>(
          amrex::MakePatchDataView(level.data[*patch.iterator]), equation_),
      tilebox);
  FlameMasterReactor& reactor = equation_.GetReactor();
  double dt_chem = std::numeric_limits<double>::infinity();
  ForEachIndex(Box<0>(states), [&, this](auto... is) {
    std::array<std::ptrdiff_t, sRank> index{is...};
    Load(state_, states, index);
    equation_.SetReactorStateFromComplete(state_);
    span<const double> production_rates = reactor.GetProductionRates();
    span<const double> X = reactor.GetMoleFractions();
    dt_chem = std::inner_product(
        X.begin(), X.end(), production_rates.begin(), dt_chem,
        [](double t1, double t2) { return std::min(t1, std::abs(t2)); },
        [](double /* X */, double /* dXdt */) {
          //          if (dXdt < 0.0) {
          //            FUB_ASSERT(X > 0.0);
          //            return 10000.0 * X / -dXdt;
          //          }
          return std::numeric_limits<double>::infinity();
        });
  });
  return Duration(dt_chem);
}

template <int Rank>
Result<void, TimeStepTooLarge>
KineticSourceTerm<Rank>::AdvanceHierarchy(Duration dt) {
  const int nlevels = GetPatchHierarchy().GetNumberOfLevels();
  for (int level = 0; level < nlevels; ++level) {
#ifdef _OPENMP
    const std::size_t n_threads =
        static_cast<std::size_t>(::omp_get_max_threads());
    std::vector<std::optional<Result<void, TimeStepTooLarge>>> results(
        n_threads);
    GetPatchHierarchy().ForEachPatch(
        execution::openmp, level,
        [&results, dt, solver = *this](amrex::PatchHandle patch) mutable {
          const std::size_t this_thread_num =
              static_cast<std::size_t>(::omp_get_thread_num());
          if (!results[this_thread_num] || *results[this_thread_num]) {
            results[this_thread_num] = solver.AdvancePatch(patch, dt);
          }
        });
    for (std::size_t thread_num = 0; thread_num < n_threads; ++thread_num) {
      if (results[thread_num] && !(*results[thread_num])) {
        return *results[thread_num];
      }
    }
#else
    std::optional<Result<void, TimeStepTooLarge>> result{};
    GetPatchHierarchy().ForEachPatch(
        level, [&result, dt, solver = *this](amrex::PatchHandle patch) mutable {
          if (!result || *result) {
            result = solver.AdvancePatch(patch, dt);
          }
        });
    if (result && !(*result)) {
      return *result;
    }
#endif
  }
  return boost::outcome_v2::success();
}

template <int Rank>
Result<void, TimeStepTooLarge>
KineticSourceTerm<Rank>::AdvancePatch(amrex::PatchHandle patch, Duration dt) {
  amrex::PatchLevel& level = GetPatchHierarchy().GetPatchLevel(patch.level);
  const IndexBox<Rank> tilebox =
      amrex::AsIndexBox<Rank>(patch.iterator->tilebox());
  const View<Complete<IdealGasMix<Rank>>> states = Subview(
      amrex::MakeView<BasicView<Complete<IdealGasMix<Rank>>>>(
          amrex::MakePatchDataView(level.data[*patch.iterator]), equation_),
      tilebox);
  FlameMasterReactor& reactor = equation_.GetReactor();
  ForEachIndex(Box<0>(states), [&, this](auto... is) {
    std::array<std::ptrdiff_t, sRank> index{is...};
    Load(state_, states, index);
    equation_.SetReactorStateFromComplete(state_);
    reactor.Advance(dt.count());
    equation_.CompleteFromReactor(state_);
    Store(states, state_, index);
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

} // namespace ideal_gas
} // namespace fub

#endif
