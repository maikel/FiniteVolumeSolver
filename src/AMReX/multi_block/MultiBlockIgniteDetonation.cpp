#include "fub/AMReX/multi_block/MultiBlockIgniteDetonation.hpp"

namespace fub::amrex {

MultiBlockIgniteDetonation::MultiBlockIgniteDetonation(
    const fub::IdealGasMix<1>& equation,
    const std::shared_ptr<MultiBlockGriddingAlgorithm>& grid,
    const IgniteDetonationOptions& opts) {
  span<const std::shared_ptr<GriddingAlgorithm>> tubes = grid->GetTubes();
  source_terms_.reserve(static_cast<std::size_t>(tubes.size()));
  std::transform(tubes.begin(), tubes.end(), std::back_inserter(source_terms_),
                 [&](const std::shared_ptr<GriddingAlgorithm>& tube) {
                   return IgniteDetonation(equation, tube, opts);
                 });
}

Duration MultiBlockIgniteDetonation::ComputeStableDt() noexcept {
  return Duration(std::numeric_limits<double>::infinity());
}

void MultiBlockIgniteDetonation::ResetHierarchyConfiguration(
    const std::shared_ptr<MultiBlockGriddingAlgorithm>& grid) {
  const std::shared_ptr<GriddingAlgorithm>* tube = grid->GetTubes().begin();
  for (IgniteDetonation& source : source_terms_) {
    source.ResetHierarchyConfiguration(*tube++);
  }
}

Result<void, TimeStepTooLarge>
MultiBlockIgniteDetonation::AdvanceLevel(int level, Duration dt) {
  for (IgniteDetonation& source : source_terms_) {
    if (level < source.GetPatchHierarchy().GetNumberOfLevels()) {
      Result<void, TimeStepTooLarge> result = source.AdvanceLevel(level, dt);
      if (!result) {
        return result;
      }
    }
  }
  return boost::outcome_v2::success();
}

}