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

Duration MultiBlockIgniteDetonation::ComputeStableDt([
    [maybe_unused]] int level) noexcept {
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

std::vector<Duration>
MultiBlockIgniteDetonation::GetLastIgnitionTimePoints() const {
  std::vector<Duration> times{};
  times.reserve(source_terms_.size());
  std::transform(source_terms_.begin(), source_terms_.end(),
                 std::back_inserter(times), [](const IgniteDetonation& ign) {
                   return ign.GetLastIgnitionTimePoint(0);
                 });
  return times;
}

void MultiBlockIgniteDetonation::SetLastIgnitionTimePoints(
    span<const Duration> timepoints) {
  int k = 0;
  for (Duration t_ign : timepoints) {
    const int nlevel = source_terms_[k].GetPatchHierarchy().GetNumberOfLevels();
    for (int ilvl = 0; ilvl < nlevel; ++ilvl) {
      source_terms_[k].SetLastIgnitionTimePoint(ilvl, t_ign);
    }
    k += 1;
  }
}

} // namespace fub::amrex
