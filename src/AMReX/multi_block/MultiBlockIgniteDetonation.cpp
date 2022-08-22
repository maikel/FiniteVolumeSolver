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
    const fub::IdealGasMix<1>& equation, std::size_t n_tubes,
    int max_number_levels, const IgniteDetonationOptions& opts)
    : source_terms_(n_tubes,
                    IgniteDetonation(equation, max_number_levels, opts)),
      max_number_levels_{max_number_levels} {}

MultiBlockIgniteDetonation::MultiBlockIgniteDetonation(
    const fub::IdealGasMix<1>& equation, std::size_t n_tubes,
    int max_number_levels, const std::vector<IgniteDetonationOptions>& opts)
    : max_number_levels_{max_number_levels} {
  source_terms_.reserve(n_tubes);
  for (std::size_t i = 0; i < n_tubes; ++i) {
    source_terms_.emplace_back(equation, max_number_levels, opts[i]);
  }
}

Duration MultiBlockIgniteDetonation::ComputeStableDt(const MultiBlockIntegratorContext2&, [
    [maybe_unused]] int level) noexcept {
  return Duration(std::numeric_limits<double>::max());
}

Result<void, TimeStepTooLarge>
MultiBlockIgniteDetonation::AdvanceLevel(MultiBlockIntegratorContext2& context,
                                         int level, Duration dt,
                                         const ::amrex::IntVect& ngrow) {
  IntegratorContext* tube = context.Tubes().begin();
  for (IgniteDetonation& source : source_terms_) {
    if (level < tube->GetPatchHierarchy().GetNumberOfLevels()) {
      Result<void, TimeStepTooLarge> result =
          source.AdvanceLevel(*tube, level, dt, ngrow);
      if (!result) {
        return result;
      }
    }
    ++tube;
  }
  return boost::outcome_v2::success();
}

} // namespace fub::amrex
