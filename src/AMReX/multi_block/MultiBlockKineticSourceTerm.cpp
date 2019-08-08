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

#include "fub/AMReX/multi_block/MultiBlockKineticSourceTerm.hpp"

namespace fub::amrex {

MultiBlockKineticSouceTerm::MultiBlockKineticSouceTerm(
    const IdealGasMix<1>& equation,
    std::shared_ptr<MultiBlockGriddingAlgorithm> gridding)
    : source_terms_() {
      source_terms_.reserve(static_cast<std::size_t>(gridding->GetTubes().size()));
  std::transform(gridding->GetTubes().begin(), gridding->GetTubes().end(),
                 std::back_inserter(source_terms_),
                 [&equation](const std::shared_ptr<GriddingAlgorithm>& grid) {
                   return ideal_gas::KineticSourceTerm<1>(equation, grid);
                 });
}

void MultiBlockKineticSouceTerm::ResetHierarchyConfiguration(
    std::shared_ptr<MultiBlockGriddingAlgorithm> gridding) {
  const std::shared_ptr<GriddingAlgorithm>* tube = gridding->GetTubes().begin();
  for (ideal_gas::KineticSourceTerm<1>& source : source_terms_) {
    source.ResetHierarchyConfiguration(*tube++);
  }
}

Duration MultiBlockKineticSouceTerm::ComputeStableDt() {
  return std::accumulate(source_terms_.begin(), source_terms_.end(),
                     Duration(std::numeric_limits<double>::infinity()),
                     [](Duration dt, ideal_gas::KineticSourceTerm<1>& source) {
                       return std::min(dt, source.ComputeStableDt());
                     });
}

Result<void, TimeStepTooLarge>
MultiBlockKineticSouceTerm::AdvanceHierarchy(Duration dt) {
  for (ideal_gas::KineticSourceTerm<1>& source : source_terms_) {
    Result<void, TimeStepTooLarge> result = source.AdvanceHierarchy(dt);
    if (!result) {
      return result;
    }
  }
  return boost::outcome_v2::success();
}

Duration MultiBlockKineticSouceTerm::GetTimePoint() const {
  return source_terms_[0].GetTimePoint();
}

std::ptrdiff_t MultiBlockKineticSouceTerm::GetCycles() const {
  return source_terms_[0].GetCycles();
}

} // namespace fub::amrex
