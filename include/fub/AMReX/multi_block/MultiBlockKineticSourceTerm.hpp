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

#ifndef FUB_AMREX_MULTI_BLOCK_KINETIC_SOURCE_TERM_HPP
#define FUB_AMREX_MULTI_BLOCK_KINETIC_SOURCE_TERM_HPP

#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp"
#include "fub/equations/ideal_gas_mix/KineticSourceTerm.hpp"

#include <vector>

namespace fub::amrex {

class MultiBlockKineticSouceTerm {
public:
    static constexpr int Rank = 1;

  MultiBlockKineticSouceTerm(
      const IdealGasMix<1>& equation,
      std::shared_ptr<MultiBlockGriddingAlgorithm> gridding);

  void ResetHierarchyConfiguration(
      std::shared_ptr<MultiBlockGriddingAlgorithm> gridding);

  Duration ComputeStableDt();

  Result<void, TimeStepTooLarge> AdvanceHierarchy(Duration dt);
  Result<void, TimeStepTooLarge> AdvanceLevel(int level, Duration dt);

  Duration GetTimePoint() const;

  std::ptrdiff_t GetCycles() const;

private:
  std::vector<ideal_gas::KineticSourceTerm<1>> source_terms_;
};

} // namespace fub::amrex

#endif
