// Copyright (c) 2020 Maikel Nadolski
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

#ifndef FUB_AMREX_AXI_SYMMETRIC_SOURCE_TERM_HPP
#define FUB_AMREX_AXI_SYMMETRIC_SOURCE_TERM_HPP

#include "fub/AMReX/cutcell/IntegratorContext.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/equations/IdealGasMix.hpp"
#include "fub/ext/outcome.hpp"

namespace fub::amrex {

/// \ingroup LevelIntegrator
class AxiSymmetricSourceTerm {
public:
  static constexpr int Rank = 2;

  AxiSymmetricSourceTerm(const IdealGasMix<2>& eq);

  Duration ComputeStableDt(int level);

  Result<void, TimeStepTooLarge>
  AdvanceLevel(cutcell::IntegratorContext& simulation_data, int level,
               Duration dt, const ::amrex::IntVect& = ::amrex::IntVect(0));

private:
  OmpLocal<IdealGasMix<2>> equation_;
};

} // namespace fub::amrex

#endif