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

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/equations/IdealGasMix.hpp"
#include "fub/ext/omp.hpp"
#include "fub/ext/outcome.hpp"
#include "fub/counter/CounterRegistry.hpp"

#include <memory>

namespace fub::ideal_gas {

template <int R> class KineticSourceTerm {
public:
  static constexpr int Rank = R;
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);

  KineticSourceTerm(const IdealGasMix<Rank>& eq);

  KineticSourceTerm(const IdealGasMix<Rank>& eq,
                    std::shared_ptr<CounterRegistry> registry);

  Duration ComputeStableDt(int level);

  Result<void, TimeStepTooLarge> AdvanceLevel(amrex::IntegratorContext& simulation_data, int level, Duration dt, const ::amrex::IntVect& ngrow = ::amrex::IntVect(0));

private:
  OmpLocal<IdealGasMix<Rank>> equation_;
  OmpLocal<Complete<IdealGasMix<Rank>>> state_;
  std::shared_ptr<CounterRegistry> registry_;
};

extern template class KineticSourceTerm<1>;
extern template class KineticSourceTerm<2>;
extern template class KineticSourceTerm<3>;

} // namespace fub::ideal_gas

#endif
