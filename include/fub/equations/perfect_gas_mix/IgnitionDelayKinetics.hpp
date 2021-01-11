// Copyright (c) 2020 Maikel Nadolski
// Copyright (c) 2020 Christian Zenker
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

#ifndef FUB_PERFECT_GAS_MIX_IGNITION_DELAY_KINETICS_HPP
#define FUB_PERFECT_GAS_MIX_IGNITION_DELAY_KINETICS_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/equations/PerfectGasMix.hpp"
#include "fub/ext/omp.hpp"
#include "fub/ext/outcome.hpp"

#include <memory>

namespace fub::perfect_gas_mix {

struct IgnitionDelayKineticsOptions {
 double Yign{0.2};
 double Yinit{1.0};
 double tau{1.0};
 double R_B{0.};
 double R_C{70.};
 double R_E{1.2};
 double R2_A{0.};
 double R2_C{5.};
 double R2_E{1.5};
 double E{0.5};
 double C{7.};
 double Tdiff{2.7};
};

/// \ingroup LevelIntegrator
template <int R> class IgnitionDelayKinetics {
public:
  static constexpr int Rank = R;
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);

  IgnitionDelayKinetics(const PerfectGasMix<Rank>& eq);

  Duration ComputeStableDt(int level);

  Result<void, TimeStepTooLarge>
  AdvanceLevel(amrex::IntegratorContext& simulation_data, int level,
               Duration dt,
               const ::amrex::IntVect& ngrow = ::amrex::IntVect(0));
  
  IgnitionDelayKineticsOptions options;

private:
  OmpLocal<PerfectGasMix<Rank>> equation_;
  OmpLocal<CompleteArray<PerfectGasMix<Rank>>> state_;
  OmpLocal<KineticStateArray<PerfectGasMix<Rank>>> kinetic_state_;
};

extern template class IgnitionDelayKinetics<1>;
extern template class IgnitionDelayKinetics<2>;
extern template class IgnitionDelayKinetics<3>;

} // namespace fub::perfect_gas_mix

#endif
