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

#ifndef FUB_PERFECT_GAS_MIX_ARRHENIUS_KINETICS_HPP
#define FUB_PERFECT_GAS_MIX_ARRHENIUS_KINETICS_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/equations/PerfectGasMix.hpp"
#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"
#include "fub/ext/omp.hpp"
#include "fub/ext/outcome.hpp"

#include <cmath>
#include <memory>

namespace fub::perfect_gas_mix {

struct ArrheniusKineticsOptions {
  ArrheniusKineticsOptions() = default; ///< default constructor
  ArrheniusKineticsOptions(
      const ProgramOptions& options); ///< constructor with options

  /// \brief simple Print function for the options
  ///
  /// \param log the logger for the output
  template <typename Log> void Print(Log& log);

  /* dimensionless thermochemical parameters */
  /* think of T1 as compressor exit temperature */
  // default values taken from SEC_C/SEC_C/Input/userdata_Combustion.c

  // Test
  double Q{3.0 / 0.4};   ///< heat of reaction
  double EA{11.0};       ///< activation energy
  double B{0.05 / 1.25}; ///< coefficient of reaction
  /* Reference: 0.11/1.25; 0.05/1.25 */
  double T_switch{1.05}; ///< switch temperature
};

/// \ingroup LevelIntegrator
template <int R> class ArrheniusKinetics {
public:
  static constexpr int Rank = R;
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);

  ArrheniusKinetics(const PerfectGasMix<Rank>& eq);

  Duration ComputeStableDt(int level);

  Result<void, TimeStepTooLarge>
  AdvanceLevel(amrex::IntegratorContext& simulation_data, int level,
               Duration dt,
               const ::amrex::IntVect& ngrow = ::amrex::IntVect(0));

  ArrheniusKineticsOptions options;

private:
  OmpLocal<PerfectGasMix<Rank>> equation_;
  OmpLocal<CompleteArray<PerfectGasMix<Rank>>> state_;
  OmpLocal<KineticStateArray<PerfectGasMix<Rank>>> kinetic_state_;
};

template <typename Log> void ArrheniusKineticsOptions::Print(Log& log) {
  BOOST_LOG(log) << fmt::format(" - Q = {}", Q);
  BOOST_LOG(log) << fmt::format(" - EA = {}", EA);
  BOOST_LOG(log) << fmt::format(" - B = {}", B);
  BOOST_LOG(log) << fmt::format(" - T_switch = {}", T_switch);
}

extern template class ArrheniusKinetics<1>;
extern template class ArrheniusKinetics<2>;
extern template class ArrheniusKinetics<3>;

} // namespace fub::perfect_gas_mix

#endif
