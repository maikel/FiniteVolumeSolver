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

/// \file
/// \brief Ths file contains all declarations for a kinetic source term
/// for the reactive euler equations.

#ifndef FUB_PERFECT_GAS_MIX_IGNITION_DELAY_KINETICS_HPP
#define FUB_PERFECT_GAS_MIX_IGNITION_DELAY_KINETICS_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/equations/PerfectGasMix.hpp"
#include "fub/ext/omp.hpp"
#include "fub/ext/outcome.hpp"

#include <memory>

namespace fub::perfect_gas_mix {

///@{
/// \brief The options for the IgnitionDelay model.
///
/// Factors determining lambda for the radical reaction:
/// - Yign
/// - Yinit,
/// - tau.
///
/// this shouldn't need any modifying ever, since the model is scaled to
/// have ignition delay times of order unity.
///
/// Factors determining the temperature dependence in the radical buildup
/// reaction:
/// - R_B,
/// - R_C,
/// - R_E.
///
/// Shouldn't need modifying neither, but can be used to make the
/// NTC region smaller (e.g. currently |1.0 - 1.4| -> 400 K)
///
/// Factors determining the secondary temperature dependence in the radical
/// buildup reaction:
/// - R2_A
/// - R2_C,
/// - R2_E.
///
/// This one is interesting. It controls the temperature dependence of the
/// ignition delay time, as a perturbation around unity. Via R2_A = 0, this
/// can be disabled for completely homogeneous reactions.
/// Modifying R2_A has a visible impact on the pressure waves' slope in the
/// heatmaps: Larger values make the lines more even.
/// It's a good idea to choose A/C small and E large to retain the overall
/// ignition delay time behaviour.
///
/// Factors determining the heat releasing reaction:
/// - E,
/// - C,
/// -Tdiff.
///
/// This controls the exitation time. If E < .9, excitation time is
/// sufficiently small that combustion is homogeneous. For E > .9,
/// it is no more (at least not reliably).
///@}
struct IgnitionDelayKineticsOptions {
  // Factors determining lambda for the radical reaction:
  double Yign{0.2};  ///< Y_fuel below which low temperature mech. kicks in
  double Yinit{1.0}; ///< Initial Y_fuel
  double tau{1.0};   ///< Ignition delay time

  // Factors determining the temperature dependence in the radical buildup
  // reaction:
  double R_B{0.};
  double R_C{70.};
  double R_E{1.2};

  // Factors determining the secondary temperature dependence in the radical
  // buildup reaction:
  double R2_A{0.};
  double R2_C{5.};
  double R2_E{1.5};

  // Factors determining the heat releasing reaction:

  /// T_ign for Arrhenius kinetics, makes excitation time of
  /// order 10^-2, i.e. ~10µs, if set to .5 (or .9)
  double E{0.5};
  /// Factor that controls how sharp the activation energy barrier
  /// is. Many global reactions have A = 10^10, dedimensionalize
  /// by multiplication with τ this gives C=16 We use 7 because
  /// this models excitation times only, and excitation times do
  /// not strongly depend on temperature.
  double C{7.};
  /// Heat release. T0 ≈ 1000 K, so this leads to T_0 + 2500 K
  double Tdiff{2.7};
};

/// \ingroup LevelIntegrator
///
/// \brief A kinetics implementation with two species and NTC/fixed tau.
///
/// This is a non-dimensionalized model with two reactions:
///  - Fuel -> Fuel radicals    -- temperature independent
///  - Fuel radicals -> Product -- Arrhenius, Heaviside-dependent on Fuel
///  existence
///
/// Further Details can be found in \cite Berndt2017 or in \cite Berndt2016PHD.
template <int R> class IgnitionDelayKinetics {
public:
  static constexpr int Rank = R;
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);

  IgnitionDelayKinetics(const PerfectGasMix<Rank>& eq);

  Duration ComputeStableDt(const amrex::IntegratorContext&, int level);

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

// We define this class only for dimensions 1 to 3.
// The definitions will be found in its source file IgnitionDelayKinetics.cpp
extern template class IgnitionDelayKinetics<1>;
extern template class IgnitionDelayKinetics<2>;
extern template class IgnitionDelayKinetics<3>;

} // namespace fub::perfect_gas_mix

#endif
