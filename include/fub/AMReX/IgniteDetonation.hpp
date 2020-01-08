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

#ifndef FUB_AMREX_IGNITE_DETONATION_HPP
#define FUB_AMREX_IGNITE_DETONATION_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/equations/IdealGasMix.hpp"
#include "fub/ext/outcome.hpp"

#include "fub/ext/ProgramOptions.hpp"

#include <boost/log/trivial.hpp>
#include <boost/serialization/access.hpp>

namespace fub::amrex {

struct IgniteDetonationOptions {
  IgniteDetonationOptions() = default;

  explicit IgniteDetonationOptions(
      const std::map<std::string, pybind11::object>& vm,
      const std::string& prefix = "ignite");

  template <typename Logger> void Print(Logger& log) const {
    // clang-format off
    BOOST_LOG(log) << fmt::format("Ignite Detonation '{}' Options:", prefix);
    BOOST_LOG(log) << fmt::format("  - measurement_position = {} [m]", measurement_position);
    BOOST_LOG(log) << fmt::format("  - equivalence_ratio_criterium = {} [-]", equivalence_ratio_criterium);
    BOOST_LOG(log) << fmt::format("  - temperature_low = {} [K]", temperature_low);
    BOOST_LOG(log) << fmt::format("  - temperature_high = {} [K]", temperature_high);
    BOOST_LOG(log) << fmt::format("  - ramp_width = {} [m]", ramp_width);
    BOOST_LOG(log) << fmt::format("  - ignite_position = {} [m]", ignite_position);
    BOOST_LOG(log) << fmt::format("  - ignite_interval = {} [s]", ignite_interval.count());
    // clang-format on
  }

  std::string prefix{"ignite"};
  double measurement_position{1.0};
  double equivalence_ratio_criterium{0.95};
  double temperature_low{300.0};
  double temperature_high{2000.0};
  double ramp_width{0.05};
  double ignite_position{0.0};
  Duration ignite_interval{0.0};
};

} // namespace fub::amrex

namespace boost::serialization {

template <typename Archive>
void serialize(Archive& ar, ::fub::amrex::IgniteDetonationOptions& opts,
               unsigned int version);

} // namespace boost::serialization

namespace fub::amrex {

class IgniteDetonation {
public:
  static constexpr int Rank = 1;

  /// \brief Constructs the source term
  ///
  /// \param[in] eq  The equation which knows the species
  /// 
  /// \param[in] max_refinement_level  This is needed to allocate storate for
  /// ignition times at each refinement level.
  ///
  /// \param[in] opts  Options which
  /// manipulate the way of igniting the the gas mixture
  IgniteDetonation(const IdealGasMix<1>& eq, int max_refinement_level,
                   const IgniteDetonationOptions& opts = {});

  /// \brief Resets internal configuration
  ///
  /// Resets the last ignition time points if the time point on grid is lower
  /// than the last recorded ignition. This happens when the CFL condition gets
  /// violated right after the ignition.
  void ResetHierarchyConfiguration(std::shared_ptr<GriddingAlgorithm> grid);

  /// \brief Returns numeric_limits<double>::max()
  ///
  /// This operator artificially ignites a detonation and has no restriction on
  /// the time step size.
  [[nodiscard]] Duration ComputeStableDt(int) const noexcept;

  /// Uses the scratch space of simulation_data to evaluate the criterion on the
  /// current equivalence ratio.
  [[nodiscard]] Result<void, TimeStepTooLarge>
  AdvanceLevel(IntegratorContext& simulation_data, int level, Duration dt,
               const ::amrex::IntVect& ngrow = ::amrex::IntVect(0));

  /// \brief Returns the options for this operator
  [[nodiscard]] const IgniteDetonationOptions& GetOptions() const noexcept {
    return options_;
  }

  /// \brief Returns the time points for the last ignition on refinement level
  /// `level`.
  [[nodiscard]] Duration GetLastIgnitionTimePoint(int level) const noexcept;

  /// \brief Set a time point for an ignition on refinement level `level`.
  void SetLastIgnitionTimePoint(int level, Duration t) noexcept;

private:
  IdealGasMix<1> equation_;
  IgniteDetonationOptions options_;
  std::vector<Duration> last_ignition_backup_{};
  std::vector<Duration> last_ignition_{};

  friend class boost::serialization::access;
  template <typename Archive> void serialize(Archive& ar, unsigned int version);
};

template <typename Archive>
void IgniteDetonation::serialize(Archive& ar, unsigned int /* version */) {
  // clang-format off
  ar & options_;
  ar & last_ignition_;
  // clang-format on
}

} // namespace fub::amrex

namespace boost::serialization {

template <typename Archive>
void serialize(Archive& ar, ::fub::amrex::IgniteDetonationOptions& opts,
               unsigned int /* version */) {
  // clang-format off
  ar & opts.equivalence_ratio_criterium;
  ar & opts.measurement_position;
  ar & opts.ramp_width;
  ar & opts.temperature_high;
  ar & opts.temperature_low;
  ar & opts.ignite_interval;
  ar & opts.ignite_position;
  // clang-format on
}

} // namespace boost::serialization

#endif
