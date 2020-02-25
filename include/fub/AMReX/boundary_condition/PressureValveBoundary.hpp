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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_PRESSURE_VALVE_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_PRESSURE_VALVE_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/boundary_condition/IsentropicPressureBoundary.hpp"
#include "fub/AMReX/boundary_condition/ReflectiveBoundary.hpp"
#include "fub/AMReX/boundary_condition/MassflowBoundary.hpp"

#include "fub/Duration.hpp"

#include <boost/log/trivial.hpp>
#include <boost/serialization/access.hpp>

#include "fub/ext/ProgramOptions.hpp"

namespace fub::amrex {

struct PressureValveOptions {
  PressureValveOptions() = default;
  PressureValveOptions(const ProgramOptions& vm);

  template <typename Logger> void Print(Logger& log) {
    BOOST_LOG(log) << fmt::format("Pressure Valve '{}' Options:", prefix);
    BOOST_LOG(log) << fmt::format(" - equivalence_ratio = {} [-]",
                                  equivalence_ratio);
    BOOST_LOG(log) << fmt::format(" - outer_pressure = {} [Pa]",
                                  outer_pressure);
    BOOST_LOG(log) << fmt::format(" - outer_temperature = {} [K]",
                                  outer_temperature);
    BOOST_LOG(log) << fmt::format(
        " - pressure_value_which_opens_boundary = {} [Pa]",
        pressure_value_which_opens_boundary);
    BOOST_LOG(log) << fmt::format(
        " - pressure_value_which_closes_boundary = {} [Pa]",
        pressure_value_which_closes_boundary);
    BOOST_LOG(log) << fmt::format(" - oxygen_measurement_position = {} [m]",
                                  oxygen_measurement_position);
    BOOST_LOG(log) << fmt::format(" - oxygen_measurement_criterium = {} [mole]",
                                  oxygen_measurement_criterium);
    BOOST_LOG(log) << fmt::format(" - fuel_measurement_position = {} [m]",
                                  fuel_measurement_position);
    BOOST_LOG(log) << fmt::format(" - fuel_measurement_criterium = {} [-]",
                                  fuel_measurement_criterium);
    BOOST_LOG(log) << fmt::format(" - valve_efficiency = {} [-]",
                                  valve_efficiency);
    BOOST_LOG(log) << fmt::format(" - open_at_interval = {} [s]",
                                  open_at_interval.count());
    BOOST_LOG(log) << fmt::format(" - offset = {} [s]", offset.count());
    BOOST_LOG(log) << fmt::format(" --- Massflow part:")
    massflow_boundary.Print(log);
  }

  std::string prefix{"PressureValve"};
  double equivalence_ratio{1.0};
  double outer_pressure{1.5 * 101325.0};
  double outer_temperature{300.0};
  double pressure_value_which_opens_boundary{101325.0};
  double pressure_value_which_closes_boundary{3.0 * 101325.0};
  double oxygen_measurement_position{1.0};
  double oxygen_measurement_criterium{0.1};
  double fuel_measurement_position{1.0};
  double fuel_measurement_criterium{0.95};
  double valve_efficiency{1.0};
  Duration open_at_interval{0.0};
  Duration offset{0.0};
  MassflowBoundaryOptions massflow_boundary{};
};

} // namespace fub::amrex

namespace boost::serialization {

template <typename Archive>
void serialize(Archive& ar, ::fub::amrex::PressureValveOptions& opts,
               unsigned int version);

}

namespace fub::amrex {

enum class PressureValveState { open_air, open_fuel, closed };

struct PressureValve {
  PressureValveState state{PressureValveState::open_air};
  Duration last_closed{std::numeric_limits<double>::lowest()};
  Duration last_fuel{std::numeric_limits<double>::lowest()};
};

class PressureValveBoundary {
public:
  PressureValveBoundary(const IdealGasMix<1>& equation,
                        PressureValveOptions options);

  PressureValveBoundary(const IdealGasMix<1>& equation,
                        const std::map<std::string, pybind11::object>& options);

  [[nodiscard]] const PressureValveOptions& GetOptions() const noexcept;

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration dt, const GriddingAlgorithm&);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration dt, const GriddingAlgorithm& grid, Direction dir);

  [[nodiscard]] const std::shared_ptr<PressureValve>& GetSharedState() const
      noexcept {
    return shared_valve_;
  }

private:
  PressureValveOptions options_;
  IdealGasMix<1> equation_;
  std::shared_ptr<PressureValve> shared_valve_;
};

} // namespace fub::amrex

namespace boost::serialization {

template <typename Archive>
void serialize(Archive& ar, ::fub::amrex::PressureValve& valve,
               unsigned int /* version */) {
  int state = static_cast<int>(valve.state);
  ar& state;
  valve.state = static_cast<::fub::amrex::PressureValveState>(state);
  ar& valve.last_closed;
  ar& valve.last_fuel;
}

template <typename Archive>
void serialize(Archive& ar, ::fub::amrex::PressureValveOptions& opts,
               unsigned int /* version */) {
  ar& opts.open_at_interval;
  ar& opts.offset;
}

} // namespace boost::serialization

#endif // FINITEVOLUMESOLVER_PRESSUREVALVE_HPP
