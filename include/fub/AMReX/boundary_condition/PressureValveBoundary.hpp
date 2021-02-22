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
#include "fub/AMReX/boundary_condition/MassflowBoundary.hpp"
#include "fub/AMReX/boundary_condition/ReflectiveBoundary.hpp"

#include "fub/Duration.hpp"

#include <boost/serialization/access.hpp>

#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include <limits>
#include <string>

namespace fub::amrex {
/// \ingroup BoundaryCondition
///
struct PressureValveOptions {
  PressureValveOptions() = default;
  PressureValveOptions(const ProgramOptions& vm);

  void Print(SeverityLogger& log);

  std::string channel{"PressureValve"};
  double equivalence_ratio{1.0};
  double pressure_value_which_closes_boundary{3.0 * 101325.0};
  Duration change_to_fuel_at_interval{0.0};
  Duration change_to_fuel_time_offset{0.0};
  MassflowBoundaryOptions massflow_boundary{};
};

enum class PressureValveState { open_air, open_fuel, closed };
/// \ingroup BoundaryCondition
///
struct PressureValve {
  PressureValveState state{PressureValveState::open_air};
  Duration last_fuel{std::numeric_limits<double>::lowest()};
};
}

namespace boost::serialization {
template <typename Archive>
void serialize(Archive& ar, ::fub::amrex::PressureValve& valve,
               unsigned int /* version */) {
  // clang-format off
  int state = static_cast<int>(valve.state);
  ar & state;
  valve.state = static_cast<::fub::amrex::PressureValveState>(state);
  ar & valve.last_fuel;
  // clang-format on
}
}

namespace fub::amrex {
/// \ingroup BoundaryCondition
///
class PressureValveBoundary {
public:
  PressureValveBoundary(const IdealGasMix<1>& equation,
                        PressureValveOptions options);

  [[nodiscard]] const PressureValveOptions& GetOptions() const noexcept;

  [[nodiscard]] const PressureValve& GetValve() const noexcept;

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);


private:
  template <typename Archive>
  friend void serialize(Archive& ar, PressureValveBoundary& boundary,
                        unsigned int /* version */) {
    // clang-format off
    ar & boundary.valve_;
    // clang-format on
  }

  PressureValveOptions options_;
  IdealGasMix<1> equation_;
  PressureValve valve_{};
};
} // namespace fub::amrex


#endif // FINITEVOLUMESOLVER_PRESSUREVALVE_HPP
