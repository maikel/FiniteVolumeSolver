// Copyright (c) 2018 Maikel Nadolski
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

#ifndef FUB_EULER_BOUNDARY_CONDITION_PRESSURE_VALVE_BOUNDARY_HPP
#define FUB_EULER_BOUNDARY_CONDITION_PRESSURE_VALVE_BOUNDARY_HPP

#include "fub/core/mdspan.hpp"
#include "fub/ideal_gas/FlameMasterKinetics.hpp"
#include "fub/ideal_gas/HyperbolicTimeIntegrator.hpp"
#include "fub/ideal_gas/boundary_condition/ConstantBoundary.hpp"
#include "fub/ideal_gas/boundary_condition/ReflectiveBoundary.hpp"
#include "fub/solver/SplitBoundaryCondition.hpp"

namespace fub {
namespace ideal_gas {

struct PressureValveOptions {
  double compressor_pressure;
  double open_boundary_pressure_limit;
  double close_boundary_pressure_limit;

  CoordinateRange ignition_range;
  double flush_air_duration;
  double flush_fuel_duration;

  struct State {
    double temperature;
    std::vector<double> fractions;
  };
  State air;
  State fuel;
};

enum class PressureValveState { air, fuel };

class PressureValveBoundary : public SplitBoundaryCondition {
public:
  // Constructor

  PressureValveBoundary(const PressureValveOptions& options,
                        const HyperbolicTimeIntegrator& i,
                        FlameMasterKinetics& kinetics);

  // Virtual Methods

  void SetPhysicalBoundaryCondition(const SAMRAI::hier::Patch& patch,
                                    const SAMRAI::hier::Box& fill_box,
                                    double fill_time, Direction dir,
                                    int side) const override;

  int GetStencilWidth() const override { return 0; }

  void
  PreFillGhostLayer(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& grid,
                    double time) override;

  // Observers

  const PressureValveOptions& Options() const noexcept { return options_; }

  double LastChangeTimepoint() const noexcept {
    return state_changed_timepoint_;
  }

  FlameMasterKinetics& GetEquation() noexcept {
    return *equation_;
  }

  const FlameMasterKinetics& GetEquation() const noexcept {
    return *equation_;
  }

  PressureValveState GetValveState() const noexcept { return valve_state_; }

  span<const double> GetMassFractions(PressureValveState state) const;

  // Modifiers

  void ChangeStateTo(PressureValveState state, double timepoint);

private:
  PressureValveOptions options_;
  FlameMasterKinetics* equation_;
  ConstantBoundary pressure_boundary_;
  ReflectiveBoundary reflective_boundary_;
  PressureValveState valve_state_{PressureValveState::fuel};
  double state_changed_timepoint_{0.0};
  double observed_mean_pressure_{options_.compressor_pressure};
  bool is_opened_{true};
};

} // namespace ideal_gas
} // namespace fub

#endif