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

#include "fub/ideal_gas/boundary_condition/PressureValveBoundary.hpp"

#include "fub/core/algorithm.hpp"

#include <numeric>

namespace fub {
namespace ideal_gas {

namespace {
void SetTemperature_(SAMRAI::hier::Patch& patch, std::ptrdiff_t i,
                     double temperature, FlameMasterKinetics& equation,
                     span<double> buffer) {
  IdealGasEquation::CompleteStatePatchData q =
      equation.GetCompleteStatePatchData(patch);
  SAMRAI::pdat::CellIndex cell(
      SAMRAI::hier::Index(SAMRAI::tbox::Dimension(1), i));
  CopyMassFractions(buffer, q.species, cell);
  FlameMasterReactor& reactor = equation.GetReactor();
  reactor.setMassFractions(buffer);
  reactor.setTemperature(temperature);
  reactor.setPressure(q.pressure(cell));
  const double u = q.momentum(cell) / q.density(cell);
  equation.UpdateCellFromReactor(q, cell, reactor, u);
}

void IgniteDetonation_(SAMRAI::hier::PatchHierarchy& hier,
                       const CoordinateRange& ignition_range,
                       FlameMasterKinetics& equation) {
  std::vector<double> buffer(equation.GetNSpecies());
  forEachPatch(hier, [&](SAMRAI::hier::Patch& patch) {
    const double x_lo = *GetCartesianPatchGeometry(patch)->getXLower();
    const double x_up = *GetCartesianPatchGeometry(patch)->getXUpper();
    const double ir_lo = ignition_range.lower[0];
    const double ir_up = ignition_range.upper[0]; 
    const double dx = *GetCartesianPatchGeometry(patch)->getDx();
    const std::ptrdiff_t i_lo = patch.getBox().lower(0);
    const std::ptrdiff_t i_up = patch.getBox().upper(0);
    const std::ptrdiff_t first = fub::clamp<std::ptrdiff_t>(
        i_lo + std::floor((ir_lo - x_lo) / dx), i_lo, i_up);
    const std::ptrdiff_t last = fub::clamp<std::ptrdiff_t>(
        i_lo + std::floor((ir_up - x_lo) / dx), i_lo, i_up);
    for (std::ptrdiff_t i = first; i < last; ++i) {
      SetTemperature_(patch, i, 2000, equation, buffer);
    }
  });
}

bool MaybeChangeState_(PressureValveBoundary& boundary,
                       SAMRAI::hier::PatchHierarchy& hier, double timepoint) {
  const double current_state_duration =
      timepoint - boundary.LastChangeTimepoint();
  PressureValveState state = boundary.GetValveState();
  if (state == PressureValveState::air) {
    const double duration_until_change = boundary.Options().flush_air_duration;
    if (duration_until_change < current_state_duration) {
      boundary.ChangeStateTo(PressureValveState::fuel, timepoint);
      return true;
    }
  } else {
    FUB_ASSERT(state == PressureValveState::fuel);
    const double duration_until_change = boundary.Options().flush_fuel_duration;
    if (duration_until_change < current_state_duration) {
      IgniteDetonation_(hier, boundary.Options().ignition_range,
                        boundary.GetEquation());
      boundary.ChangeStateTo(PressureValveState::air, timepoint);
      return true;
    }
  }
  return false;
}

std::vector<double> MakeState_(const PressureValveOptions& options,
                               PressureValveState state,
                               FlameMasterKinetics& equation) {
  using Variable = IdealGasEquation::Variable;
  std::vector<double> buffer(equation.GetNVariables());
  FlameMasterReactor& reactor = equation.GetReactor();
  reactor.setTemperature(state == PressureValveState::air
                             ? options.air.temperature
                             : options.fuel.temperature);
  reactor.setMassFractions(state == PressureValveState::air
                               ? options.air.fractions
                               : options.fuel.fractions);
  reactor.setPressure(options.compressor_pressure);
  equation.UpdateStateFromReactor(buffer, reactor, 45.0);
  return buffer;
}
} // namespace

PressureValveBoundary::PressureValveBoundary(
    const PressureValveOptions& options, const HyperbolicTimeIntegrator& i,
    FlameMasterKinetics& kinetics)
    : options_(options), equation_{&kinetics},
      pressure_boundary_(
          MakeState_(options, PressureValveState::fuel, kinetics), i),
      reflective_boundary_(i) {}

span<const double>
PressureValveBoundary::GetMassFractions(PressureValveState state) const {
  if (state == PressureValveState::air) {
    return options_.air.fractions;
  }
  return options_.fuel.fractions;
}

void PressureValveBoundary::ChangeStateTo(PressureValveState state,
                                          double timepoint) {
  if (valve_state_ != state) {
    valve_state_ = state;
    state_changed_timepoint_ = timepoint;
    pressure_boundary_.SetStates(MakeState_(options_, state, GetEquation()));
  }
}

void PressureValveBoundary::SetPhysicalBoundaryCondition(
    const SAMRAI::hier::Patch& patch, const SAMRAI::hier::Box& fill_box,
    double fill_time, Direction dir, int side) const {
  if (!is_opened_) {
    reflective_boundary_.SetPhysicalBoundaryCondition(patch, fill_box,
                                                      fill_time, dir, side);
  } else {
    pressure_boundary_.SetPhysicalBoundaryCondition(patch, fill_box, fill_time,
                                                    dir, side);
  }
}

void PressureValveBoundary::PreFillGhostLayer(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hier,
    double timepoint) {
  MaybeChangeState_(*this, *hier, timepoint);
  struct MeanPressure {
    double pressure;
    double volume;
  };
  MeanPressure local = ReducePatches(
      *hier, MeanPressure{},
      [&](const MeanPressure& value, const SAMRAI::hier::Patch& patch) {
        const double x_lower = *GetCartesianPatchGeometry(patch)->getXLower();
        const double x_upper = *GetCartesianPatchGeometry(patch)->getXUpper();
        const double volume = x_upper - x_lower;
        const SAMRAI::pdat::CellData<double>& data = GetEquation().GetCellData(
            patch, IdealGasEquation::Variable::pressure);
        const std::ptrdiff_t size = data.getBox().size();
        const double pressure =
            std::accumulate(data.getPointer(), data.getPointer() + size, 0.0) /
            size;
        return MeanPressure{value.pressure + volume * pressure,
                            value.volume + volume};
      });
  MPI_Comm comm = hier->getMPI().getCommunicator();
  MeanPressure global;
  MPI_Allreduce(&local.pressure, &global.pressure, 2, MPI_DOUBLE, MPI_SUM,
                comm);
  observed_mean_pressure_ = global.pressure / global.volume;
  if (observed_mean_pressure_ > options_.close_boundary_pressure_limit) {
    is_opened_ = false;
  } else if (observed_mean_pressure_ < options_.open_boundary_pressure_limit) {
    is_opened_ = true;
  }
}

} // namespace ideal_gas
} // namespace fub
