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

#include "fub/ideal_gas/KineticDriver.hpp"
#include "fub/SAMRAI/ScopeGuard.hpp"
#include "fub/geometry/Halfspace.hpp"
#include "fub/geometry/PolymorphicGeometry.hpp"
#include "fub/ideal_gas/FlameMasterKinetics.hpp"
#include "fub/ideal_gas/boundary_condition/IsentropicExpansionBoundary.hpp"
#include "fub/ideal_gas/boundary_condition/PressureValveBoundary.hpp"
#include "fub/ideal_gas/initial_data/Ambient.hpp"
#include "fub/ideal_gas/mechanism/Burke2012.hpp"

#include <cmath>

using namespace fub;

int main(int argc, char** argv) {
  ScopeGuard scope(argc, argv);
  ideal_gas::Burke2012 mechanism;

  // Setup Equation
  const SAMRAI::tbox::Dimension dim(1);
  auto equation = std::make_shared<ideal_gas::FlameMasterKinetics>(
      "IdealGas", dim, mechanism);

  // Setup Domain Sizes
  const CoordinateRange coordinates{Coordinates(1, 0.0), Coordinates(1, 1.5)};
  const int n_cells = 300;

  // Create Driver
  ideal_gas::KineticDriver driver(equation, coordinates, n_cells);

  // Set Boundary Conditions
  const double plenum_pressure = 101325.0; // [Pa]

  ideal_gas::PressureValveOptions valve;
  valve.compressor_pressure = 121325.0; // [Pa]
  valve.ignition_range.lower = Coordinates(1, 0.45);
  valve.ignition_range.upper = Coordinates(1, 0.5);
  valve.flush_air_duration = 2e-3;  // [s]
  valve.flush_fuel_duration = 1e-3; // [s]
  valve.air.temperature = 300;      // [K]
  valve.fuel.temperature = 300;     // [K]

  valve.air.fractions.resize(equation->GetNSpecies());
  ideal_gas::FlameMasterReactor& reactor = equation->GetReactor();
  reactor.setMoleFractions("N2:79,O2:21");
  std::copy(reactor.getMassFractions().begin(),
            reactor.getMassFractions().end(), valve.air.fractions.begin());

  valve.fuel.fractions.resize(equation->GetNSpecies());
  reactor.setMoleFractions("N2:79,O2:21,H2:42");
  std::copy(reactor.getMassFractions().begin(),
            reactor.getMassFractions().end(), valve.fuel.fractions.begin());

  driver.SetLeftBoundaryCondition(
      std::make_shared<ideal_gas::PressureValveBoundary>(
          valve, driver.GetHyperbolicTimeIntegrator(), *equation));

  driver.SetRightBoundaryCondition(
      std::make_shared<ideal_gas::IsentropicExpansionBoundary>(
          plenum_pressure, driver.GetHyperbolicTimeIntegrator(), *equation));

  // Initialize Data
  ideal_gas::Ambient::PrimState state;
  state.temperature = 300.0; // [K]
  state.velocity = 0.0;      // [m/s]
  state.pressure = 101325.0; // [Pa]
  state.species.resize(mechanism.getNSpecies());
  state.species[ideal_gas::Burke2012::sN2] = 79.0; // [-]
  state.species[ideal_gas::Burke2012::sO2] = 21.0; // [-]
  driver.InitializeHierarchy(ideal_gas::Ambient(state, *driver.GetEquation()));

  // Print Grid
  std::vector<double> buffer(equation->GetNVariables() * n_cells);
  using Variable = ideal_gas::IdealGasEquation::Variable;
  auto tock = std::chrono::steady_clock::now();
  // driver.Advance(30e-3, [&](const ideal_gas::KineticDriver& drv, double dt) {
  //   auto tick = std::chrono::steady_clock::now();
  //   std::chrono::duration<double> dur = tick - tock;
  //   if (drv.GetPatchHierarchy()->getMPI().getRank() == 0) {
  //     std::cout << "t = " << drv.TimePoint() << "s, dt = " << dt
  //               << "s, duration = " << dur.count() << "s\n";
  //   }
  //   tock = std::chrono::steady_clock::now();
  // });
  while (driver.TimePoint() < 1.0) {
    driver.Advance(1e-4, [&](const ideal_gas::KineticDriver& drv, double dt) {
      auto tick = std::chrono::steady_clock::now();
      std::chrono::duration<double> dur = tick - tock;
      if (drv.GetPatchHierarchy()->getMPI().getRank() == 0) {
        std::cout << "# t = " << drv.TimePoint() << "s, dt = " << dt
                  << "s, duration = " << dur.count() << "s\n";
      }
      tock = std::chrono::steady_clock::now();
    });
    auto grid = GatherGrid(*driver.GetPatchHierarchy(), *equation, buffer);
    if (grid) {
      const double dx = grid->Dx();
      for (int i = 0; i < grid->Extent(1); ++i) {
        std::cout << ' ' << std::setw(16) << 0.0 + 0.5 * dx + i * dx;
        for (int j = 0; j < grid->Extent(0) - 11; ++j) {
          std::cout << ' ' << std::setw(16) << (*grid)(Variable(j), i);
        }
        constexpr int s0 = static_cast<int>(Variable::species);
        const double rho = (*grid)(Variable::density, i);
        std::cout << ' ' << std::setw(16)
                  << (*grid)(Variable(s0 + ideal_gas::Burke2012::sO2), i) / rho;
        std::cout << ' ' << std::setw(16)
                  << (*grid)(Variable(s0 + ideal_gas::Burke2012::sH2), i) / rho;
        std::cout << '\n';
      }
      std::cout << "\n\n";
    }
  }
}