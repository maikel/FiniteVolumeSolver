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

#include "fub/core/assert.hpp"

#include "fub/SAMRAI/ScopeGuard.hpp"
#include "fub/SAMRAI/utility.hpp"
#include "fub/geometry/Halfspace.hpp"
#include "fub/geometry/PolymorphicGeometry.hpp"
#include "fub/SAMRAI/ideal_gas/FlameMasterKinetics.hpp"
#include "fub/SAMRAI/ideal_gas/DimensionalSplitTimeIntegrator.hpp"
#include "fub/SAMRAI/ideal_gas/KineticSourceTerm.hpp"
#include "fub/SAMRAI/ideal_gas/PerfectGasEquation.hpp"
#include "fub/SAMRAI/ideal_gas/boundary_condition/ReflectiveBoundary.hpp"
#include "fub/ideal_gas/mechanism/Zhao2008Dme.hpp"
#include "fub/SAMRAI/initial_data/RiemannProblem.hpp"
#include "fub/SAMRAI/output/GnuplotWriter.hpp"
#include "fub/SAMRAI/BoundaryCondition.hpp"
#include "fub/SAMRAI/DimensionalSplitSystemSolver.hpp"
#include "fub/SAMRAI/GodunovSplitting.hpp"
#include "fub/SAMRAI/DimensionalSplitSystemSourceSolver.hpp"
#include "fub/SAMRAI/InitialCondition.hpp"
#include "fub/SAMRAI/SourceTermIntegrator.hpp"

#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"

#include <array>
#include <chrono>
#include <cstdio>

struct PrimState {
  double temperature;
  double momentum;
  double pressure;
  std::vector<double> species;
};

struct CompleteState {
  double density;
  double momentum;
  double energy;
  double pressure;
  double temperature;
  double speed_of_sound;
  std::vector<double> species;
};

fub::ideal_gas::IdealGasEquation::PrimStateSpan<double> MakeSpan(PrimState& w) {
  return {{&w.temperature, 1}, {&w.momentum, 1}, {&w.pressure, 1}, w.species};
}

fub::ideal_gas::IdealGasEquation::PrimStateSpan<const double>
MakeSpan(const PrimState& w) {
  return {{&w.temperature, 1}, {&w.momentum, 1}, {&w.pressure, 1}, w.species};
}

class RiemannProblem : public fub::RiemannProblem {
public:
  using CompletePatchData =
      fub::ideal_gas::IdealGasEquation::CompletePatchData;
  const fub::ideal_gas::IdealGasEquation* ideal_gas;
  PrimState left;
  PrimState right;

  RiemannProblem(fub::PolymorphicGeometry geom,
                 const fub::ideal_gas::IdealGasEquation& equation)
      : fub::RiemannProblem(std::move(geom)), ideal_gas{&equation} {
    using fub::ideal_gas::Zhao2008Dme;
    left.temperature = 1100.0;
    left.momentum = 0.0;
    left.pressure = 101325.0;
    right.temperature = 300.0;
    right.momentum = 0.0;
    right.pressure = 101325.0;
    const int n_species = equation.GetNSpecies();
    if (n_species) {
      left.species.resize(n_species);
      right.species.resize(n_species);
      left.species[Zhao2008Dme::sH2] = 0.5;
      left.species[Zhao2008Dme::sO2] = 0.5;
      right.species[Zhao2008Dme::sO2] = 0.5;
      right.species[Zhao2008Dme::sH2] = 0.5;
    }
  }

private:
  void FillPrimState(const SAMRAI::hier::Patch& patch,
                     const SAMRAI::hier::Index& index,
                     const PrimState& state) const {
    fub::ideal_gas::IdealGasEquation::PrimPatchData prim =
        ideal_gas->GetPrimPatchData(patch);
    SAMRAI::pdat::CellIndex i(index);
    prim.momentum(i) = state.momentum;
    prim.pressure(i) = state.pressure;
    prim.temperature(i) = state.temperature;
    for (int s = 0; s < state.species.size(); ++s) {
      prim.species(i, s) = state.species[s];
    }
  }

  void FillLeftState(const SAMRAI::hier::Patch& patch,
                     const SAMRAI::hier::Index& index) const override {
    FillPrimState(patch, index, left);
  }

  void FillRightState(const SAMRAI::hier::Patch& patch,
                      const SAMRAI::hier::Index& index) const override {
    FillPrimState(patch, index, right);
  }

  void PostInitialize(const SAMRAI::hier::Patch& patch) const override {
    ideal_gas->FillFromPrim(ideal_gas->GetCompletePatchData(patch),
                            ideal_gas->GetPrimPatchData(patch));
  }
};

int main(int argc, char** argv) {
  auto wall_start = std::chrono::steady_clock::now();
  fub::ScopeGuard guard(argc, argv);
  const SAMRAI::tbox::Dimension dim(1);

  fub::ideal_gas::Zhao2008Dme mechanism;
  auto equation = std::make_shared<fub::ideal_gas::FlameMasterKinetics>(
      "IdealGas", dim, mechanism);

  fub::ideal_gas::DimensionalSplitTimeIntegrator integrator(equation);
  fub::ideal_gas::KineticSourceTerm source_term(equation);
  fub::GodunovSplitting splitting{};
  fub::DimensionalSplitSystemSolver hyperbolic(integrator, splitting);
  fub::DimensionalSplitSystemSourceSolver solver(hyperbolic, source_term, splitting);

  fub::ideal_gas::ReflectiveBoundary boundary_condition{integrator};

  fub::Halfspace geometry(fub::Coordinates(1, 1.0), 0.0);
  RiemannProblem initial_data(geometry, *integrator.GetEquation());

  fub::IndexRange indices{SAMRAI::hier::Index(dim, 0),
                          SAMRAI::hier::Index(dim, 200 - 1)};
  fub::CoordinateRange coordinates{fub::Coordinates(1, -1.0),
                                   fub::Coordinates(1, +1.0)};

  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy =
      fub::MakeCartesianPatchHierarchy(indices, coordinates);

  fub::InitializePatchHierarchy(hierarchy, integrator, initial_data);

  fub::GnuplotWriter writer("output", equation->GetPatchDataIds());
  // Do the output
  writer.writePlotData(hierarchy, 0, 0.0);

  double t = 0;
  const int rank = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank();
  for (int step = 0; step < 200; ++step) {
    auto start = std::chrono::steady_clock::now();

    // Estimate time step size
    const double dt =
        0.5 * solver.ComputeStableDt(hierarchy, boundary_condition, t);

    // Do one time step in X (use ghost cells from previous work)
    solver.AdvanceTime(hierarchy, boundary_condition, t, dt);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::chrono::duration<double> wall_time = end - wall_start;
    if (rank == 0) {
      std::printf("t = %fs, dt = %.3es, wall-time %fs, step-time %fs.\n", t, dt,
                  wall_time.count(), duration.count());
    }
    // Do the output
    t += dt;

    // writer.writePlotData(hierarchy, step + 1, t);
  }
}
