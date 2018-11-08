#include "fub/core/assert.hpp"

#include "fub/SAMRAI/ScopeGuard.hpp"
#include "fub/SAMRAI/utility.hpp"
#include "fub/geometry/Halfspace.hpp"
#include "fub/geometry/PolymorphicGeometry.hpp"
#include "fub/ideal_gas/ForwardEulerTimeIntegrator.hpp"
#include "fub/ideal_gas/KineticSourceTerm.hpp"
#include "fub/ideal_gas/PerfectGasEquation.hpp"
#include "fub/ideal_gas/TChemKinetics.hpp"
#include "fub/ideal_gas/boundary_condition/ReflectiveCondition.hpp"
#include "fub/initial_data/RiemannProblem.hpp"
#include "fub/output/GnuplotWriter.hpp"
#include "fub/solver/BoundaryCondition.hpp"
#include "fub/solver/DimensionalSplitSystemSolver.hpp"
#include "fub/solver/GodunovSplitting.hpp"
#include "fub/solver/HyperbolicSystemSourceSolver.hpp"
#include "fub/solver/InitialCondition.hpp"
#include "fub/solver/SourceTermIntegrator.hpp"

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
  using CompleteStatePatchData =
      fub::ideal_gas::IdealGasEquation::CompleteStatePatchData;
  const fub::ideal_gas::IdealGasEquation* ideal_gas;
  PrimState left;
  PrimState right;

  RiemannProblem(fub::PolymorphicGeometry geom,
                 const fub::ideal_gas::IdealGasEquation& equation)
      : fub::RiemannProblem(std::move(geom)), ideal_gas{&equation} {
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
      left.species[3] = 1.0;
      left.species[0] = 1.0;
      right.species[0] = 1.0;
      right.species[3] = 1.0;
    }
  }

private:
  void FillPrimState(const SAMRAI::hier::Patch& patch,
                     const SAMRAI::hier::Index& index,
                     const PrimState& state) const {
    fub::ideal_gas::IdealGasEquation::PrimStatePatchData prim =
        ideal_gas->GetPrimStatePatchData(patch);
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
    ideal_gas->FillFromPrim(ideal_gas->GetCompleteStatePatchData(patch),
                            ideal_gas->GetPrimStatePatchData(patch));
  }
};

int main(int argc, char** argv) {
  if (argc < 3) {
    return -1;
  }
  auto wall_start = std::chrono::steady_clock::now();
  fub::ScopeGuard guard(argc, argv);
  const SAMRAI::tbox::Dimension dim(2);

  fub::ideal_gas::TChemKineticsOptions options;
  options.chemfile = argv[1];
  options.thermofile = argv[2];

  auto equation =
      std::make_shared<fub::ideal_gas::TChemKinetics>("IdealGas", dim, options);

  fub::ideal_gas::ForwardEulerTimeIntegrator integrator(equation);
  fub::ideal_gas::KineticSourceTerm source_term(equation);
  fub::GodunovSplitting splitting{};
  fub::DimensionalSplitSystemSolver hyperbolic(integrator, splitting);
  fub::HyperbolicSystemSourceSolver solver(hyperbolic, source_term, splitting);

  fub::ideal_gas::ReflectiveCondition boundary_condition{integrator};

  fub::Halfspace geometry(fub::Coordinates(dim.getValue(), 1.0), 0.0);
  RiemannProblem initial_data(geometry, *integrator.GetEquation());

  fub::IndexRange indices{SAMRAI::hier::Index(dim, 0),
                          SAMRAI::hier::Index(dim, 200 - 1)};
  fub::CoordinateRange coordinates{fub::Coordinates(dim.getValue(), -1.0), fub::Coordinates(dim.getValue(), +1.0)};

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
        solver.computeStableDt(hierarchy, boundary_condition, t);

    // Do one time step in X (use ghost cells from previous work)
    solver.advanceTime(hierarchy, boundary_condition, t, dt);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::chrono::duration<double> wall_time = end - wall_start;
    if (rank == 0) {
      std::printf("t = %fs, dt = %.3es, wall-time %fs, step-time %fs.\n", t, dt,
                  wall_time.count(), duration.count());
    }
    // Do the output
    t += dt;

    writer.writePlotData(hierarchy, step + 1, t);
  }
}
