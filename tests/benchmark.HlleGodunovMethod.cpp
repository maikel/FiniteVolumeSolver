#include <benchmark/benchmark.h>

#include "fub/SAMRAI/ideal_gas/HlleGodunovMethod.hpp"
#include "fub/ideal_gas/flux_method/HlleGodunovMethod.hpp"

#include "fub/SAMRAI/ScopeGuard.hpp"
#include "fub/SAMRAI/utility.hpp"
#include "fub/geometry/Halfspace.hpp"
#include "fub/geometry/PolymorphicGeometry.hpp"
#include "fub/ideal_gas/mechanism/Zhao2008Dme.hpp"

#include "fub/SAMRAI/ideal_gas/initial_data/RiemannProblem.hpp"

#include "fub/SAMRAI/BoundaryCondition.hpp"
#include "fub/SAMRAI/DimensionalSplitSystemSolver.hpp"
#include "fub/SAMRAI/DimensionalSplitSystemSourceSolver.hpp"
#include "fub/SAMRAI/GodunovSplitting.hpp"
#include "fub/SAMRAI/InitialCondition.hpp"
#include "fub/SAMRAI/SourceTermIntegrator.hpp"
#include "fub/SAMRAI/ideal_gas/DimensionalSplitTimeIntegrator.hpp"
#include "fub/SAMRAI/ideal_gas/FlameMasterKinetics.hpp"
#include "fub/SAMRAI/ideal_gas/KineticSourceTerm.hpp"
#include "fub/SAMRAI/ideal_gas/PerfectGasEquation.hpp"
#include "fub/SAMRAI/ideal_gas/boundary_condition/ReflectiveBoundary.hpp"
#include "fub/SAMRAI/output/GnuplotWriter.hpp"

static fub::ideal_gas::Zhao2008Dme mechanism;
static std::shared_ptr<fub::ideal_gas::FlameMasterKinetics> equation{};
static std::unique_ptr<fub::ideal_gas::DimensionalSplitTimeIntegrator>
    integrator{};
static std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy{};

static void ComputeFluxes(benchmark::State& state) {
  fub::samrai::ideal_gas::HlleGodunovMethod method;
  for (auto&& _ : state) {
    fub::ForEachPatch(*hierarchy, [&](SAMRAI::hier::Patch& patch) {
      using Scratch = fub::ideal_gas::DimensionalSplitTimeIntegrator::Scratch;
      auto fluxes = integrator->GetFluxPatchData(patch);
      auto states = integrator->GetCompletePatchData(patch, Scratch::X);
      method.ComputeFluxesOnPatch(fluxes, states, patch, 0.5,
                                  fub::Direction::X);
    });
  }
}
BENCHMARK(ComputeFluxes);

static void ComputeFluxes_Eigen(benchmark::State& state) {
  fub::ideal_gas::HlleGodunovMethod method;
  for (auto&& _ : state) {
    fub::ForEachPatch(*hierarchy, [&](SAMRAI::hier::Patch& patch) {
      using Scratch = fub::ideal_gas::DimensionalSplitTimeIntegrator::Scratch;
      using Direction = fub::Direction;
      auto fluxes = integrator->GetFluxSpans(patch, Direction::X);
      auto states = integrator->GetCompleteSpans(patch, Scratch::X);
      method.ComputeNumericFluxesOnSpans(fluxes, states);
    });
  }
}
BENCHMARK(ComputeFluxes_Eigen);

int main(int argc, char** argv) {
  fub::ScopeGuard guard(argc, argv);
  const SAMRAI::tbox::Dimension dim(3);
  equation = std::make_shared<fub::ideal_gas::FlameMasterKinetics>(
      "IdealGas", SAMRAI::tbox::Dimension(3), mechanism);
  integrator = std::make_unique<fub::ideal_gas::DimensionalSplitTimeIntegrator>(
      equation);
  fub::ideal_gas::RiemannProblem::PrimState left;
  left.temperature = 345.0;     // [K]
  left.pressure = 8 * 101325.0; // [Pa]
  left.velocity = 0.0;          // [m/s]
  left.species = std::vector<double>(mechanism.getNSpecies());
  left.species[0] = 1.0;

  fub::ideal_gas::RiemannProblem::PrimState right;
  right.temperature = 300.0; // [K]
  right.pressure = 101325.0; // [Pa]
  right.velocity = 0.0;      // [m/s]
  right.species = std::vector<double>(mechanism.getNSpecies());
  right.species[0] = 1.0;

  fub::Halfspace geometry(fub::Coordinates(dim.getValue(), 1.0), 0.0);

  fub::ideal_gas::RiemannProblem initial_data(left, right, geometry,
                                              *integrator->GetEquation());

  fub::IndexRange indices{SAMRAI::hier::Index(dim, 0),
                          SAMRAI::hier::Index(dim, 50 - 1)};

  fub::CoordinateRange coordinates{fub::Coordinates(dim.getValue(), -1.0),
                                   fub::Coordinates(dim.getValue(), +1.0)};

  hierarchy = fub::MakeCartesianPatchHierarchy(indices, coordinates);

  fub::InitializePatchHierarchy(hierarchy, *integrator, initial_data);
  fub::ideal_gas::ReflectiveBoundary boundary_condition{*integrator};
  integrator->FillGhostLayer(hierarchy, boundary_condition, 0.0,
                             fub::Direction::X);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}