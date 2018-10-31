#include "fub/core/assert.hpp"

#include "fub/SAMRAI/ScopeGuard.hpp"
#include "fub/SAMRAI/utility.hpp"
#include "fub/solver/BoundaryCondition.hpp"
#include "fub/solver/GodunovSplitting.hpp"
#include "fub/solver/HyperbolicSystemSolver.hpp"
#include "fub/solver/InitialCondition.hpp"
#include "fub/solver/euler/ForwardEulerTimeIntegrator.hpp"

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

struct CircleData : fub::InitialCondition {
  const fub::euler::ForwardEulerTimeIntegrator* integrator;

  explicit CircleData(const fub::euler::ForwardEulerTimeIntegrator& i)
      : integrator{&i} {}

  struct State {
    double rho;
    double p;
    double E;
    double a;
  };

  std::array<State, 2> getStates(fub::FlameMasterReactor& reactor) const {
    std::vector<double> Y(reactor.GetNSpecies());
    Y[0] = 1.0;
    reactor.SetMassFractions(Y);
    reactor.SetTemperature(2000.);
    reactor.SetPressure(8 * 101325.0);
    double gamma = reactor.GetCp() / reactor.GetCv();

    State inner;
    inner.rho = reactor.GetDensity();
    inner.p = reactor.GetPressure();
    inner.E = reactor.GetInternalEnergy();
    inner.a = std::sqrt(gamma * inner.p / inner.rho);

    reactor.SetMassFractions(Y);
    reactor.SetTemperature(300.);
    reactor.SetPressure(1 * 101325.0);
    gamma = reactor.GetCp() / reactor.GetCv();

    State outer;
    outer.rho = reactor.GetDensity();
    outer.p = reactor.GetPressure();
    outer.E = reactor.GetInternalEnergy();
    outer.a = std::sqrt(gamma * outer.p / outer.rho);

    return {inner, outer};
  }

  static double Dot_(const fub::Coordinates& x, const fub::Coordinates& y) {
    const int size = x.size();
    double total = 0.0;
    for (int i = 0; i < size; ++i) {
      total += x[i] * y[i];
    }
    return total;
  }

  void initializeDataOnPatch(const SAMRAI::hier::Patch& patch) const override {
    auto state = integrator->getCompleteState(patch);
    double* density = state.density.getPointer();
    double* momentum = state.momentum.getPointer();
    double* energy = state.energy.getPointer();
    double* pressure = state.pressure.getPointer();
    double* speed_of_sound = state.speed_of_sound.getPointer();
    const auto& geom = *fub::getCartesianPatchGeometry(patch);
    const SAMRAI::hier::Box& box = patch.getBox();
    const double sqrt_ = std::sqrt(1.4 * 8.0);
    fub::FlameMasterReactor reactor = integrator->ideal_gas()->GetReactor();
#ifdef __cpp_structured_bindings
    const auto [inner, outer] = getStates(reactor);
#else
    const std::array<State, 2> states = getStates(reactor);
    const State& inner = states[0];
    const State& outer = states[1];
#endif
    for (SAMRAI::hier::Index i : box) {
      fub::Coordinates x = fub::computeCellCoordinates(geom, box, i);
      SAMRAI::pdat::CellIndex cell(i);
      const double radius2 = Dot_(x, x);
      if (radius2 < 0.025) {
        state.density(cell) = inner.rho;
        state.pressure(cell) = inner.p;
        state.energy(cell) = inner.E;
        state.speed_of_sound(cell) = inner.a;
        state.temperature(cell) = 300.0;
      } else {
        state.density(cell) = outer.rho;
        state.pressure(cell) = outer.p;
        state.energy(cell) = outer.E;
        state.speed_of_sound(cell) = outer.a;
        state.temperature(cell) = 300.0;
      }
    }
  }
};

struct ConstantBoundary : public fub::BoundaryCondition {
  const fub::euler::ForwardEulerTimeIntegrator* integrator;

  explicit ConstantBoundary(const fub::euler::ForwardEulerTimeIntegrator& i)
      : integrator{&i} {}

  struct State {
    double rho;
    double p;
    double E;
    double a;
  };

  void setPhysicalBoundaryConditions(
      const SAMRAI::hier::Patch& patch, double fill_time,
      const SAMRAI::hier::IntVector& ghost_width_to_fill) const override {
    const SAMRAI::hier::PatchGeometry& patch_geom = *patch.getPatchGeometry();
    const std::vector<SAMRAI::hier::BoundaryBox>& faces =
        patch_geom.getCodimensionBoundaries(1);
    using Variable = fub::euler::IdealGas::Variable;
    int dim = 0;
    for (int d = 0; d < ghost_width_to_fill.getDim().getValue(); ++d) {
      if (ghost_width_to_fill[d]) {
        dim = d;
        break;
      }
    }
    using Scratch = fub::euler::ForwardEulerTimeIntegrator::Scratch;
    auto state = integrator->getCompleteState(patch, Scratch(dim));
    SAMRAI::pdat::CellData<double>& density = state.density;
    SAMRAI::pdat::CellData<double>& momentum = state.momentum;
    SAMRAI::pdat::CellData<double>& energy = state.energy;
    SAMRAI::pdat::CellData<double>& pressure = state.pressure;
    SAMRAI::pdat::CellData<double>& temperature = state.temperature;
    SAMRAI::pdat::CellData<double>& speed_of_sound = state.speed_of_sound;

    fub::FlameMasterReactor reactor = integrator->ideal_gas()->GetReactor();
    std::vector<double> Y(reactor.GetNSpecies());
    Y[0] = 1.0;
    reactor.SetMassFractions(Y);
    reactor.SetTemperature(300.);
    reactor.SetPressure(101325.0);
    const double gamma = reactor.GetCp() / reactor.GetCv();

    State outer;
    outer.rho = reactor.GetDensity();
    outer.p = reactor.GetPressure();
    outer.E = reactor.GetInternalEnergy();
    outer.a = std::sqrt(gamma * outer.p / outer.rho);
    for (const SAMRAI::hier::BoundaryBox& face : faces) {
      SAMRAI::hier::Box box = patch_geom.getBoundaryFillBox(
          face, patch.getBox(), ghost_width_to_fill);
      for (SAMRAI::hier::Index i : box) {
        SAMRAI::pdat::CellIndex cell(i);
        density(cell) = outer.rho;
        momentum(cell) = 0.0;
        pressure(cell) = outer.p;
        energy(cell) = outer.E;
        speed_of_sound(cell) = outer.a;
        temperature(cell) = 300.0;
      }
    }
  }

  SAMRAI::hier::IntVector
  getStencilWidth(const SAMRAI::tbox::Dimension& dim) const override {
    return SAMRAI::hier::IntVector::getZero(dim);
  }
};

int main(int argc, char** argv) {
  fub::ScopeGuard guard(argc, argv);
  const SAMRAI::tbox::Dimension dim(2);

  fub::euler::ForwardEulerTimeIntegrator integrator(
      std::make_shared<fub::euler::IdealGas>("IdealGas", dim));
  ConstantBoundary boundary_cond(integrator);
  CircleData initial_data(integrator);

  fub::IndexRange indices{SAMRAI::hier::Index(dim, 0),
                          SAMRAI::hier::Index(dim, 200 - 1)};
  fub::CoordinateRange coordinates{{-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}};

  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy =
      makeCartesianPatchHierarchy(indices, coordinates);

  fub::initializePatchHierarchy(hierarchy, integrator, initial_data);

  std::unique_ptr<SAMRAI::appu::VisItDataWriter> writer = nullptr;
  if (dim.getValue() == 2 || dim.getValue() == 3) {
    writer = std::make_unique<SAMRAI::appu::VisItDataWriter>(dim, "VisItWriter",
                                                             "output");
    using Var = fub::euler::IdealGas::Variable;
    writer->registerPlotQuantity("IdealGas_Density", "SCALAR",
                                 integrator.getPatchDataId(Var::density));
    writer->registerPlotQuantity("IdealGas_Momentum", "VECTOR",
                                 integrator.getPatchDataId(Var::momentum));
    writer->registerPlotQuantity("IdealGas_Pressure", "SCALAR",
                                 integrator.getPatchDataId(Var::pressure));
    writer->registerPlotQuantity("IdealGas_Temperature", "SCALAR",
                                 integrator.getPatchDataId(Var::temperature));
    writer->registerPlotQuantity(
        "IdealGas_SpeedOfSound", "SCALAR",
        integrator.getPatchDataId(Var::speed_of_sound));
    writer->registerPlotQuantity("IdealGas_Species", "VECTOR",
                                 integrator.getPatchDataId(Var::species));
    // Do the output
    writer->writePlotData(hierarchy, 0, 0.0);
  }

  double t = 0;
  const int rank = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank();
  fub::GodunovSplitting splitting{};
  fub::HyperbolicSystemSolver solver(integrator, splitting);
  for (int step = 0; step < 10; ++step) {
    auto start = std::chrono::steady_clock::now();

    // Estimate time step size
    const double dt = solver.ComputeStableDt(hierarchy, boundary_cond, t);
    // Do one time step in X (use ghost cells from previous work)
    solver.AdvanceTime(hierarchy, boundary_cond, t, dt);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;
    if (rank == 0) {
      std::printf("Time Step Counter: %fs.\n", duration.count());
    }
    // Do the output
    t += dt;
    if (writer) {
      writer->writePlotData(hierarchy, step + 1, t);
    }
  }
}
