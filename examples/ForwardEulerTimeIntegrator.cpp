#include "fub/core/assert.hpp"

#include "fub/solver/euler/ForwardEulerTimeIntegrator.hpp"

#include "fub/SAMRAI/ScopeGuard.hpp"
#include "fub/SAMRAI/utility.hpp"
#include "fub/solver/BoundaryCondition.hpp"
#include "fub/solver/InitialCondition.hpp"

#include <array>

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"

#include "SAMRAI/appu/VisItDataWriter.h"

struct CircleData : fub::InitialCondition {
  const fub::euler::ForwardEulerTimeIntegrator* integrator;

  explicit CircleData(const fub::euler::ForwardEulerTimeIntegrator& i)
      : integrator{&i} {}

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
    for (SAMRAI::hier::Index i : box) {
      fub::Coordinates x = fub::computeCellCoordinates(geom, box, i);
      const double radius2 = x[0] * x[0] + x[1] * x[1];
      *density++ = (radius2 <= 0.025) * 1.0 + (radius2 > 0.025) * 0.125;
      *momentum++ = 0.0;
      *pressure++ = (radius2 <= 0.025) * 8.0 + (radius2 > 0.025) * 1.0;
      *energy++ = 20.0;
      *speed_of_sound++ = sqrt_;
    }
  }
};

struct ConstantBoundary : public fub::BoundaryCondition {
  const fub::euler::ForwardEulerTimeIntegrator* integrator;

  explicit ConstantBoundary(const fub::euler::ForwardEulerTimeIntegrator& i)
      : integrator{&i} {}

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
    SAMRAI::pdat::CellData<double>& speed_of_sound = state.speed_of_sound;
    for (const SAMRAI::hier::BoundaryBox& face : faces) {
      SAMRAI::hier::Box box = patch_geom.getBoundaryFillBox(
          face, patch.getBox(), ghost_width_to_fill);
      for (SAMRAI::hier::Index i : box) {
        SAMRAI::pdat::CellIndex cell(i);
        density(cell) = 0.125;
        momentum(cell) = 0.0;
        pressure(cell) = 1.0;
        energy(cell) = 20.0;
        speed_of_sound(cell) = std::sqrt(1.4 * 8.0);
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
      fub::euler::IdealGas("IdealGas", dim));

  fub::IndexRange indices{SAMRAI::hier::Index(0, 0),
                          SAMRAI::hier::Index(127, 127)};

  fub::CoordinateRange coordinates{{-1.0, -1.0}, {1.0, 1.0}};

  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy =
      makeCartesianPatchHierarchy(indices, coordinates);

  CircleData initial_data{integrator};
  fub::initializePatchHierarchy(hierarchy, integrator, initial_data);

  SAMRAI::appu::VisItDataWriter writer(dim, "VisItWriter", "output");
  using Var = fub::euler::IdealGas::Variable;
  writer.registerPlotQuantity("IdealGas_Density", "SCALAR",
                              integrator.getPatchDataId(Var::density));
  writer.registerPlotQuantity("IdealGas_Momentum", "VECTOR",
                              integrator.getPatchDataId(Var::momentum));
  writer.registerPlotQuantity("IdealGas_Pressure", "SCALAR",
                              integrator.getPatchDataId(Var::pressure));

  // Do the output
  writer.writePlotData(hierarchy, 0, 0.0);

  ConstantBoundary boundary_cond{integrator};

  // Estimate time step size
  integrator.fillGhostLayer(hierarchy, boundary_cond, 0.0, fub::Direction::X);
  const double dt = integrator.computeStableDt(*hierarchy, 0.0);

  // Do one time step in X (use ghost cells from previous work)
  integrator.advanceTime(hierarchy, 0.0, dt, fub::Direction::X);

  // Do one time step in Y (transfer ghost cell layer for this direction)
  integrator.fillGhostLayer(hierarchy, boundary_cond, 0.0, fub::Direction::Y);
  integrator.advanceTime(hierarchy, 0.0, dt, fub::Direction::Y);

  // Do the output
  writer.writePlotData(hierarchy, 1, dt);
}