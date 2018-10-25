#include "fub/core/assert.hpp"

#include "fub/solver/euler/IdealGasSplitIntegrator.hpp"
// #include "fub/solver/hyperbolic_system_solver.hpp"

#include "fub/SAMRAI/ScopeGuard.hpp"
#include "fub/SAMRAI/utility.hpp"

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
  const fub::euler::IdealGasSplitIntegrator* integrator;

  explicit CircleData(const fub::euler::IdealGasSplitIntegrator* i)
      : integrator{i} {}

  void initializeDataOnPatch(const SAMRAI::hier::Patch& patch) const override {
    using Variable = fub::euler::IdealGas::Variable;
    const auto& geom = *fub::getCartesianPatchGeometry(patch);
    auto getPointer = [&](Variable var) {
      return integrator->getCurrentData(patch, var).getPointer();
    };
    double* density = getPointer(Variable::density);
    double* momentum = getPointer(Variable::momentum);
    double* energy = getPointer(Variable::energy);
    double* pressure = getPointer(Variable::pressure);
    double* speed_of_sound = getPointer(Variable::speed_of_sound);
    const SAMRAI::hier::Box& box = patch.getBox();
    for (SAMRAI::hier::Index i : box) {
      fub::Coordinates x = fub::computeCellCoordinates(geom, box, i);
      const double radius2 = x[0] * x[0] + x[1] * x[1];
      if (radius2 < 0.025) {
        *density++ = 1.0;
        *momentum++ = 0.0;
        *pressure++ = 8.0;
        *energy++ = 20.0;
        *speed_of_sound++ = std::sqrt(1.4 * 8.0);
      } else {
        *density++ = 0.125;
        *momentum++ = 0.0;
        *pressure++ = 1.0;
        *energy++ = 20.0;
        *speed_of_sound++ = std::sqrt(1.4 * 8.0);
      }
    }
  }
};

struct ConstantBoundary : public fub::BoundaryCondition {
  const fub::euler::IdealGasSplitIntegrator* integrator;

  explicit ConstantBoundary(const fub::euler::IdealGasSplitIntegrator* integ)
      : integrator{integ} {}

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
    auto getData = [&](Variable var) -> SAMRAI::pdat::CellData<double>& {
      return integrator->getScratchData(patch, var, dim);
    };
    SAMRAI::pdat::CellData<double>& density = getData(Variable::density);
    SAMRAI::pdat::CellData<double>& momentum = getData(Variable::momentum);
    SAMRAI::pdat::CellData<double>& energy = getData(Variable::energy);
    SAMRAI::pdat::CellData<double>& pressure = getData(Variable::pressure);
    SAMRAI::pdat::CellData<double>& speed_of_sound = getData(Variable::speed_of_sound);
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
  fub::SAMRAI::ScopeGuard guard(argc, argv);
  const SAMRAI::tbox::Dimension dim(2);

  fub::euler::IdealGasSplitIntegrator integrator(
      fub::euler::IdealGas("IdealGas", dim));

  integrator.setInitialCondition(std::make_shared<CircleData>(&integrator));
  integrator.setBoundaryCondition(
      std::make_shared<ConstantBoundary>(&integrator));

  fub::IndexRange indices{SAMRAI::hier::Index(0, 0),
                          SAMRAI::hier::Index(127, 127)};

  fub::CoordinateRange coordinates{{-1.0, -1.0}, {1.0, 1.0}};

  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy =
      integrator.makePatchHierarchy({indices, coordinates, 2});

  SAMRAI::appu::VisItDataWriter writer(dim, "VisItWriter", "output");
  using Var = fub::euler::IdealGas::Variable;
  writer.registerPlotQuantity("Density", "SCALAR",
                              integrator.current(Var::density));
  writer.registerPlotQuantity("Pressure", "SCALAR",
                              integrator.current(Var::pressure));

  // Do the output
  writer.writePlotData(hierarchy, 0, 0.0);

  // Estimate time step size
  integrator.fillGhostLayer(hierarchy, 0.0, 0);
  const double dt = integrator.estimateHierarchyTimeStepSize(*hierarchy, 0.0);

  // Do one time step in X (use ghost cells from previous work)
  integrator.integratePatchHierarchy(hierarchy, 0.0, dt, 0);

  // Do one time step in Y (transfer ghost cell layer for this direction)
  integrator.fillGhostLayer(hierarchy, 0.0, 1);
  integrator.integratePatchHierarchy(hierarchy, 0.0, dt, 1);

  // Do the output
  writer.writePlotData(hierarchy, 1, dt);
}