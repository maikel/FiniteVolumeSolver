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

#include "fub/solver/euler/IdealGasSplitIntegrator.hpp"

#include "src/solver/euler/hll_flux_method.hpp"
#include "src/solver/euler/hlle_signal_velocities.hpp"

#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/mesh/CascadePartitioner.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/TileClustering.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/OuterfaceData.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"

namespace fub {
namespace euler {
namespace {
std::string getScratchName_(int dir) {
  static constexpr std::string_view scratch = "scratch";
  std::string result{};
  result.reserve(scratch.size() + 3);
  result += scratch;
  result += "_";
  result += std::to_string(dir);
  return result;
}

SAMRAI::hier::VariableDatabase& getDatabase_() {
  return *SAMRAI::hier::VariableDatabase::getDatabase();
}

int reassignContextOfDataId_(
    int index, const std::shared_ptr<SAMRAI::hier::VariableContext> scratch,
    const SAMRAI::hier::IntVector& gcw) {
  SAMRAI::hier::VariableDatabase& database = getDatabase_();
  std::shared_ptr<SAMRAI::hier::Variable> variable;
  database.mapIndexToVariable(index, variable);
  FUB_ASSERT(variable != nullptr);
  return database.registerVariableAndContext(variable, scratch, gcw);
}

std::shared_ptr<SAMRAI::hier::Variable>
getVariableOr_(const std::string& name,
               const std::shared_ptr<SAMRAI::hier::Variable>& alt) {
  std::shared_ptr<SAMRAI::hier::Variable> variable =
      getDatabase_().getVariable(name);
  if (variable) {
    return variable;
  }
  return alt;
}

template <typename VariableType>
std::shared_ptr<SAMRAI::hier::Variable>
getVariable_(const SAMRAI::tbox::Dimension& dim, const std::string& name,
             int depth = 1) {
  std::shared_ptr<SAMRAI::hier::Variable> alt =
      std::make_shared<VariableType>(dim, name, depth);
  return getVariableOr_(name, alt);
}

std::string getPrefixedFluxName_(const std::string& name,
                                 IdealGasSplitIntegrator::FluxVariable var) {
  static constexpr std::array<const char*,
                              IdealGasSplitIntegrator::flux_variables_size>
      names{"/DensityFlux", "/MomentumFlux", "/EnergyFlux"};
  return name + names[static_cast<int>(var)];
}
} // namespace

IdealGasSplitIntegrator::IdealGasSplitIntegrator(const IdealGas& ideal_gas)
    : ideal_gas_{ideal_gas} {
  SAMRAI::hier::VariableDatabase& database = getDatabase_();
  const SAMRAI::tbox::Dimension& dim = ideal_gas_.getDim();
  const int dim_value = dim.getValue();
  for (int dir = 0; dir < dim_value; ++dir) {
    std::shared_ptr<SAMRAI::hier::VariableContext> scratch =
        database.getContext(getScratchName_(dir));
    SAMRAI::hier::IntVector ghost_width = SAMRAI::hier::IntVector::getZero(dim);
    ghost_width[dir] = 1;
    for (int v = 0; v < variables_size; ++v) {
      scratch_[v + dir * variables_size] =
          reassignContextOfDataId_(current(Variable(v)), scratch, ghost_width);
    }
  }
  const std::string& name = ideal_gas_.name();
  const SAMRAI::hier::IntVector zeros = SAMRAI::hier::IntVector::getZero(dim);
  std::shared_ptr<SAMRAI::hier::VariableContext> current =
      database.getContext("current");
  for (int v = 0; v < flux_variables_size; ++v) {
    int depth =
        (FluxVariable(v) == FluxVariable::momentum) ? dim.getValue() : 1;
    auto flux = getVariable_<SAMRAI::pdat::FaceVariable<double>>(
        dim, getPrefixedFluxName_(name, FluxVariable(v)), depth);
    face_[v] = database.registerVariableAndContext(flux, current, zeros);
  }
}

namespace {
void selectVariablesToAllocate_(const IdealGasSplitIntegrator& integrator,
                                SAMRAI::hier::ComponentSelector& variables) {
  for (int id : integrator.current()) {
    variables.setFlag(id);
  }
  for (int id : integrator.scratch()) {
    variables.setFlag(id);
  }
  for (int id : integrator.face()) {
    variables.setFlag(id);
  }
}
} // namespace

void IdealGasSplitIntegrator::initializePatchHierarchy(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    double time_point) const {
  struct InitializeStrategy : public SAMRAI::mesh::StandardTagAndInitStrategy {
    const IdealGasSplitIntegrator* that;

    explicit InitializeStrategy(const IdealGasSplitIntegrator& solv)
        : that{&solv} {}

    void initializeLevelData(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
        int level_number, double init_data_time, bool can_be_refined,
        bool initial_time,
        const std::shared_ptr<SAMRAI::hier::PatchLevel>& old_level,
        bool allocate_data) override {
      const std::shared_ptr<SAMRAI::hier::PatchLevel>& patch_level =
          hierarchy->getPatchLevel(level_number);
      if (allocate_data) {
        SAMRAI::hier::ComponentSelector variables_to_allocate;
        selectVariablesToAllocate_(*that, variables_to_allocate);
        patch_level->allocatePatchData(variables_to_allocate, initial_time);
      }
      if (that->initial_condition_) {
        for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : *patch_level) {
          that->initial_condition_->initializeDataOnPatch(*patch);
        }
      }
    }

    void resetHierarchyConfiguration(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
        int coarsest_level, int finest_level) override {}
  };
  InitializeStrategy initialize(*this);
  const SAMRAI::tbox::Dimension& dim = hierarchy->getDim();
  SAMRAI::mesh::GriddingAlgorithm algorithm(
      hierarchy, "", nullptr,
      std::make_shared<SAMRAI::mesh::StandardTagAndInitialize>("", &initialize),
      std::make_shared<SAMRAI::mesh::TileClustering>(dim),
      std::make_shared<SAMRAI::mesh::CascadePartitioner>(dim, ""));
  algorithm.makeCoarsestLevel(time_point);
}

namespace {
struct BoundaryConditionPatchStrategy_ : SAMRAI::xfer::RefinePatchStrategy {
  const BoundaryCondition* boundary_condition_;

  explicit BoundaryConditionPatchStrategy_(const BoundaryCondition* bc)
      : boundary_condition_{bc} {}

  void setPhysicalBoundaryConditions(
      SAMRAI::hier::Patch& patch, double fill_time,
      const SAMRAI::hier::IntVector& ghost_width_to_fill) override {
    if (boundary_condition_) {
      boundary_condition_->setPhysicalBoundaryConditions(patch, fill_time,
                                                         ghost_width_to_fill);
    }
  }

  SAMRAI::hier::IntVector
  getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const override {
    if (boundary_condition_) {
      return boundary_condition_->getStencilWidth(dim);
    }
    return SAMRAI::hier::IntVector::getZero(dim);
  }

  void preprocessRefine(SAMRAI::hier::Patch& fine,
                        const SAMRAI::hier::Patch& coarse,
                        const SAMRAI::hier::Box& fine_box,
                        const SAMRAI::hier::IntVector& ratio) override {}

  void postprocessRefine(SAMRAI::hier::Patch& fine,
                         const SAMRAI::hier::Patch& coarse,
                         const SAMRAI::hier::Box& fine_box,
                         const SAMRAI::hier::IntVector& ratio) override {}
};
} // namespace

void IdealGasSplitIntegrator::doFillGhostLayer(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    double time_point, int dir) const {
  SAMRAI::xfer::RefineAlgorithm algorithm{};
  using Var = Variable;
  algorithm.registerRefine(scratch(Var::density, dir), current(Var::density),
                           scratch(Var::density, dir), nullptr);
  algorithm.registerRefine(scratch(Var::momentum, dir), current(Var::momentum),
                           scratch(Var::momentum, dir), nullptr);
  algorithm.registerRefine(scratch(Var::energy, dir), current(Var::energy),
                           scratch(Var::energy, dir), nullptr);
  algorithm.registerRefine(scratch(Var::pressure, dir), current(Var::pressure),
                           scratch(Var::pressure, dir), nullptr);
  algorithm.registerRefine(scratch(Var::speed_of_sound, dir),
                           current(Var::speed_of_sound),
                           scratch(Var::speed_of_sound, dir), nullptr);
  BoundaryConditionPatchStrategy_ boundary_condition{boundary_condition_.get()};
  const int max_level_number = hierarchy->getNumberOfLevels();
  for (int level_number = 0; level_number < max_level_number; ++level_number) {
    algorithm
        .createSchedule(hierarchy->getPatchLevel(level_number),
                        &boundary_condition)
        ->fillData(time_point);
  }
}

double IdealGasSplitIntegrator::estimateStableTimeStepSize(
    const SAMRAI::hier::Patch& patch, double time_point) const {
  using Var = IdealGas::Variable;
  SAMRAI::pdat::CellData<double>& density =
      getScratchData(patch, Var::density, 0);
  SAMRAI::pdat::CellData<double>& momentum =
      getScratchData(patch, Var::momentum, 0);
  SAMRAI::pdat::CellData<double>& speed_of_sound =
      getScratchData(patch, Var::speed_of_sound, 0);
  double velocity = 0.0;
  const SAMRAI::hier::Box& patch_box = patch.getBox();
  auto first = SAMRAI::pdat::FaceGeometry::begin(patch_box, 0);
  auto last = SAMRAI::pdat::FaceGeometry::end(patch_box, 0);
  while (first != last) {
    using FaceIndex = SAMRAI::pdat::FaceIndex;
    FaceIndex face = *first++;
    SAMRAI::pdat::CellIndex left(face.toCell(FaceIndex::Lower));
    SAMRAI::pdat::CellIndex right(face.toCell(FaceIndex::Upper));
    const double rhoL = density(left);
    const double rhoR = density(right);
    const double rhouL = momentum(left, 0);
    const double rhouR = momentum(right, 0);
    const double aL = speed_of_sound(left);
    const double aR = speed_of_sound(right);

    auto [sL, sR] =
        computeHlleSignalVelocities({rhoL, rhouL, aL}, {rhoR, rhouR, aR});
    velocity = std::max({velocity, std::abs(sL), std::abs(sR)});
  }
  const double* dx = getCartesianPatchGeometry(patch)->getDx();
  return 0.5 * dx[0] / velocity;
}

namespace {
void computeFluxes_(const IdealGasSplitIntegrator& integrator,
                    const SAMRAI::hier::Patch& patch, int dir) {
  using Var = IdealGas::Variable;
  SAMRAI::pdat::CellData<double>& density =
      integrator.getScratchData(patch, Var::density, dir);
  SAMRAI::pdat::CellData<double>& momentum =
      integrator.getScratchData(patch, Var::momentum, dir);
  SAMRAI::pdat::CellData<double>& energy =
      integrator.getScratchData(patch, Var::energy, dir);
  SAMRAI::pdat::CellData<double>& pressure =
      integrator.getScratchData(patch, Var::pressure, dir);
  SAMRAI::pdat::CellData<double>& speed_of_sound =
      integrator.getScratchData(patch, Var::speed_of_sound, dir);

  using Flux = IdealGasSplitIntegrator::FluxVariable;
  SAMRAI::pdat::FaceData<double>& density_flux =
      integrator.getFaceData(patch, Flux::density);
  SAMRAI::pdat::FaceData<double>& momentum_flux =
      integrator.getFaceData(patch, Flux::momentum);
  SAMRAI::pdat::FaceData<double>& energy_flux =
      integrator.getFaceData(patch, Flux::energy);

  const SAMRAI::hier::Box& patch_box = patch.getBox();
  auto first = SAMRAI::pdat::FaceGeometry::begin(patch_box, dir);
  auto last = SAMRAI::pdat::FaceGeometry::end(patch_box, dir);
  while (first != last) {
    using FaceIndex = SAMRAI::pdat::FaceIndex;
    FaceIndex face = *first++;
    SAMRAI::pdat::CellIndex left(face.toCell(FaceIndex::Lower));
    SAMRAI::pdat::CellIndex right(face.toCell(FaceIndex::Upper));

    const double rhoL = density(left);
    const double rhoR = density(right);
    const double rhouL = momentum(left, dir);
    const double rhouR = momentum(right, dir);
    const double aL = speed_of_sound(left);
    const double aR = speed_of_sound(right);

    HllSignals<double> signals =
        computeHlleSignalVelocities({rhoL, rhouL, aL}, {rhoR, rhouR, aR});

    const double rhoeL = energy(left);
    const double rhoeR = energy(right);
    const double pL = pressure(left);
    const double pR = pressure(right);

    auto [fRho, fRhoU, fRhoE] = computeHllFlux(
        signals, {rhoL, rhouL, rhoeL, pL}, {rhoR, rhouR, rhoeR, pR});

    density_flux(face) = fRho;
    const int dim = integrator.ideal_gas().getDim().getValue();
    for (int d = 0; d < dim; ++d) {
      const double sL = signals.left;
      const double sR = signals.right;
      const double fRhoUL = rhouL / rhoL * momentum(left, d);
      const double fRhoUR = rhouR / rhoR * momentum(right, d);
      momentum_flux(face, d) =
          (d == dir) ? fRhoU
                     : computeHllFlux(sL, sR, momentum(left, d),
                                      momentum(right, d), fRhoUL, fRhoUR);
    }
    energy_flux(face) = fRhoE;
  }
}

void advanceConservatively_(const IdealGasSplitIntegrator& integrator,
                            const SAMRAI::hier::Patch& patch, double dt,
                            int dir) {
  using Var = IdealGas::Variable;
  SAMRAI::pdat::CellData<double>& density =
      integrator.getScratchData(patch, Var::density, dir);
  SAMRAI::pdat::CellData<double>& momentum =
      integrator.getScratchData(patch, Var::momentum, dir);
  SAMRAI::pdat::CellData<double>& energy =
      integrator.getScratchData(patch, Var::energy, dir);

  using Flux = IdealGasSplitIntegrator::FluxVariable;
  SAMRAI::pdat::FaceData<double>& density_flux =
      integrator.getFaceData(patch, Flux::density);
  SAMRAI::pdat::FaceData<double>& momentum_flux =
      integrator.getFaceData(patch, Flux::momentum);
  SAMRAI::pdat::FaceData<double>& energy_flux =
      integrator.getFaceData(patch, Flux::energy);

  const double* dx = getCartesianPatchGeometry(patch)->getDx();
  const double lambda = dt / dx[dir];
  const SAMRAI::hier::Box& patch_box = patch.getBox();
  for (const SAMRAI::hier::Index& index : patch_box) {
    SAMRAI::pdat::CellIndex cell(index);
    using SAMRAI::pdat::FaceIndex;
    FaceIndex left(index, dir, FaceIndex::Lower);
    FaceIndex right(index, dir, FaceIndex::Upper);

    density(cell) += lambda * (density_flux(left) - density_flux(right));
    const int dim = integrator.ideal_gas().getDim().getValue();
    for (int d = 0; d < dim; ++d) {
      momentum(cell, d) +=
          lambda * (momentum_flux(left, d) - momentum_flux(right, d));
    }
    energy(cell) += lambda * (energy_flux(left) - energy_flux(right));
  }
}

void reconstructStates_(const IdealGasSplitIntegrator& integrator,
                        const SAMRAI::hier::Patch& patch, int dir) {
  using Var = IdealGas::Variable;
  SAMRAI::pdat::CellData<double>& density =
      integrator.getScratchData(patch, Var::density, dir);
  SAMRAI::pdat::CellData<double>& momentum =
      integrator.getScratchData(patch, Var::momentum, dir);
  SAMRAI::pdat::CellData<double>& energy =
      integrator.getScratchData(patch, Var::energy, dir);
  SAMRAI::pdat::CellData<double>& pressure =
      integrator.getScratchData(patch, Var::pressure, dir);
  SAMRAI::pdat::CellData<double>& speed_of_sound =
      integrator.getScratchData(patch, Var::speed_of_sound, dir);

  const SAMRAI::hier::Box& patch_box = patch.getBox();
  for (const SAMRAI::hier::Index& index : patch_box) {
    SAMRAI::pdat::CellIndex cell(index);
    IdealGas::Conservative<double> cons;
    cons.density = density(cell);
    cons.momentum = momentum(cell, dir);
    cons.energy = energy(cell);
    const double p = integrator.ideal_gas().computePressure(cons);
    const double a =
        integrator.ideal_gas().computeSpeedOfSound(cons.density, p);
    pressure(cell) = p;
    speed_of_sound(cell) = a;
  }
}
} // namespace

void IdealGasSplitIntegrator::integratePatch(const SAMRAI::hier::Patch& patch,
                                             double time_point,
                                             double time_step_size,
                                             int dir) const {
  computeFluxes_(*this, patch, dir);
  advanceConservatively_(*this, patch, time_step_size, dir);
}

void IdealGasSplitIntegrator::coarsenFluxOnCoarseFineInteraces(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy) const {}

void IdealGasSplitIntegrator::fixPatchConservationOnCoarseFineBoundary(
    const SAMRAI::hier::Patch& patch, const SAMRAI::hier::BoundaryBox& boundary,
    double time_step_size) const {}

void IdealGasSplitIntegrator::postIntegratePatchHierarchy(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    double time_point, double time_step_size, int dir) const {
  forEachPatch(*hierarchy, [&](const SAMRAI::hier::Patch& patch) {
    reconstructStates_(*this, patch, dir);
  });
  SAMRAI::xfer::RefineAlgorithm algorithm{};
  using Var = Variable;
  algorithm.registerRefine(current(Var::density), scratch(Var::density, dir),
                           current(Var::density), nullptr);
  algorithm.registerRefine(current(Var::momentum), scratch(Var::momentum, dir),
                           current(Var::momentum), nullptr);
  algorithm.registerRefine(current(Var::energy), scratch(Var::energy, dir),
                           current(Var::energy), nullptr);
  algorithm.registerRefine(current(Var::pressure), scratch(Var::pressure, dir),
                           current(Var::pressure), nullptr);
  algorithm.registerRefine(current(Var::speed_of_sound),
                           scratch(Var::speed_of_sound, dir),
                           current(Var::speed_of_sound), nullptr);
  const int max_level_number = hierarchy->getNumberOfLevels();
  for (int level_number = 0; level_number < max_level_number; ++level_number) {
    algorithm.createSchedule(hierarchy->getPatchLevel(level_number))
        ->fillData(time_point + time_step_size);
  }
}

} // namespace euler
} // namespace fub