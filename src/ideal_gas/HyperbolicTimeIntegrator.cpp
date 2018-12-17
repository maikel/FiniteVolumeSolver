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

#include "fub/ideal_gas/HyperbolicTimeIntegrator.hpp"
#include "fub/ideal_gas/HlleRiemannSolver.hpp"

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
namespace ideal_gas {
namespace {
#ifdef __cpp_lib_string_view
std::string getScratchName_(int dir) {
  static constexpr std::string_view scratch = "scratch";
  std::string result{};
  result.reserve(scratch.size() + 3);
  result += scratch;
  result += "_";
  result += std::to_string(dir);
  return result;
}
#else
std::string getScratchName_(int dir) {
  static constexpr const char* scratch = "scratch";
  std::string result{};
  result.reserve(std::strlen(scratch) + 3);
  result += scratch;
  result += "_";
  result += std::to_string(dir);
  return result;
}
#endif

int reassignContextOfDataId_(
    int index, const std::shared_ptr<SAMRAI::hier::VariableContext> scratch,
    const SAMRAI::hier::IntVector& gcw) {
  SAMRAI::hier::VariableDatabase& database =
      *SAMRAI::hier::VariableDatabase::getDatabase();
  std::shared_ptr<SAMRAI::hier::Variable> variable;
  database.mapIndexToVariable(index, variable);
  FUB_ASSERT(variable != nullptr);
  return database.registerVariableAndContext(variable, scratch, gcw);
}

std::shared_ptr<SAMRAI::hier::Variable>
getVariableOr_(const std::string& name,
               const std::shared_ptr<SAMRAI::hier::Variable>& alt) {
  std::shared_ptr<SAMRAI::hier::Variable> variable =
      SAMRAI::hier::VariableDatabase::getDatabase()->getVariable(name);
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
                                 HyperbolicTimeIntegrator::FluxVariable var) {
  static constexpr std::array<const char*,
                              HyperbolicTimeIntegrator::flux_variables_size>
      names{"/DensityFlux", "/MomentumFlux", "/EnergyFlux", "/SpeciesFlux"};
  return name + names[static_cast<int>(var)];
}

int getDepth_(HyperbolicTimeIntegrator::FluxVariable var, int dim,
              int specs) {
  switch (var) {
  case HyperbolicTimeIntegrator::FluxVariable::momentum:
    return dim;
  case HyperbolicTimeIntegrator::FluxVariable::species:
    return specs;
  default:
    return 1;
  }
}
} // namespace

HyperbolicTimeIntegrator::HyperbolicTimeIntegrator(
    std::shared_ptr<const IdealGasEquation> ideal_gas)
    : ideal_gas_{std::move(ideal_gas)},
      flux_method_{std::make_shared<HlleRiemannSolver>(ideal_gas_)} {
  SAMRAI::hier::VariableDatabase& database =
      *SAMRAI::hier::VariableDatabase::getDatabase();
  const SAMRAI::tbox::Dimension& dim = ideal_gas_->GetDimension();
  const int dim_value = dim.getValue();
  const SAMRAI::hier::IntVector ghost_cell_width =
      flux_method_->getStencilWidth(dim);
  for (int dir = 0; dir < dim_value; ++dir) {
    std::shared_ptr<SAMRAI::hier::VariableContext> scratch =
        database.getContext(getScratchName_(dir));
    SAMRAI::hier::IntVector ghost = SAMRAI::hier::IntVector::getZero(dim);
    ghost[dir] = ghost_cell_width[dir];
    for (int v = 0; v < variables_size; ++v) {
      scratch_[v + dir * variables_size] =
          reassignContextOfDataId_(getPatchDataId(Variable(v)), scratch, ghost);
    }
  }
  const std::string& name = ideal_gas_->GetName();
  const SAMRAI::hier::IntVector zeros = SAMRAI::hier::IntVector::getZero(dim);
  std::shared_ptr<SAMRAI::hier::VariableContext> current =
      database.getContext("current");
  for (int v = 0; v < flux_variables_size; ++v) {
    const int n_species = ideal_gas_->GetNSpecies();
    const int depth = getDepth_(FluxVariable(v), dim_value, n_species);
    auto flux = getVariable_<SAMRAI::pdat::FaceVariable<double>>(
        dim, getPrefixedFluxName_(name, FluxVariable(v)), depth);
    face_[v] = database.registerVariableAndContext(flux, current, zeros);
  }
}

double HyperbolicTimeIntegrator::computeStableDtOnPatch(
    const SAMRAI::hier::Patch& patch, double time_point, Direction dir) const {
  CompleteStatePatchData state = getCompleteState(patch, Scratch(dir));
  return flux_method_->computeStableDtOnPatch(state, patch, dir);
}

namespace {
void AdvanceTimeOnPatch_(
    const HyperbolicTimeIntegrator::CompleteStatePatchData& state,
    const HyperbolicTimeIntegrator::FluxStatePatchData& flux,
    const SAMRAI::hier::Patch& patch, double dt, Direction dir) {
  const double* dx = GetCartesianPatchGeometry(patch)->getDx();
  const double lambda = dt / dx[static_cast<int>(dir)];
  const SAMRAI::hier::Box& patch_box = patch.getBox();
  for (const SAMRAI::hier::Index& index : patch_box) {
    SAMRAI::pdat::CellIndex cell(index);
    using SAMRAI::pdat::FaceIndex;
    FaceIndex left(index, static_cast<int>(dir), FaceIndex::Lower);
    FaceIndex right(index, static_cast<int>(dir), FaceIndex::Upper);

    state.density(cell) += lambda * (flux.density(left) - flux.density(right));
    FUB_ASSERT(state.density(cell) > 0);
    const int dim = patch.getDim().getValue();
    for (int d = 0; d < dim; ++d) {
      state.momentum(cell, d) +=
          lambda * (flux.momentum(left, d) - flux.momentum(right, d));
    }
    state.energy(cell) += lambda * (flux.energy(left) - flux.energy(right));
    const int species_size = state.species.getDepth();
    for (int s = 0; s < species_size; ++s) {
      state.species(cell, s) +=
          lambda * (flux.species(left, s) - flux.species(right, s));
      // FUB_ASSERT(state.species(cell, s) >= 0);
    }
  }
}
} // namespace

void HyperbolicTimeIntegrator::advanceTimeOnPatch(
    const SAMRAI::hier::Patch& patch, double time_point, double dt,
    Direction dir) const {
  FluxStatePatchData flux = getFluxState(patch);
  CompleteStatePatchData state = getCompleteState(patch, Scratch(dir));
  flux_method_->computeFluxesOnPatch(flux, state, patch, dt, dir);
  AdvanceTimeOnPatch_(state, flux, patch, dt, dir);
  CompleteStatePatchData complete = getCompleteState(patch);
  ConsStatePatchData cons = getConsState(patch, Scratch(dir));
  ideal_gas_->FillFromCons(complete, cons);
}

namespace {
SAMRAI::pdat::CellData<double>& getCellData_(const SAMRAI::hier::Patch& patch,
                                             int data_id) {
  return *static_cast<SAMRAI::pdat::CellData<double>*>(
      patch.getPatchData(data_id).get());
}

SAMRAI::pdat::FaceData<double>& getFaceData_(const SAMRAI::hier::Patch& patch,
                                             int data_id) {
  return *static_cast<SAMRAI::pdat::FaceData<double>*>(
      patch.getPatchData(data_id).get());
}
} // namespace

HyperbolicTimeIntegrator::CompleteStatePatchData
HyperbolicTimeIntegrator::getCompleteState(const SAMRAI::hier::Patch& patch,
                                             Scratch scratch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, getPatchDataId(Variable::density, scratch));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, getPatchDataId(Variable::momentum, scratch));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, getPatchDataId(Variable::energy, scratch));
  SAMRAI::pdat::CellData<double>& pressure = getCellData_(patch, getPatchDataId(Variable::pressure, scratch));
  SAMRAI::pdat::CellData<double>& temperature = getCellData_(patch, getPatchDataId(Variable::temperature, scratch));
  SAMRAI::pdat::CellData<double>& speed_of_sound = getCellData_(patch, getPatchDataId(Variable::speed_of_sound, scratch));
  SAMRAI::pdat::CellData<double>& species = getCellData_(patch, getPatchDataId(Variable::species, scratch));
  // clang-format on
  return {density,     momentum,       energy, pressure,
          temperature, speed_of_sound, species};
}

HyperbolicTimeIntegrator::CompleteStatePatchData
HyperbolicTimeIntegrator::getCompleteState(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, getPatchDataId(Variable::density));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, getPatchDataId(Variable::momentum));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, getPatchDataId(Variable::energy));
  SAMRAI::pdat::CellData<double>& pressure = getCellData_(patch, getPatchDataId(Variable::pressure));
  SAMRAI::pdat::CellData<double>& temperature = getCellData_(patch, getPatchDataId(Variable::temperature));
  SAMRAI::pdat::CellData<double>& speed_of_sound = getCellData_(patch, getPatchDataId(Variable::speed_of_sound));
  SAMRAI::pdat::CellData<double>& species = getCellData_(patch, getPatchDataId(Variable::species));
  // clang-format on
  return {density,     momentum,       energy, pressure,
          temperature, speed_of_sound, species};
}

HyperbolicTimeIntegrator::ConsStatePatchData
HyperbolicTimeIntegrator::getConsState(const SAMRAI::hier::Patch& patch,
                                         Scratch scratch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, getPatchDataId(Variable::density, scratch));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, getPatchDataId(Variable::momentum, scratch));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, getPatchDataId(Variable::energy, scratch));
  SAMRAI::pdat::CellData<double>& species = getCellData_(patch, getPatchDataId(Variable::species, scratch));
  // clang-format on
  return {density, momentum, energy, species};
}

HyperbolicTimeIntegrator::ConsStatePatchData
HyperbolicTimeIntegrator::getConsState(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, getPatchDataId(Variable::density));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, getPatchDataId(Variable::momentum));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, getPatchDataId(Variable::energy));
  SAMRAI::pdat::CellData<double>& species = getCellData_(patch, getPatchDataId(Variable::species));
  // clang-format on
  return {density, momentum, energy, species};
}

HyperbolicTimeIntegrator::FluxStatePatchData
HyperbolicTimeIntegrator::getFluxState(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::FaceData<double>& density = getFaceData_(patch, getPatchDataId(FluxVariable::density));
  SAMRAI::pdat::FaceData<double>& momentum = getFaceData_(patch, getPatchDataId(FluxVariable::momentum));
  SAMRAI::pdat::FaceData<double>& energy = getFaceData_(patch, getPatchDataId(FluxVariable::energy));
  SAMRAI::pdat::FaceData<double>& species = getFaceData_(patch, getPatchDataId(FluxVariable::species));
  // clang-format on
  return {density, momentum, energy, species};
}

void HyperbolicTimeIntegrator::flagPatchDataIdsToAllocate(
    SAMRAI::hier::ComponentSelector& which_to_allocate) const {
  for (int id : ideal_gas_->GetPatchDataIds()) {
    which_to_allocate.setFlag(id);
  }
  for (int id : getScratch()) {
    which_to_allocate.setFlag(id);
  }
  for (int id : face_) {
    which_to_allocate.setFlag(id);
  }
}

std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>
HyperbolicTimeIntegrator::getFillGhostLayerRefineAlgorithm(
    Direction dir) const {
  auto algorithm = std::make_shared<SAMRAI::xfer::RefineAlgorithm>();
  for (int v = 0; v < variables_size; ++v) {
    algorithm->registerRefine(
        getPatchDataId(Variable(v), Scratch(dir)), getPatchDataId(Variable(v)),
        getPatchDataId(Variable(v), Scratch(dir)), nullptr);
  }
  return algorithm;
}

std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
HyperbolicTimeIntegrator::getFillGhostLayerCoarsenAlgorithm(
    Direction dir) const {
  auto algorithm = std::make_shared<SAMRAI::xfer::CoarsenAlgorithm>(
      ideal_gas_->GetDimension());
  for (int v = 0; v < variables_size; ++v) {
    algorithm->registerCoarsen(getPatchDataId(Variable(v), Scratch(dir)),
                               getPatchDataId(Variable(v)), nullptr);
  }
  return algorithm;
}

} // namespace ideal_gas
} // namespace fub