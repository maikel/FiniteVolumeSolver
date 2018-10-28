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

#include "fub/solver/euler/ForwardEulerTimeIntegrator.hpp"
#include "fub/solver/euler/HlleRiemannSolver.hpp"

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
                                 ForwardEulerTimeIntegrator::FluxVariable var) {
  static constexpr std::array<const char*,
                              ForwardEulerTimeIntegrator::flux_variables_size>
      names{"/DensityFlux", "/MomentumFlux", "/EnergyFlux"};
  return name + names[static_cast<int>(var)];
}
} // namespace

ForwardEulerTimeIntegrator::ForwardEulerTimeIntegrator(
    const IdealGas& ideal_gas)
    : ideal_gas_{ideal_gas}, flux_method_{
                                 std::make_shared<HlleRiemannSolver>()} {
  SAMRAI::hier::VariableDatabase& database =
      *SAMRAI::hier::VariableDatabase::getDatabase();
  const SAMRAI::tbox::Dimension& dim = ideal_gas_.getDim();
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

double ForwardEulerTimeIntegrator::computeStableDtOnPatch(
    const SAMRAI::hier::Patch& patch, double time_point) const {
  CompleteState state = getCompleteState(patch, Scratch(0));
  return flux_method_->computeStableDtOnPatch(state, patch);
}

namespace {
void advanceTimeOnPatch_(const ForwardEulerTimeIntegrator::CompleteState& state,
                         const ForwardEulerTimeIntegrator::FluxState& flux,
                         const SAMRAI::hier::Patch& patch, double dt,
                         Direction dir) {
  const double* dx = getCartesianPatchGeometry(patch)->getDx();
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
    FUB_ASSERT(state.energy(cell) > 0);
  }
}
} // namespace

void ForwardEulerTimeIntegrator::advanceTimeOnPatch(
    const SAMRAI::hier::Patch& patch, double time_point, double dt,
    Direction dir) const {
  FluxState flux = getFluxState(patch);
  CompleteState state = getCompleteState(patch, Scratch(dir));
  flux_method_->computeFluxesOnPatch(flux, state, patch, dt, dir);
  advanceTimeOnPatch_(state, flux, patch, dt, dir);
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

ForwardEulerTimeIntegrator::CompleteState
ForwardEulerTimeIntegrator::getCompleteState(const SAMRAI::hier::Patch& patch,
                                             Scratch scratch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, getPatchDataId(Variable::density, scratch));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, getPatchDataId(Variable::momentum, scratch));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, getPatchDataId(Variable::energy, scratch));
  SAMRAI::pdat::CellData<double>& pressure = getCellData_(patch, getPatchDataId(Variable::pressure, scratch));
  SAMRAI::pdat::CellData<double>& speed_of_sound = getCellData_(patch, getPatchDataId(Variable::speed_of_sound, scratch));
  // clang-format on
  return {density, momentum, energy, pressure, speed_of_sound};
}

ForwardEulerTimeIntegrator::CompleteState
ForwardEulerTimeIntegrator::getCompleteState(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, getPatchDataId(Variable::density));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, getPatchDataId(Variable::momentum));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, getPatchDataId(Variable::energy));
  SAMRAI::pdat::CellData<double>& pressure = getCellData_(patch, getPatchDataId(Variable::pressure));
  SAMRAI::pdat::CellData<double>& speed_of_sound = getCellData_(patch, getPatchDataId(Variable::speed_of_sound));
  // clang-format on
  return {density, momentum, energy, pressure, speed_of_sound};
}

ForwardEulerTimeIntegrator::ConsState
ForwardEulerTimeIntegrator::getConsState(const SAMRAI::hier::Patch& patch,
                                         Scratch scratch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, getPatchDataId(Variable::density, scratch));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, getPatchDataId(Variable::momentum, scratch));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, getPatchDataId(Variable::energy, scratch));
  // clang-format on
  return {density, momentum, energy};
}

ForwardEulerTimeIntegrator::ConsState ForwardEulerTimeIntegrator::getConsState(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, getPatchDataId(Variable::density));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, getPatchDataId(Variable::momentum));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, getPatchDataId(Variable::energy));
  // clang-format on
  return {density, momentum, energy};
}

ForwardEulerTimeIntegrator::FluxState ForwardEulerTimeIntegrator::getFluxState(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::FaceData<double>& density = getFaceData_(patch, getPatchDataId(FluxVariable::density));
  SAMRAI::pdat::FaceData<double>& momentum = getFaceData_(patch, getPatchDataId(FluxVariable::momentum));
  SAMRAI::pdat::FaceData<double>& energy = getFaceData_(patch, getPatchDataId(FluxVariable::energy));
  // clang-format on
  return {density, momentum, energy};
}

void ForwardEulerTimeIntegrator::flagPatchDataIdsToAllocate(
    SAMRAI::hier::ComponentSelector& which_to_allocate) const {
  for (int id : ideal_gas_.getDataIds()) {
    which_to_allocate.setFlag(id);
  }
  for (int id : scratch_) {
    which_to_allocate.setFlag(id);
  }
  for (int id : face_) {
    which_to_allocate.setFlag(id);
  }
}

std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>
ForwardEulerTimeIntegrator::getFillGhostLayerRefineAlgorithm(
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
ForwardEulerTimeIntegrator::getFillGhostLayerCoarsenAlgorithm(
    Direction dir) const {
  auto algorithm =
      std::make_shared<SAMRAI::xfer::CoarsenAlgorithm>(ideal_gas_.getDim());
  for (int v = 0; v < variables_size; ++v) {
    algorithm->registerCoarsen(getPatchDataId(Variable(v), Scratch(dir)),
                               getPatchDataId(Variable(v)), nullptr);
  }
  return algorithm;
}

} // namespace euler
} // namespace fub