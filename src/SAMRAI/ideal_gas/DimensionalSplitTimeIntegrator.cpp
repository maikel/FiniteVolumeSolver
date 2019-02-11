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

#include "fub/SAMRAI/ideal_gas/DimensionalSplitTimeIntegrator.hpp"
#include "fub/SAMRAI/ideal_gas/HlleGodunovMethod.hpp"

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
std::string GetScratchName_(int dir) {
  static constexpr std::string_view scratch = "scratch";
  std::string result{};
  result.reserve(scratch.size() + 3);
  result += scratch;
  result += "_";
  result += std::to_string(dir);
  return result;
}
#else
std::string GetScratchName_(int dir) {
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

std::string
getPrefixedFluxName_(const std::string& name,
                     DimensionalSplitTimeIntegrator::FluxVariable var) {
  static constexpr std::array<
      const char*, DimensionalSplitTimeIntegrator::kFluxVariablesSize>
      names{"/DensityFlux", "/MomentumFlux", "/EnergyFlux", "/SpeciesFlux"};
  return name + names[static_cast<int>(var)];
}

int getDepth_(DimensionalSplitTimeIntegrator::FluxVariable var, int dim,
              int specs) {
  switch (var) {
  case DimensionalSplitTimeIntegrator::FluxVariable::momentum:
    return dim;
  case DimensionalSplitTimeIntegrator::FluxVariable::species:
    return specs;
  default:
    return 1;
  }
}
} // namespace

DimensionalSplitTimeIntegrator::DimensionalSplitTimeIntegrator(
    std::shared_ptr<const IdealGasEquation> ideal_gas)
    : ideal_gas_{std::move(ideal_gas)},
      flux_method_{std::make_shared<samrai::ideal_gas::HlleGodunovMethod>()} {
  SAMRAI::hier::VariableDatabase& database =
      *SAMRAI::hier::VariableDatabase::getDatabase();
  const SAMRAI::tbox::Dimension& dim = ideal_gas_->GetDimension();
  const int dim_value = dim.getValue();
  const SAMRAI::hier::IntVector ghost_cell_width =
      flux_method_->GetStencilWidth(dim);
  for (int dir = 0; dir < dim_value; ++dir) {
    std::shared_ptr<SAMRAI::hier::VariableContext> scratch =
        database.getContext(GetScratchName_(dir));
    SAMRAI::hier::IntVector ghost = SAMRAI::hier::IntVector::getZero(dim);
    ghost[dir] = ghost_cell_width[dir];
    for (int v = 0; v < kVariablesSize; ++v) {
      scratch_[v + dir * kVariablesSize] =
          reassignContextOfDataId_(GetPatchDataId(Variable(v)), scratch, ghost);
    }
  }
  const std::string& name = ideal_gas_->GetName();
  const SAMRAI::hier::IntVector zeros = SAMRAI::hier::IntVector::getZero(dim);
  std::shared_ptr<SAMRAI::hier::VariableContext> current =
      database.getContext("current");
  for (int v = 0; v < kFluxVariablesSize; ++v) {
    const int n_species = ideal_gas_->GetNSpecies();
    const int depth = getDepth_(FluxVariable(v), dim_value, n_species);
    auto flux = getVariable_<SAMRAI::pdat::FaceVariable<double>>(
        dim, getPrefixedFluxName_(name, FluxVariable(v)), depth);
    face_[v] = database.registerVariableAndContext(flux, current, zeros);
  }
}

double DimensionalSplitTimeIntegrator::ComputeStableDtOnPatch(
    const SAMRAI::hier::Patch& patch, double /* time_point */,
    Direction dir) const {
  CompletePatchData state = GetCompletePatchData(patch, Scratch(dir));
  return flux_method_->ComputeStableDtOnPatch(state, patch, dir);
}

namespace {
void AdvanceTimeOnPatch_(
    const DimensionalSplitTimeIntegrator::CompletePatchData& state,
    const DimensionalSplitTimeIntegrator::FluxPatchData& flux,
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

void DimensionalSplitTimeIntegrator::AdvanceTimeOnPatch(
    const SAMRAI::hier::Patch& patch, double /* time_point */, double dt,
    Direction dir) const {
  FluxPatchData flux = GetFluxPatchData(patch);
  CompletePatchData state = GetCompletePatchData(patch, Scratch(dir));
  flux_method_->ComputeFluxesOnPatch(flux, state, patch, dt, dir);
  AdvanceTimeOnPatch_(state, flux, patch, dt, dir);
  CompletePatchData complete = GetCompletePatchData(patch);
  ConsPatchData cons = GetConsPatchData(patch, Scratch(dir));
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

DimensionalSplitTimeIntegrator::CompletePatchData
DimensionalSplitTimeIntegrator::GetCompletePatchData(
    const SAMRAI::hier::Patch& patch, Scratch scratch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, GetPatchDataId(Variable::density, scratch));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, GetPatchDataId(Variable::momentum, scratch));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, GetPatchDataId(Variable::energy, scratch));
  SAMRAI::pdat::CellData<double>& pressure = getCellData_(patch, GetPatchDataId(Variable::pressure, scratch));
  SAMRAI::pdat::CellData<double>& temperature = getCellData_(patch, GetPatchDataId(Variable::temperature, scratch));
  SAMRAI::pdat::CellData<double>& speed_of_sound = getCellData_(patch, GetPatchDataId(Variable::speed_of_sound, scratch));
  SAMRAI::pdat::CellData<double>& species = getCellData_(patch, GetPatchDataId(Variable::species, scratch));
  // clang-format on
  return {density,     momentum,       energy, pressure,
          temperature, speed_of_sound, species};
}

DimensionalSplitTimeIntegrator::CompleteSpans
DimensionalSplitTimeIntegrator::GetCompleteSpans(
    const SAMRAI::hier::Patch& patch, Scratch scratch) const {
  CompletePatchData state = GetCompletePatchData(patch, scratch);
  return boost::hana::make<StateTag<IdealGasEquation::Complete>>(
      MakeMdSpan<double, 3, layout_left>(state.density),
      MakeMdSpan<double, 4, layout_left>(state.momentum),
      MakeMdSpan<double, 3, layout_left>(state.energy),
      MakeMdSpan<double, 3, layout_left>(state.pressure),
      MakeMdSpan<double, 3, layout_left>(state.temperature),
      MakeMdSpan<double, 3, layout_left>(state.speed_of_sound),
      MakeMdSpan<double, 4, layout_left>(state.species));
}

DimensionalSplitTimeIntegrator::CompletePatchData
DimensionalSplitTimeIntegrator::GetCompletePatchData(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, GetPatchDataId(Variable::density));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, GetPatchDataId(Variable::momentum));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, GetPatchDataId(Variable::energy));
  SAMRAI::pdat::CellData<double>& pressure = getCellData_(patch, GetPatchDataId(Variable::pressure));
  SAMRAI::pdat::CellData<double>& temperature = getCellData_(patch, GetPatchDataId(Variable::temperature));
  SAMRAI::pdat::CellData<double>& speed_of_sound = getCellData_(patch, GetPatchDataId(Variable::speed_of_sound));
  SAMRAI::pdat::CellData<double>& species = getCellData_(patch, GetPatchDataId(Variable::species));
  // clang-format on
  return {density,     momentum,       energy, pressure,
          temperature, speed_of_sound, species};
}

DimensionalSplitTimeIntegrator::CompleteSpans
DimensionalSplitTimeIntegrator::GetCompleteSpans(
    const SAMRAI::hier::Patch& patch) const {
  CompletePatchData state = GetCompletePatchData(patch);
  return boost::hana::make<StateTag<IdealGasEquation::Complete>>(
      MakeMdSpan<double, 3, layout_left>(state.density),
      MakeMdSpan<double, 4, layout_left>(state.momentum),
      MakeMdSpan<double, 3, layout_left>(state.energy),
      MakeMdSpan<double, 3, layout_left>(state.pressure),
      MakeMdSpan<double, 3, layout_left>(state.temperature),
      MakeMdSpan<double, 3, layout_left>(state.speed_of_sound),
      MakeMdSpan<double, 4, layout_left>(state.species));
}

DimensionalSplitTimeIntegrator::ConsPatchData
DimensionalSplitTimeIntegrator::GetConsPatchData(
    const SAMRAI::hier::Patch& patch, Scratch scratch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, GetPatchDataId(Variable::density, scratch));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, GetPatchDataId(Variable::momentum, scratch));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, GetPatchDataId(Variable::energy, scratch));
  SAMRAI::pdat::CellData<double>& species = getCellData_(patch, GetPatchDataId(Variable::species, scratch));
  // clang-format on
  return {density, momentum, energy, species};
}

DimensionalSplitTimeIntegrator::ConsPatchData
DimensionalSplitTimeIntegrator::GetConsPatchData(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = getCellData_(patch, GetPatchDataId(Variable::density));
  SAMRAI::pdat::CellData<double>& momentum = getCellData_(patch, GetPatchDataId(Variable::momentum));
  SAMRAI::pdat::CellData<double>& energy = getCellData_(patch, GetPatchDataId(Variable::energy));
  SAMRAI::pdat::CellData<double>& species = getCellData_(patch, GetPatchDataId(Variable::species));
  // clang-format on
  return {density, momentum, energy, species};
}

DimensionalSplitTimeIntegrator::FluxPatchData
DimensionalSplitTimeIntegrator::GetFluxPatchData(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::FaceData<double>& density = getFaceData_(patch, GetPatchDataId(FluxVariable::density));
  SAMRAI::pdat::FaceData<double>& momentum = getFaceData_(patch, GetPatchDataId(FluxVariable::momentum));
  SAMRAI::pdat::FaceData<double>& energy = getFaceData_(patch, GetPatchDataId(FluxVariable::energy));
  SAMRAI::pdat::FaceData<double>& species = getFaceData_(patch, GetPatchDataId(FluxVariable::species));
  // clang-format on
  return {density, momentum, energy, species};
}

DimensionalSplitTimeIntegrator::FluxSpans
DimensionalSplitTimeIntegrator::GetFluxSpans(const SAMRAI::hier::Patch& patch,
                                             Direction face_normal) const {
  FluxPatchData fluxes = GetFluxPatchData(patch);
  const int face_n = int(face_normal);
  return boost::hana::make<StateTag<IdealGasEquation::Conservative>>(
      MakeMdSpan<double, 3, layout_left>(fluxes.density, face_n),
      MakeMdSpan<double, 4, layout_left>(fluxes.momentum, face_n),
      MakeMdSpan<double, 3, layout_left>(fluxes.energy, face_n),
      MakeMdSpan<double, 4, layout_left>(fluxes.species, face_n));
}

void DimensionalSplitTimeIntegrator::FlagPatchDataIdsToAllocate(
    SAMRAI::hier::ComponentSelector& which_to_allocate) const {
  for (int id : ideal_gas_->GetPatchDataIds()) {
    which_to_allocate.setFlag(id);
  }
  for (int id : GetScratch()) {
    which_to_allocate.setFlag(id);
  }
  for (int id : face_) {
    which_to_allocate.setFlag(id);
  }
}

std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>
DimensionalSplitTimeIntegrator::GetFillGhostLayerRefineAlgorithm(
    Direction dir) const {
  auto algorithm = std::make_shared<SAMRAI::xfer::RefineAlgorithm>();
  for (int v = 0; v < kVariablesSize; ++v) {
    algorithm->registerRefine(
        GetPatchDataId(Variable(v), Scratch(dir)), GetPatchDataId(Variable(v)),
        GetPatchDataId(Variable(v), Scratch(dir)), nullptr);
  }
  return algorithm;
}

std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
DimensionalSplitTimeIntegrator::GetFillGhostLayerCoarsenAlgorithm(
    Direction dir) const {
  auto algorithm = std::make_shared<SAMRAI::xfer::CoarsenAlgorithm>(
      ideal_gas_->GetDimension());
  for (int v = 0; v < kVariablesSize; ++v) {
    algorithm->registerCoarsen(GetPatchDataId(Variable(v), Scratch(dir)),
                               GetPatchDataId(Variable(v)), nullptr);
  }
  return algorithm;
}

} // namespace ideal_gas
} // namespace fub