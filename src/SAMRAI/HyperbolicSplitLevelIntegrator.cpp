// Copyright (c) 2019 Maikel Nadolski
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

#include "fub/SAMRAI/HyperbolicSplitLevelIntegrator.hpp"

#include "fub/SAMRAI/CartesianPatchHierarchy.hpp"

#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/pdat/OutersideVariable.h>
#include <SAMRAI/pdat/SideGeometry.h>
#include <SAMRAI/pdat/SideVariable.h>

#include <SAMRAI/geom/CartesianOutersideDoubleWeightedAverage.h>

namespace fub {
namespace samrai {
void BoundaryCondition::setPhysicalBoundaryConditions(
    SAMRAI::hier::Patch& patch, double fill_time,
    const SAMRAI::hier::IntVector& ghost_width_to_fill) {
  const SAMRAI::hier::PatchGeometry& geometry = *patch.getPatchGeometry();
  const std::vector<SAMRAI::hier::BoundaryBox>& boundaries =
      geometry.getCodimensionBoundaries(1);
  for (const SAMRAI::hier::BoundaryBox& bbox : boundaries) {
    SAMRAI::hier::Box fill_box =
        geometry.getBoundaryFillBox(bbox, patch.getBox(), ghost_width_to_fill);
    if (fill_box.size() > 0) {
      const int direction = bbox.getLocationIndex() / 2;
      const int side = bbox.getLocationIndex() % 2;
      const int index = 2 * direction + side;
      base_->SetPhysicalBoundaryCondition(patch, fill_box, fill_time,
                                          Direction(direction), side);
    }
  }
}

namespace {
template <typename T>
constexpr span<const T*> AsConst(span<T*> pointers) noexcept {
  void* p = static_cast<void*>(pointers.data());
  return span<const T*>(static_cast<const T**>(p), pointers.size());
}

template <typename PatchData>
auto GetPatchData(const SAMRAI::hier::Patch& patch, span<const int> ids) {
  static constexpr int MaxVariables =
      HyperbolicSplitLevelIntegrator::MaxVariables;
  boost::container::static_vector<PatchData, MaxVariables> datas;
  // std::vector<PatchData> datas;
  for (int id : ids) {
    datas.push_back(static_cast<PatchData>(patch.getPatchData(id).get()));
  }
  return datas;
}

template <typename PatchData>
auto GetPatchData(const SAMRAI::hier::Patch& patch, span<const int> ids,
                  span<const int> subset) {
  static constexpr int MaxVariables =
      HyperbolicSplitLevelIntegrator::MaxVariables;
  boost::container::static_vector<PatchData, MaxVariables> datas;
  for (int pos : subset) {
    datas.push_back(static_cast<PatchData>(patch.getPatchData(ids[pos]).get()));
  }
  return datas;
}

SAMRAI::hier::IntVector UnitVector(SAMRAI::tbox::Dimension dim, int dir) {
  SAMRAI::hier::IntVector unit = SAMRAI::hier::IntVector::getZero(dim);
  unit[dir] = 1;
  return unit;
}

int RegisterIntermediateResult(int id) {
  SAMRAI::hier::VariableDatabase* vardb =
      SAMRAI::hier::VariableDatabase::getDatabase();
  std::shared_ptr<SAMRAI::hier::Variable> variable;
  if (!vardb->mapIndexToVariable(id, variable)) {
    throw std::invalid_argument(
        "Specified Index is not in SAMRAIs Variable Database.");
  }
  std::shared_ptr<SAMRAI::hier::VariableContext> context =
      vardb->getContext("intermediate");
  return vardb->registerVariableAndContext(
      variable, context, SAMRAI::hier::IntVector::getZero(variable->getDim()));
}

int RegisterScratch(int id, int dir, const SAMRAI::hier::IntVector& ghosts) {
  SAMRAI::hier::VariableDatabase* vardb =
      SAMRAI::hier::VariableDatabase::getDatabase();
  std::shared_ptr<SAMRAI::hier::Variable> variable;
  if (!vardb->mapIndexToVariable(id, variable)) {
    throw std::invalid_argument(
        "Specified Index is not in SAMRAIs Variable Database.");
  }
  const std::string name = "scratch" + std::to_string(dir);
  std::shared_ptr<SAMRAI::hier::VariableContext> context =
      vardb->getContext(name);
  return vardb->registerVariableAndContext(variable, context, ghosts);
}

int RegisterFlux(int id, int dir, const SAMRAI::hier::IntVector& ghosts) {
  SAMRAI::hier::VariableDatabase* vardb =
      SAMRAI::hier::VariableDatabase::getDatabase();
  std::shared_ptr<SAMRAI::hier::Variable> variable;
  if (!vardb->mapIndexToVariable(id, variable)) {
    throw std::invalid_argument(
        "Specified Index is not in SAMRAIs Variable Database.");
  }
  std::string variable_name = variable->getName();
  variable_name += "_flux";
  if (vardb->checkVariableExists(variable_name)) {
    variable = vardb->getVariable(variable_name);
  } else {
    variable = std::make_shared<SAMRAI::pdat::SideVariable<double>>(
        variable->getDim(), variable_name, dir);
  }
  const std::string name = "flux" + std::to_string(dir);
  std::shared_ptr<SAMRAI::hier::VariableContext> context =
      vardb->getContext(name);
  return vardb->registerVariableAndContext(variable, context, ghosts);
}

int RegisterOuterFlux(int id) {
  SAMRAI::hier::VariableDatabase* vardb =
      SAMRAI::hier::VariableDatabase::getDatabase();
  std::shared_ptr<SAMRAI::hier::Variable> variable;
  if (!vardb->mapIndexToVariable(id, variable)) {
    throw std::invalid_argument(
        "Specified Index is not in SAMRAIs Variable Database.");
  }
  std::string variable_name = variable->getName();
  variable_name += "_outer_flux";
  if (vardb->checkVariableExists(variable_name)) {
    variable = vardb->getVariable(variable_name);
  } else {
    variable = std::make_shared<SAMRAI::pdat::OutersideVariable<double>>(
        variable->getDim(), variable_name);
  }
  const std::string name = "outer_flux";
  std::shared_ptr<SAMRAI::hier::VariableContext> context =
      vardb->getContext(name);
  return vardb->registerVariableAndContext(
      variable, context, SAMRAI::hier::IntVector::getZero(variable->getDim()));
}

HyperbolicSplitLevelIntegrator::InternalDataIds
RegisterInternalDataIds(const PatchDataIdSet& desc,
                        const SAMRAI::hier::IntVector& stencil_width,
                        SAMRAI::tbox::Dimension dim) {
  // Register intermediate result data ids
  HyperbolicSplitLevelIntegrator::InternalDataIds data_ids{};
  data_ids.intermediate.reserve(desc.data_ids.size());
  for (int id : desc.data_ids) {
    data_ids.intermediate.push_back(RegisterIntermediateResult(id));
  }
  // Register scratch data ids (with ghost cell width)
  for (int dir = 0; dir < dim.getValue(); ++dir) {
    data_ids.scratch[dir].reserve(desc.data_ids.size());
    SAMRAI::hier::IntVector n_ghosts_cells = stencil_width;
    n_ghosts_cells *= 2;
    for (int id : desc.data_ids) {
      const int ghost_layer_width = n_ghosts_cells[dir];
      SAMRAI::hier::IntVector unit = UnitVector(dim, dir);
      unit *= ghost_layer_width;
      data_ids.scratch[dir].push_back(RegisterScratch(id, dir, unit));
    }
  }
  // Register outer flux data ids
  data_ids.outerside_fluxes.reserve(desc.conservative.size());
  for (int cons : desc.conservative) {
    data_ids.outerside_fluxes.push_back(RegisterOuterFlux(desc.data_ids[cons]));
  }
  // Register flux data ids
  for (int dir = 0; dir < dim.getValue(); ++dir) {
    data_ids.fluxes[dir].reserve(desc.conservative.size());
    for (int cons : desc.conservative) {
      SAMRAI::hier::IntVector unit = UnitVector(dim, dir);
      unit *= stencil_width[dir];
      data_ids.fluxes[dir].push_back(
          RegisterFlux(desc.data_ids[cons], dir, unit));
    }
  }
  return data_ids;
}

std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm> CoarsenOuterSideAlgorithm(
    SAMRAI::tbox::Dimension dim,
    const HyperbolicSplitLevelIntegrator::InternalDataIds& data_ids) {
  std::shared_ptr algorithm =
      std::make_shared<SAMRAI::xfer::CoarsenAlgorithm>(dim);
  for (int id : data_ids.outerside_fluxes) {
    algorithm->registerCoarsen(
        id, id,
        std::make_shared<
            SAMRAI::geom::CartesianOutersideDoubleWeightedAverage>());
  }
  return algorithm;
}

template <typename DataIds>
std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
CoarsenInnerRegionAlgorithm(SAMRAI::tbox::Dimension dim,
                            const DataIds& data_ids,
                            span<const int> conservative) {
  std::shared_ptr algorithm =
      std::make_shared<SAMRAI::xfer::CoarsenAlgorithm>(dim);
  for (int cons : conservative) {
    algorithm->registerCoarsen(
        data_ids.intermediate[cons], data_ids.intermediate[cons],
        std::make_shared<
            SAMRAI::geom::CartesianOutersideDoubleWeightedAverage>());
  }
  return algorithm;
}

std::vector<std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>>
MakeCoarsenSchedules(
    SAMRAI::xfer::CoarsenAlgorithm& algorithm,
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy) {
  const int n_levels = hierarchy->getNumberOfLevels();
  std::vector<std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> schedules;
  schedules.reserve(n_levels);
  schedules.push_back(nullptr);
  for (int lvl = 1; lvl < n_levels; ++lvl) {
    std::shared_ptr<SAMRAI::hier::PatchLevel> coarse =
        hierarchy->getPatchLevel(lvl - 1);
    std::shared_ptr<SAMRAI::hier::PatchLevel> fine =
        hierarchy->getPatchLevel(lvl);
    schedules.push_back(algorithm.createSchedule(coarse, fine));
  }
  return schedules;
}

std::array<std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>, 3>
FillGhostLayerAlgorithms(
    SAMRAI::tbox::Dimension dim,
    const HyperbolicSplitLevelIntegrator::InternalDataIds& data_ids) {
  std::array<std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>, 3> algorithms;
  const int n_components = data_ids.intermediate.size();
  for (int dir = 0; dir < dim.getValue(); ++dir) {
    algorithms[dir] = std::make_shared<SAMRAI::xfer::RefineAlgorithm>();
    for (int component = 0; component < n_components; ++component) {
      // note: nullptr means to only fill ghost cells and perform no refine
      // interpolation.
      algorithms[dir]->registerRefine(
          data_ids.scratch[dir][component], data_ids.intermediate[component],
          data_ids.scratch[dir][component], nullptr);
    }
  }
  return algorithms;
}

std::array<std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>, 3>
MakeRefineSchedules(
    const std::array<std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>, 3>&
        algorithms,
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    SAMRAI::xfer::RefinePatchStrategy& boundary_condition) {
  std::array<std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>, 3>
      schedules{};
  const int n_levels = hierarchy->getNumberOfLevels();
  for (int dir = 0; dir < hierarchy->getDim().getValue(); ++dir) {
    for (int level = 0; level < n_levels; ++level) {
      std::shared_ptr patch_level = hierarchy->getPatchLevel(level);
      schedules[dir].push_back(
          algorithms[dir]->createSchedule(patch_level, &boundary_condition));
    }
  }
  return schedules;
}

} // namespace

HyperbolicSplitLevelIntegrator::HyperbolicSplitLevelIntegrator(
    SAMRAI::tbox::Dimension dim, PatchDataIdSet description,
    HyperbolicSplitPatchIntegrator& patch_integrator,
    DimensionalSplitFluxMethod& flux_method,
    DimensionalSplitBoundaryCondition& boundary_condition)
    : description_{std::move(description)}, hierarchy_{nullptr},
      patch_integrator_{&patch_integrator}, flux_method_{&flux_method},
      boundary_condition_{&boundary_condition},
      data_ids_{RegisterInternalDataIds(
          description_, flux_method_->GetStencilWidth(dim), dim)},
      coarsen_outerside_algorithm_{CoarsenOuterSideAlgorithm(dim, data_ids_)},
      coarsen_inner_region_algorithm_{CoarsenInnerRegionAlgorithm(
          dim, data_ids_, description_.conservative)},
      fill_ghost_layer_algorithm_{FillGhostLayerAlgorithms(dim, data_ids_)} {}

void HyperbolicSplitLevelIntegrator::ResetHierarchyConfiguration(
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy) {
  hierarchy_ = std::move(hierarchy);
  coarsen_outerside_ =
      MakeCoarsenSchedules(*coarsen_outerside_algorithm_, hierarchy_);
  coarsen_inner_region_ =
      MakeCoarsenSchedules(*coarsen_inner_region_algorithm_, hierarchy_);
  fill_ghost_layer_ = MakeRefineSchedules(fill_ghost_layer_algorithm_,
                                          hierarchy_, boundary_condition_);
}

void HyperbolicSplitLevelIntegrator::AdvanceLevel(
    const SAMRAI::hier::PatchLevel& level, int direction, double dt) {
  // Do the normal conservative update on each patch of this level
  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
    // Compute fluxes in the specified direction
    auto scratch = GetScratch(*patch, direction);
    StaticVector<SAMRAI::pdat::SideData<double>*> fluxes =
        GetFluxes(*patch, direction);
    flux_method_->ComputeFluxesOnPatch(fluxes, AsConst(span(scratch)), *patch,
                                       direction, dt);

    // Use these fluxes to update cons variables at the "SCRATCH" context and
    // write their updated values to the "NEW" context.
    StaticVector<SAMRAI::pdat::CellData<double>*> scratch_cons =
        GetScratchCons(*patch, direction);
    StaticVector<SAMRAI::pdat::CellData<double>*> new_cons =
        GetNextCons(*patch);
    patch_integrator_->ConservativeUpdateOnPatch(
        new_cons, AsConst(span(fluxes)), AsConst(span(scratch_cons)), *patch,
        direction, dt);

    // Add the fluxes on outer sides to a specific storage
    StaticVector<SAMRAI::pdat::OutersideData<double>*> outerside_fluxes =
        GetOutersideFluxes(*patch);
    AccumulateFluxesOnOuterside(outerside_fluxes, AsConst(span(fluxes)), *patch,
                                direction, dt);
  }
  // If a finer level exists in the hierarchy, advance that finer level and use
  // the fine fluxes on coarse-fine interfaces
  const int level_num = level.getLevelNumber();
  if (hierarchy_->levelExists(level_num + 1)) {
    const SAMRAI::hier::PatchLevel& next_level =
        *hierarchy_->getPatchLevel(level_num + 1);
    ClearOutersideFluxes(next_level);
    const int refine_ratio = next_level.getRatioToCoarserLevel()[direction];
    for (int r = 0; r < refine_ratio; ++r) {
      FillGhostLayer(next_level, direction);
      AdvanceLevel(next_level, direction, dt / refine_ratio);
    }
    ConservativelyCoarsenOutersideFluxes(next_level, level, direction);
    ApplyCoarsenedFluxesOnLevel(next_level, level, dt, direction);
    ConservativelyCoarsenInnerRegions(next_level, level, direction);
  }
  // Interpolation happened conservatively on conservative variables.
  // Reconstruct a full state on each patch of the level.
  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
    StaticVector<SAMRAI::pdat::CellData<double>*> state = GetNextState(*patch);
    StaticVector<SAMRAI::pdat::CellData<double>*> cons = GetNextCons(*patch);
    patch_integrator_->ReconstructStatesFromCons(state, AsConst(span(cons)),
                                                 *patch, direction);
  }
}

HyperbolicSplitLevelIntegrator::StaticVector<SAMRAI::pdat::CellData<double>*>
HyperbolicSplitLevelIntegrator::GetScratch(const SAMRAI::hier::Patch& patch,
                                           int direction) const {
  return GetPatchData<SAMRAI::pdat::CellData<double>*>(
      patch, data_ids_.scratch[direction]);
}

HyperbolicSplitLevelIntegrator::StaticVector<SAMRAI::pdat::CellData<double>*>
HyperbolicSplitLevelIntegrator::GetScratchCons(const SAMRAI::hier::Patch& patch,
                                               int direction) const {
  return GetPatchData<SAMRAI::pdat::CellData<double>*>(
      patch, data_ids_.scratch[direction], description_.conservative);
}

HyperbolicSplitLevelIntegrator::StaticVector<SAMRAI::pdat::CellData<double>*>
HyperbolicSplitLevelIntegrator::GetNextState(
    const SAMRAI::hier::Patch& patch) const {
  return GetPatchData<SAMRAI::pdat::CellData<double>*>(patch,
                                                       data_ids_.intermediate);
}

HyperbolicSplitLevelIntegrator::StaticVector<SAMRAI::pdat::CellData<double>*>
HyperbolicSplitLevelIntegrator::GetNextCons(
    const SAMRAI::hier::Patch& patch) const {
  return GetPatchData<SAMRAI::pdat::CellData<double>*>(
      patch, data_ids_.intermediate, description_.conservative);
}

HyperbolicSplitLevelIntegrator::StaticVector<SAMRAI::pdat::SideData<double>*>
HyperbolicSplitLevelIntegrator::GetFluxes(const SAMRAI::hier::Patch& patch,
                                          int direction) const {
  return GetPatchData<SAMRAI::pdat::SideData<double>*>(
      patch, data_ids_.fluxes[direction]);
}

HyperbolicSplitLevelIntegrator::StaticVector<
    SAMRAI::pdat::OutersideData<double>*>
HyperbolicSplitLevelIntegrator::GetOutersideFluxes(
    const SAMRAI::hier::Patch& patch) const {
  return GetPatchData<SAMRAI::pdat::OutersideData<double>*>(
      patch, data_ids_.outerside_fluxes);
}

void HyperbolicSplitLevelIntegrator::AccumulateFluxesOnOuterside(
    span<SAMRAI::pdat::OutersideData<double>*> outerside,
    span<const SAMRAI::pdat::SideData<double>*> side,
    const SAMRAI::hier::Patch& patch, int direction, double dt) {
  FUB_ASSERT(outerside.size() == side.size());
  for (int i = 0; i < side.size(); ++i) {
    SAMRAI::pdat::ArrayData<double>& outer_lower =
        outerside[i]->getArrayData(direction, 0);
    SAMRAI::pdat::ArrayData<double>& outer_upper =
        outerside[i]->getArrayData(direction, 0);
    const SAMRAI::pdat::ArrayData<double>& all =
        side[i]->getArrayData(direction);
    SAMRAI::hier::Box intersection = outer_lower.getBox() * all.getBox();
    outer_lower.sum(all, intersection);
    intersection = outer_upper.getBox() * all.getBox();
    outer_upper.sum(all, intersection);
  }
}

void HyperbolicSplitLevelIntegrator::FillGhostLayer(
    const SAMRAI::hier::PatchLevel& level, int direction) {
  const int level_num = level.getLevelNumber();
  fill_ghost_layer_[direction][level_num]->fillData(0, 0.0);
}

void HyperbolicSplitLevelIntegrator::ClearOutersideFluxes(
    const SAMRAI::hier::PatchLevel& level) {
  // Loop over all patches
  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
    // Fetch outerside patch data objects for all flux variables
    StaticVector<SAMRAI::pdat::OutersideData<double>*> outersides =
        GetOutersideFluxes(*patch);
    // zero out outerside data for each variable
    for (SAMRAI::pdat::OutersideData<double>* outerside : outersides) {
      outerside->fill(0.0);
    }
  }
}

void HyperbolicSplitLevelIntegrator::ConservativelyCoarsenOutersideFluxes(
    const SAMRAI::hier::PatchLevel& fine_level,
    const SAMRAI::hier::PatchLevel& /* coarse_level */, int /* direction */) {
  FUB_ASSERT(fine_level.getLevelNumber() > 0);
  coarsen_outerside_[fine_level.getLevelNumber()]->coarsenData();
}

void HyperbolicSplitLevelIntegrator::ConservativelyCoarsenInnerRegions(
    const SAMRAI::hier::PatchLevel& fine_level,
    const SAMRAI::hier::PatchLevel& /* coarse_level */, int /* direction */) {
  FUB_ASSERT(fine_level.getLevelNumber() > 0);
  coarsen_inner_region_[fine_level.getLevelNumber()]->coarsenData();
}

namespace {
void ApplyCoarsenedFluxesOnPatch(
    SAMRAI::pdat::CellData<double>& states,
    const SAMRAI::pdat::SideData<double>& coarse_flux,
    const SAMRAI::pdat::OutersideData<double>& fine_flux, double lambda,
    int side, SAMRAI::pdat::SideIterator first,
    SAMRAI::pdat::SideIterator last) {
  while (first != last) {
    const SAMRAI::pdat::SideIndex index = *first;
    const SAMRAI::pdat::CellIndex cell(index.toCell((side + 1) % 2));
    const int sign = (side == 0) - (side != 0);
    states(cell) +=
        sign * lambda * (fine_flux(index, side) - coarse_flux(index));
    ++first;
  }
}

void ApplyCoarsenedFluxesOnPatch(
    SAMRAI::pdat::CellData<double>& states,
    const SAMRAI::pdat::SideData<double>& flux,
    const SAMRAI::pdat::OutersideData<double>& fine_flux, double lambda,
    int direction) {
  // Fix the lower side of the outerside data
  const SAMRAI::hier::Box lower_box =
      fine_flux.getArrayData(direction, 0).getBox();
  ApplyCoarsenedFluxesOnPatch(
      states, flux, fine_flux, lambda, 0,
      SAMRAI::pdat::SideGeometry::begin(lower_box, direction),
      SAMRAI::pdat::SideGeometry::end(lower_box, direction));
  // Fix the uper side of the outerside data
  const SAMRAI::hier::Box upper_box =
      fine_flux.getArrayData(direction, 1).getBox();
  ApplyCoarsenedFluxesOnPatch(
      states, flux, fine_flux, lambda, 1,
      SAMRAI::pdat::SideGeometry::begin(upper_box, direction),
      SAMRAI::pdat::SideGeometry::end(upper_box, direction));
}
} // namespace

void HyperbolicSplitLevelIntegrator::ApplyCoarsenedFluxesOnLevel(
    const SAMRAI::hier::PatchLevel& /* fine_level */,
    const SAMRAI::hier::PatchLevel& coarse_level, double dt, int direction) {
  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : coarse_level) {
    StaticVector<SAMRAI::pdat::OutersideData<double>*> outer_flux =
        GetOutersideFluxes(*patch);
    StaticVector<SAMRAI::pdat::SideData<double>*> flux =
        GetFluxes(*patch, direction);
    StaticVector<SAMRAI::pdat::CellData<double>*> states = GetNextState(*patch);
    const double dx = GetCartesianPatchGeometry(*patch)->getDx()[direction];
    const double lambda = dt / dx;
    for (int i = 0; i < states.size(); ++i) {
      ApplyCoarsenedFluxesOnPatch(*states[i], *flux[i], *outer_flux[i], lambda,
                                  direction);
    }
  }
}

} // namespace samrai
} // namespace fub