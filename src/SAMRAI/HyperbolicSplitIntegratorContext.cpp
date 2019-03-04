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

#include "fub/grid/SAMRAI/HyperbolicSplitIntegratorContext.hpp"

#include "fub/grid/SAMRAI/CartesianPatchHierarchy.hpp"

#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/pdat/OutersideVariable.h>
#include <SAMRAI/pdat/SideGeometry.h>
#include <SAMRAI/pdat/SideVariable.h>

#include <SAMRAI/geom/CartesianCellDoubleWeightedAverage.h>
#include <SAMRAI/geom/CartesianOutersideDoubleWeightedAverage.h>
#include <SAMRAI/pdat/CellDoubleConstantRefine.h>

#include <fmt/format.h>

namespace fub {
namespace samrai {
namespace {
template <typename PatchData>
auto GetPatchData(const SAMRAI::hier::Patch& patch, span<const int> ids) {
  static constexpr int MaxVariables =
      HyperbolicSplitIntegratorContext::MaxVariables;
  boost::container::static_vector<PatchData, MaxVariables> datas;
  for (int id : ids) {
    FUB_ASSERT(dynamic_cast<PatchData>(patch.getPatchData(id).get()));
    datas.push_back(static_cast<PatchData>(patch.getPatchData(id).get()));
  }
  return datas;
}

template <typename PatchData>
auto GetPatchData(const SAMRAI::hier::Patch& patch, span<const int> ids,
                  span<const int> subset) {
  static constexpr int MaxVariables =
      HyperbolicSplitIntegratorContext::MaxVariables;
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

// int RegisterIntermediateResult(int id) {
//   SAMRAI::hier::VariableDatabase* vardb =
//       SAMRAI::hier::VariableDatabase::getDatabase();
//   std::shared_ptr<SAMRAI::hier::Variable> variable;
//   if (!vardb->mapIndexToVariable(id, variable)) {
//     throw std::invalid_argument(
//         "Specified Index is not in SAMRAIs Variable Database.");
//   }
//   std::shared_ptr<SAMRAI::hier::VariableContext> context =
//       vardb->getContext("intermediate");
//   return vardb->registerVariableAndContext(
//       variable, context,
//       SAMRAI::hier::IntVector::getZero(variable->getDim()));
// }

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
  const int depth =
      static_cast<SAMRAI::pdat::CellVariable<double>*>(variable.get())
          ->getDepth();
  variable_name += "_flux" + std::to_string(dir);
  SAMRAI::tbox::Dimension dim = variable->getDim();
  const SAMRAI::hier::IntVector directions = UnitVector(dim, dir);
  if (vardb->checkVariableExists(variable_name)) {
    variable = vardb->getVariable(variable_name);
  } else {
    variable = std::make_shared<SAMRAI::pdat::SideVariable<double>>(
        dim, variable_name, directions, depth);
  }
  const std::string name = "flux";
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
  SAMRAI::tbox::Dimension dim = variable->getDim();
  std::string variable_name = variable->getName();
  const int depth =
      static_cast<SAMRAI::pdat::CellVariable<double>*>(variable.get())
          ->getDepth();
  variable_name += "_outer_flux";
  variable = std::make_shared<SAMRAI::pdat::OutersideVariable<double>>(
      dim, variable_name, depth);
  const std::string context_name = "flux";
  std::shared_ptr<SAMRAI::hier::VariableContext> context =
      vardb->getContext(context_name);
  return vardb->registerVariableAndContext(
      variable, context, SAMRAI::hier::IntVector::getZero(dim));
}

HyperbolicSplitIntegratorContext::InternalDataIds
RegisterInternalDataIds(const DataDescription& desc,
                        const SAMRAI::hier::IntVector& stencil_width,
                        SAMRAI::tbox::Dimension dim) {
  // Register intermediate result data ids
  HyperbolicSplitIntegratorContext::InternalDataIds data_ids{};
  data_ids.intermediate = desc.data_ids;
  // for (int id : desc.data_ids) {
  //   data_ids.intermediate.push_back(RegisterIntermediateResult(id));
  // }
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
  // Register outerside flux data ids for coarse fine interfaces
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

std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm> CoarsenOutersideAlgorithm(
    SAMRAI::tbox::Dimension dim,
    const HyperbolicSplitIntegratorContext::InternalDataIds& data_ids) {
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
                            span<const int> conservative, int dir) {
  std::shared_ptr algorithm =
      std::make_shared<SAMRAI::xfer::CoarsenAlgorithm>(dim);
  for (int cons : conservative) {
    algorithm->registerCoarsen(
        data_ids.scratch[dir][cons], data_ids.scratch[dir][cons],
        std::make_shared<SAMRAI::geom::CartesianCellDoubleWeightedAverage>());
  }
  return algorithm;
}

template <typename DataIds>
std::array<std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>,
           SAMRAI_MAXIMUM_DIMENSION>
CoarsenInnerRegionAlgorithmArray(SAMRAI::tbox::Dimension dim,
                                 const DataIds& data_ids,
                                 span<const int> conservative) {
  std::array<std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>, 3> algorithm{};
  for (int dir = 0; dir < dim.getValue(); ++dir) {
    algorithm[dir] =
        CoarsenInnerRegionAlgorithm(dim, data_ids, conservative, dir);
  }
  return algorithm;
}

void UpdateCoarsenSchedules(
    span<std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> schedules, int level,
    SAMRAI::xfer::CoarsenAlgorithm& algorithm,
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy) {
  const int n_levels = hierarchy->getNumberOfLevels();
  level = std::max(level, 1);
  // FUB_ASSERT(level < n_levels);
  for (int lvl = level; lvl < n_levels; ++lvl) {
    std::shared_ptr<SAMRAI::hier::PatchLevel> coarse =
        hierarchy->getPatchLevel(lvl - 1);
    std::shared_ptr<SAMRAI::hier::PatchLevel> fine =
        hierarchy->getPatchLevel(lvl);
    schedules[lvl] = algorithm.createSchedule(coarse, fine);
  }
}

std::array<std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>,
           SAMRAI_MAXIMUM_DIMENSION>
FillGhostLayerAlgorithms(
    SAMRAI::tbox::Dimension dim,
    const HyperbolicSplitIntegratorContext::InternalDataIds& data_ids) {
  std::array<std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>,
             SAMRAI_MAXIMUM_DIMENSION>
      algorithms;
  const int n_components = data_ids.intermediate.size();
  for (int dir = 0; dir < dim.getValue(); ++dir) {
    algorithms[dir] = std::make_shared<SAMRAI::xfer::RefineAlgorithm>();
    for (int component = 0; component < n_components; ++component) {
      algorithms[dir]->registerRefine(
          data_ids.scratch[dir][component], data_ids.intermediate[component],
          data_ids.scratch[dir][component],
          std::make_shared<SAMRAI::pdat::CellDoubleConstantRefine>());
    }
  }
  return algorithms;
}

std::array<std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>,
           SAMRAI_MAXIMUM_DIMENSION>
MakeRefineSchedulesTwoLevels(
    const std::array<std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>,
                     SAMRAI_MAXIMUM_DIMENSION>& algorithms,
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    SAMRAI::xfer::RefinePatchStrategy& boundary_condition) {
  std::array<std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>,
             SAMRAI_MAXIMUM_DIMENSION>
      schedules{};
  const int n_levels = hierarchy->getNumberOfLevels();
  for (int dir = 0; dir < hierarchy->getDim().getValue(); ++dir) {
    schedules[dir].reserve(n_levels);
    for (int level = 0; level < n_levels; ++level) {
      std::shared_ptr patch_level = hierarchy->getPatchLevel(level);
      schedules[dir].push_back(algorithms[dir]->createSchedule(
          patch_level, level - 1, hierarchy, &boundary_condition));
    }
  }
  return schedules;
}

std::array<std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>,
           SAMRAI_MAXIMUM_DIMENSION>
MakeRefineSchedulesSingleLevel(
    const std::array<std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>,
                     SAMRAI_MAXIMUM_DIMENSION>& algorithms,
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    SAMRAI::xfer::RefinePatchStrategy& boundary_condition) {
  std::array<std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>,
             SAMRAI_MAXIMUM_DIMENSION>
      schedules{};
  const int n_levels = hierarchy->getNumberOfLevels();
  for (int dir = 0; dir < hierarchy->getDim().getValue(); ++dir) {
    schedules[dir].reserve(n_levels);
    for (int level = 0; level < n_levels; ++level) {
      std::shared_ptr patch_level = hierarchy->getPatchLevel(level);
      schedules[dir].push_back(
          algorithms[dir]->createSchedule(patch_level, &boundary_condition));
    }
  }
  return schedules;
}

void SelectComponents(SAMRAI::hier::ComponentSelector& selector,
                      const std::vector<int>& data_ids) {
  for (int id : data_ids) {
    selector.setFlag(id);
  }
}

SAMRAI::hier::ComponentSelector SelectComponents(
    SAMRAI::tbox::Dimension dim,
    const HyperbolicSplitIntegratorContext::InternalDataIds& data_ids) {
  SAMRAI::hier::ComponentSelector selector;
  SelectComponents(selector, data_ids.intermediate);
  SelectComponents(selector, data_ids.outerside_fluxes);
  for (int dir = 0; dir < dim.getValue(); ++dir) {
    SelectComponents(selector, data_ids.scratch[dir]);
    SelectComponents(selector, data_ids.fluxes[dir]);
  }
  return selector;
}
} // namespace

void HyperbolicSplitIntegratorContext::AdaptBoundaryCondition::
    setPhysicalBoundaryConditions(
        SAMRAI::hier::Patch& /* patch */, double /* fill_time */,
        const SAMRAI::hier::IntVector& /* ghost_width_to_fill */) {
  // if (condition_) {
  //   const SAMRAI::hier::PatchGeometry& geometry = *patch.getPatchGeometry();
  //   const std::vector<SAMRAI::hier::BoundaryBox>& boundaries =
  //       geometry.getCodimensionBoundaries(1);
  //   for (const SAMRAI::hier::BoundaryBox& bbox : boundaries) {
  //     SAMRAI::hier::Box fill_box = geometry.getBoundaryFillBox(
  //         bbox, patch.getBox(), ghost_width_to_fill);
  //     if (fill_box.size() > 0) {
  //       const int direction = bbox.getLocationIndex() / 2;
  //       const int side = bbox.getLocationIndex() % 2;
  //       const int index = 2 * direction + side;
  //       auto states = context_->GetScratch(&patch, Direction(direction));
  //       (*condition_)(&patch, fill_box, fill_time, Direction(direction), side);
  //     }
  //   }
  // }
}

HyperbolicSplitIntegratorContext::HyperbolicSplitIntegratorContext(
    const GriddingAlgorithm& gridding, DataDescription description, int gcw)
    : gridding_{gridding}, description_{std::move(description)},
      data_ids_{RegisterInternalDataIds(
          description_,
          SAMRAI::hier::IntVector(gridding.hierarchy->getDim(), gcw),
          gridding.hierarchy->getDim())},
      ghost_cell_width_(GetPatchHierarchy()->getDim(), 2 * gcw),
      coarsen_outerside_algorithm_{
          CoarsenOutersideAlgorithm(GetPatchHierarchy()->getDim(), data_ids_)},
      coarsen_inner_region_algorithm_{CoarsenInnerRegionAlgorithmArray(
          GetPatchHierarchy()->getDim(), data_ids_, description_.conservative)},
      fill_ghost_layer_algorithm_{
          FillGhostLayerAlgorithms(GetPatchHierarchy()->getDim(), data_ids_)},
      adapted_boundary_{this} {
  ResetHierarchyConfiguration();
}

void HyperbolicSplitIntegratorContext::ResetHierarchyConfiguration(
    int level_num) {
  // allocate memory on all levels
  const int n_levels = GetPatchHierarchy()->getNumberOfLevels();
  SAMRAI::hier::ComponentSelector components =
      SelectComponents(GetPatchHierarchy()->getDim(), data_ids_);
  for (int lvl = level_num; lvl < n_levels; ++lvl) {
    SAMRAI::hier::PatchLevel& level = *GetPatchHierarchy()->getPatchLevel(lvl);
    // level.allocatePatchData(components);
    for (int comp = 0; comp < components.getSize(); ++comp) {
      if (components.isSet(comp) && !level.checkAllocated(comp)) {
        level.allocatePatchData(comp);
      }
    }
  }
  // Recreate Schedules.
  coarsen_outerside_.resize(GetPatchHierarchy()->getNumberOfLevels());
  UpdateCoarsenSchedules(coarsen_outerside_, level_num,
                         *coarsen_outerside_algorithm_, GetPatchHierarchy());
  for (int dir = 0; dir < GetPatchHierarchy()->getDim().getValue(); ++dir) {
    coarsen_inner_region_[dir].resize(GetPatchHierarchy()->getNumberOfLevels());
    UpdateCoarsenSchedules(coarsen_inner_region_[dir], level_num,
                           *coarsen_inner_region_algorithm_[dir],
                           GetPatchHierarchy());
  }
  fill_ghost_single_level_ = MakeRefineSchedulesSingleLevel(
      fill_ghost_layer_algorithm_, GetPatchHierarchy(), adapted_boundary_);
  fill_ghost_two_levels_ = MakeRefineSchedulesTwoLevels(
      fill_ghost_layer_algorithm_, GetPatchHierarchy(), adapted_boundary_);
  for (int i = 0; i < 3; ++i) {
    time_points_[i].resize(GetPatchHierarchy()->getNumberOfLevels());
    regrid_time_points_[i].resize(GetPatchHierarchy()->getNumberOfLevels());
    cycles_[i].resize(GetPatchHierarchy()->getNumberOfLevels());
  }
}

bool HyperbolicSplitIntegratorContext::LevelExists(int level) const noexcept {
  return level < GetPatchHierarchy()->getNumberOfLevels();
}

Duration HyperbolicSplitIntegratorContext::GetTimePoint(int level,
                                                        Direction dir) const {
  return time_points_[int(dir)][level];
}

void HyperbolicSplitIntegratorContext::SetTimePoint(Duration t, int level,
                                                    Direction dir) {
  time_points_[int(dir)][level] = t;
}

std::ptrdiff_t
HyperbolicSplitIntegratorContext::GetCycles(int level, Direction dir) const {
  return cycles_[int(dir)][level];
}

void HyperbolicSplitIntegratorContext::SetCycles(std::ptrdiff_t n, int level,
                                                 Direction dir) {
  cycles_[int(dir)][level] = n;
}

int HyperbolicSplitIntegratorContext::GetRatioToCoarserLevel(int level) const
    noexcept {
  if (level) {
    return 2;
  }
  return 1;
}

int HyperbolicSplitIntegratorContext::GetGhostCellWidth(PatchHandle,
                                                        Direction dir) {
  return ghost_cell_width_[int(dir)];
}

MPI_Comm HyperbolicSplitIntegratorContext::GetMpiCommunicator() const noexcept {
  return GetPatchHierarchy()->getMPI().getCommunicator();
}

HyperbolicSplitIntegratorContext::vector<SAMRAI::pdat::CellData<double>*>
HyperbolicSplitIntegratorContext::GetScratch(PatchHandle patch,
                                             Direction direction) {
  const int d = static_cast<int>(direction);
  return GetPatchData<SAMRAI::pdat::CellData<double>*>(*patch,
                                                       data_ids_.scratch[d]);
}

HyperbolicSplitIntegratorContext::vector<SAMRAI::pdat::CellData<double>*>
HyperbolicSplitIntegratorContext::GetData(PatchHandle patch) {
  return GetPatchData<SAMRAI::pdat::CellData<double>*>(*patch,
                                                       data_ids_.intermediate);
}

HyperbolicSplitIntegratorContext::vector<SAMRAI::pdat::SideData<double>*>
HyperbolicSplitIntegratorContext::GetFluxes(PatchHandle patch,
                                            Direction direction) {
  const int d = static_cast<int>(direction);
  return GetPatchData<SAMRAI::pdat::SideData<double>*>(*patch,
                                                       data_ids_.fluxes[d]);
}

HyperbolicSplitIntegratorContext::vector<SAMRAI::pdat::OutersideData<double>*>
HyperbolicSplitIntegratorContext::GetOutersideFluxes(PatchHandle patch) const {
  return GetPatchData<SAMRAI::pdat::OutersideData<double>*>(
      *patch, data_ids_.outerside_fluxes);
}

void HyperbolicSplitIntegratorContext::AccumulateCoarseFineFluxes(int level_num,
                                                                  Direction dir,
                                                                  Duration) {
  const SAMRAI::hier::PatchLevel& level =
      *GetPatchHierarchy()->getPatchLevel(level_num);
  const int direction = static_cast<int>(dir);
  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
    vector<SAMRAI::pdat::OutersideData<double>*> outerside =
        GetOutersideFluxes(patch.get());
    vector<SAMRAI::pdat::SideData<double>*> face = GetFluxes(patch.get(), dir);
    FUB_ASSERT(outerside.size() == face.size());
    for (std::size_t i = 0; i < face.size(); ++i) {
      SAMRAI::pdat::ArrayData<double>& outer_lower =
          outerside[i]->getArrayData(direction, 0);
      SAMRAI::pdat::ArrayData<double>& outer_upper =
          outerside[i]->getArrayData(direction, 1);
      const SAMRAI::pdat::ArrayData<double>& all_faces =
          face[i]->getArrayData(direction);
      const int ncomp = all_faces.getDepth();
      for (int d = 0; d < ncomp; ++d) {
        for (const SAMRAI::hier::Index& index : outer_lower.getBox()) {
          outer_lower(index, d) += 0.5 * all_faces(index, d);
        }
        for (const SAMRAI::hier::Index& index : outer_upper.getBox()) {
          outer_upper(index, d) += 0.5 * all_faces(index, d);
        }
      }
    }
  }
}

void HyperbolicSplitIntegratorContext::FillGhostLayerTwoLevels(
    int fine, int /* coarse */, Direction direction, BoundaryCondition boundary) {
  const int d = static_cast<int>(direction);
  adapted_boundary_.SetBoundaryCondition(boundary);
  fill_ghost_two_levels_[d][fine]->fillData(0.0);
}

void HyperbolicSplitIntegratorContext::FillGhostLayerSingleLevel(
    int level, Direction direction, BoundaryCondition boundary) {
  const int d = static_cast<int>(direction);
  adapted_boundary_.SetBoundaryCondition(boundary);
  fill_ghost_single_level_[d][level]->fillData(0.0);
}

void HyperbolicSplitIntegratorContext::ResetCoarseFineFluxes(int fine,
                                                             int /* coarse */,
                                                             Direction /* dir */) {
  // Loop over all patches
  const SAMRAI::hier::PatchLevel& level = *GetPatchHierarchy()->getPatchLevel(fine);
  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
    // Fetch outerside patch data objects for all flux variables
    vector<SAMRAI::pdat::OutersideData<double>*> outersides =
        GetOutersideFluxes(patch.get());
    // zero out outerside data for each variable
    for (SAMRAI::pdat::OutersideData<double>* outerside : outersides) {
      outerside->fill(0.0);
    }
  }
}

void HyperbolicSplitIntegratorContext::CoarsenConservatively(int fine, int,
                                                             Direction dir) {
  const int d = static_cast<int>(dir);
  coarsen_inner_region_[d][fine]->coarsenData();
}

const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&
HyperbolicSplitIntegratorContext::GetPatchHierarchy() const noexcept {
  return gridding_.hierarchy;
}

double HyperbolicSplitIntegratorContext::GetDx(PatchHandle patch,
                                               Direction dir) const {
  SAMRAI::geom::CartesianPatchGeometry* geom =
      GetCartesianPatchGeometry(*patch);
  const double* dx = geom->getDx();
  return dx[int(dir)];
}

CartesianCoordinates HyperbolicSplitIntegratorContext::GetCartesianCoordinates(
    PatchHandle patch) const {
  return ::fub::samrai::GetCartesianCoordinates(*patch);
}

void HyperbolicSplitIntegratorContext::ApplyFluxCorrection(int fine, int,
                                                           Duration,
                                                           Direction) {
  coarsen_outerside_[fine]->coarsenData();
}

void HyperbolicSplitIntegratorContext::PreAdvanceLevel(int level_num,
                                                       Direction dir,
                                                       Duration /* dt */,
                                                       int subcycle) {
  const int d = int(dir);
  if (subcycle == 0 && level_num > 0 &&
      regrid_time_points_[d][level_num] != time_points_[d][level_num]) {
    gridding_.RegridAllFinerLevels(level_num - 1, GetCycles(level_num - 1, dir),
                                   GetTimePoint(level_num - 1, dir));
    for (std::size_t lvl = level_num; lvl < regrid_time_points_.size(); ++lvl) {
      regrid_time_points_[d][level_num] = time_points_[d][level_num];
    }
    ResetHierarchyConfiguration(level_num);
    ResetCoarseFineFluxes(level_num, level_num - 1, dir);
  }
}

void HyperbolicSplitIntegratorContext::PostAdvanceLevel(int level_num,
                                                        Direction dir,
                                                        Duration dt, int) {
  SetCycles(GetCycles(level_num, dir) + 1, level_num, dir);
  SetTimePoint(GetTimePoint(level_num, dir) + dt, level_num, dir);
}

} // namespace samrai
} // namespace fub