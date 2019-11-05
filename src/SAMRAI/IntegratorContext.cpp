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

#include "fub/SAMRAI/IntegratorContext.hpp"

#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/pdat/SideDataFactory.h>

#include <SAMRAI/geom/CartesianCellDoubleWeightedAverage.h>
#include <SAMRAI/geom/CartesianSideDoubleWeightedAverage.h>
#include <SAMRAI/pdat/CellDoubleConstantRefine.h>

namespace fub::samrai {
std::shared_ptr<SAMRAI::hier::PatchDescriptor>
MakeScratchDescriptor(const DataDescription& data_description, int gcw) {
  SAMRAI::tbox::Dimension dim(data_description.dim);
  const SAMRAI::hier::IntVector ghosts(dim, gcw);
  std::shared_ptr descriptor =
      std::make_shared<SAMRAI::hier::PatchDescriptor>();
  SAMRAI::hier::VariableDatabase* vardb =
      SAMRAI::hier::VariableDatabase::getDatabase();
  for (int id : data_description.data_ids) {
    std::shared_ptr<SAMRAI::hier::Variable> variable{};
    if (vardb->mapIndexToVariable(id, variable)) {
      const std::string& name = variable->getName();
      const int depth =
          static_cast<SAMRAI::pdat::CellVariable<double>*>(variable.get())
              ->getDepth();
      descriptor->definePatchDataComponent(
          name, std::make_shared<SAMRAI::pdat::CellDataFactory<double>>(
                    depth, ghosts));
    }
  }
  return descriptor;
}

std::shared_ptr<SAMRAI::hier::PatchDescriptor>
MakeFluxDescriptor(const DataDescription& data_description, int gcw) {
  SAMRAI::tbox::Dimension dim(data_description.dim);
  const SAMRAI::hier::IntVector ghosts(dim, gcw);
  std::shared_ptr descriptor =
      std::make_shared<SAMRAI::hier::PatchDescriptor>();
  SAMRAI::hier::VariableDatabase* vardb =
      SAMRAI::hier::VariableDatabase::getDatabase();
  for (int id : data_description.data_ids) {
    std::shared_ptr<SAMRAI::hier::Variable> variable{};
    if (vardb->mapIndexToVariable(id, variable)) {
      const std::string& name = variable->getName();
      const int depth =
          static_cast<SAMRAI::pdat::CellVariable<double>*>(variable.get())
              ->getDepth();
      descriptor->definePatchDataComponent(
          name, std::make_shared<SAMRAI::pdat::SideDataFactory<double>>(
                    depth, ghosts, true));
    }
  }
  return descriptor;
}

IntegratorContext::LevelData::LevelData(
    std::shared_ptr<SAMRAI::hier::PatchLevel> level_in_hier,
    std::shared_ptr<SAMRAI::hier::PatchDescriptor> scratch_desc,
    std::shared_ptr<SAMRAI::hier::PatchDescriptor> flux_desc)
    : data(std::move(level_in_hier)),
      scratch(std::make_shared<SAMRAI::hier::PatchLevel>(
          data->getBoxLevel(), data->getGridGeometry(), scratch_desc)),
      fluxes(std::make_shared<SAMRAI::hier::PatchLevel>(
          data->getBoxLevel(), data->getGridGeometry(), flux_desc)) {
  scratch->allocatePatchData(SelectComponents(*scratch_desc));
  fluxes->allocatePatchData(SelectComponents(*flux_desc));
}

IntegratorContext::LevelData::LevelData(
    std::shared_ptr<SAMRAI::hier::PatchLevel> level_in_hier,
    const DataDescription& desc, int gcw)
    : LevelData(std::move(level_in_hier), MakeScratchDescriptor(desc, gcw),
                MakeFluxDescriptor(desc, gcw - 1)) {}

IntegratorContext::LevelData::LevelData(const LevelData& other)
    : LevelData(other.data, other.scratch->getPatchDescriptor(),
                other.fluxes->getPatchDescriptor()) {}

IntegratorContext::LevelData& IntegratorContext::LevelData::
operator=(const LevelData& other) {
  LevelData tmp(other);
  std::swap(*this, tmp);
  return *this;
}

IntegratorContext::IntegratorContext(std::shared_ptr<GriddingAlgorithm> grid,
                                     HyperbolicMethod method)
    : ghost_cell_width_(method.flux_method.GetStencilWidth()),
      gridding_(std::move(grid)), method_(std::move(method)),
      scratch_desc_(MakeScratchDescriptor(
          gridding_->GetPatchHierarchy().GetDataDescription(),
          ghost_cell_width_)),
      flux_desc_(MakeFluxDescriptor(
          gridding_->GetPatchHierarchy().GetDataDescription(),
          ghost_cell_width_)) {
  const std::size_t nlevel =
      std::size_t(gridding_->GetPatchHierarchy().GetMaxNumberOfLevels());
  level_data_.resize(nlevel);
  fill_scratch_ = std::make_shared<SAMRAI::xfer::RefineAlgorithm>();
  DataDescription desc = gridding_->GetPatchHierarchy().GetDataDescription();
  SAMRAI::tbox::Dimension dim(desc.dim);
  coarsen_scratch_ = std::make_shared<SAMRAI::xfer::CoarsenAlgorithm>(dim);
  coarsen_fluxes_ = std::make_shared<SAMRAI::xfer::CoarsenAlgorithm>(dim);

  int i = 0;
  for (int id : desc.data_ids) {
    fill_scratch_->registerRefine(
        i, id, i, std::make_shared<SAMRAI::pdat::CellDoubleConstantRefine>());
    coarsen_scratch_->registerCoarsen(
        i, i,
        std::make_shared<SAMRAI::geom::CartesianCellDoubleWeightedAverage>());
    i = i + 1;
  }
  for (int icons = 0; icons < desc.n_cons_variables; ++icons) {
    coarsen_fluxes_->registerCoarsen(
        icons, icons,
        std::make_shared<SAMRAI::geom::CartesianSideDoubleWeightedAverage>());
  }

  fill_scratch_two_level_schedule_.resize(nlevel);
  fill_scratch_one_level_schedule_.resize(nlevel);
  coarsen_scratch_schedule_.resize(nlevel);
  coarsen_fluxes_schedule_.resize(nlevel);
  ResetHierarchyConfiguration();
}

void IntegratorContext::ResetHierarchyConfiguration(int level) {
  const PatchHierarchy& hier = gridding_->GetPatchHierarchy();
  const int nlevel = hier.GetMaxNumberOfLevels();
  for (int ilvl = level; ilvl < nlevel; ++ilvl) {
    std::size_t lvl = std::size_t(ilvl);
    std::shared_ptr<SAMRAI::hier::PatchLevel> patch_level =
        hier.GetPatchLevel(ilvl);
    level_data_[lvl] = LevelData(patch_level, scratch_desc_, flux_desc_);
    fill_scratch_one_level_schedule_[lvl] =
        fill_scratch_->createSchedule(patch_level, &GetBoundaryCondition(ilvl));
    if (ilvl > 0) {
      fill_scratch_two_level_schedule_[lvl] = fill_scratch_->createSchedule(
          patch_level, ilvl - 1, hier.GetNative(), &GetBoundaryCondition(ilvl));
      std::shared_ptr<SAMRAI::hier::PatchLevel> coarse_level =
          hier.GetPatchLevel(ilvl - 1);
      coarsen_scratch_schedule_[lvl] =
          coarsen_scratch_->createSchedule(patch_level, coarse_level);
      coarsen_fluxes_schedule_[lvl] =
          coarsen_fluxes_->createSchedule(patch_level, coarse_level);
    }
  }
}

void IntegratorContext::PreAdvanceLevel(int level_num,
                                        [[maybe_unused]] Duration dt,
                                        int subcycle) {
  if (subcycle == 0) {
    gridding_->RegridAllFinerLevels(level_num - 1, gridding_->GetCycles(),
                                    gridding_->GetTimePoint());
    ResetHierarchyConfiguration(level_num);
  }
}

Result<void, TimeStepTooLarge>
IntegratorContext::PostAdvanceLevel(int level_num, Duration dt,
                                    [[maybe_unused]] int subcycle) {
  SetCycles(GetCycles(level_num) + 1, level_num);
  double timepoint = (GetTimePoint(level_num) + dt).count();
  ::MPI_Bcast(&timepoint, 1, MPI_DOUBLE, 0, GetMpiCommunicator());
  SetTimePoint(Duration(timepoint), level_num);
  const std::shared_ptr<SAMRAI::hier::PatchLevel>& level =
      GetPatchHierarchy().GetPatchLevel(level_num);
  level->setTime(GetTimePoint(level_num).count());
  return boost::outcome_v2::success();
}

void IntegratorContext::FillGhostLayerTwoLevels(int level,
                                                [[maybe_unused]] int coarse) {
  FUB_ASSERT(coarse == level - 1);
  const double fill_time = level_data_[level].time_point.count();
  fill_scratch_two_level_schedule_[level]->fillData(fill_time);
}

void IntegratorContext::FillGhostLayerSingleLevel(int level) {
  const double fill_time = level_data_[level].time_point.count();
  fill_scratch_one_level_schedule_[level]->fillData(fill_time);
}

void IntegratorContext::CoarsenConservatively(int fine_level,
                                              [[maybe_unused]] int coarse_level) {
  coarsen_scratch_schedule_[fine_level]->coarsenData();
}

const BoundaryCondition&
IntegratorContext::GetBoundaryCondition(int level) const {
  return gridding_->GetBoundaryCondition(level);
}

void IntegratorContext::ResetCoarseFineFluxes(int fine, int /* coarse */) {
  SAMRAI::hier::PatchLevel& fluxes = GetFluxes(fine);
  int ncomps = flux_desc_->getMaxNumberRegisteredComponents();
  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : fluxes) {
    for (int component = ncomps / 2; component < ncomps; ++component) {
      auto& data = static_cast<SAMRAI::pdat::SideData<double>&>(
          *patch->getPatchData(component));
      data.fillAll(0.0);
    }
  }
}

void IntegratorContext::AccumulateCoarseFineFluxes(int level, double time_scale,
                                                   Direction dir) {
  SAMRAI::hier::PatchLevel& fluxes = GetFluxes(level);
  const int ncomps = flux_desc_->getMaxNumberRegisteredComponents() / 2;
  const int d = static_cast<int>(dir);
  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : fluxes) {
    for (int component = 0; component < ncomps; ++component) {
      auto* pointer_to_src =
          dynamic_cast<const SAMRAI::pdat::SideData<double>*>(
              patch->getPatchData(component).get());
      auto* pointer_to_dest = dynamic_cast<SAMRAI::pdat::SideData<double>*>(
          patch->getPatchData(component + ncomps).get());
      FUB_ASSERT(pointer_to_src && pointer_to_dest);
      span<const double> src =
          MakeMdSpan<4>(pointer_to_src->getArrayData(d)).get_span();
      span<double> dest =
          MakeMdSpan<4>(pointer_to_dest->getArrayData(d)).get_span();
      std::transform(
          src.begin(), src.end(), dest.begin(), dest.begin(),
          [time_scale](double x, double y) { return y + time_scale * x; });
    }
  }
}

void IntegratorContext::ComputeNumericFluxes(int level, Duration dt,
                                             Direction dir) {
  method_.flux_method.ComputeNumericFluxes(*this, level, dt, dir);
}

void IntegratorContext::CompleteFromCons(int level, Duration dt) {
  method_.reconstruction.CompleteFromCons(*this, level, dt);
}

void IntegratorContext::UpdateConservatively(int level, Duration dt,
                                             Direction dir) {
  method_.time_integrator.UpdateConservatively(*this, level, dt, dir);
}

BoundaryCondition& IntegratorContext::GetBoundaryCondition(int level) {
  return gridding_->GetBoundaryCondition(level);
}

const std::shared_ptr<GriddingAlgorithm>&
IntegratorContext::GetGriddingAlgorithm() const noexcept {
  return gridding_;
}

} // namespace fub::samrai
