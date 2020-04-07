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
#include <SAMRAI/pdat/SideVariable.h>

#include <SAMRAI/geom/CartesianCellDoubleWeightedAverage.h>
#include <SAMRAI/geom/CartesianSideDoubleWeightedAverage.h>
#include <SAMRAI/pdat/CellDoubleConstantRefine.h>

#include <range/v3/view/zip.hpp>

#include <numeric>

namespace fub::samrai {

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

IntegratorContext::AuxialiaryDataDescription
IntegratorContext::RegisterVariables(const DataDescription& desc, int gcw,
                                     int flux_gcw) {
  AuxialiaryDataDescription aux_desc;
  aux_desc.scratch_ids.reserve(desc.data_ids.size());
  aux_desc.flux_ids.reserve(desc.n_cons_variables);
  aux_desc.coarse_fine_ids.reserve(desc.n_cons_variables);
  SAMRAI::tbox::Dimension dim(desc.dim);
  const SAMRAI::hier::IntVector scratch_gcws(dim, gcw);
  const SAMRAI::hier::IntVector flux_gcws(dim, flux_gcw);
  using SAMRAI::hier::VariableContext;
  std::shared_ptr scratch = std::make_shared<VariableContext>("scratch");
  std::shared_ptr flux = std::make_shared<VariableContext>("flux");
  std::shared_ptr coarse_fine =
      std::make_shared<VariableContext>("coarse_fine_flux");
  SAMRAI::hier::VariableDatabase& vardb =
      *SAMRAI::hier::VariableDatabase::getDatabase();
  for (int id : desc.data_ids) {
    std::shared_ptr<SAMRAI::hier::Variable> variable{};
    if (vardb.mapIndexToVariable(id, variable)) {
      aux_desc.scratch_ids.push_back(
          vardb.registerVariableAndContext(variable, scratch, scratch_gcws)); 
    }
  }
  span<const int> cons_ids =
      span{desc.data_ids}.subspan(0, desc.n_cons_variables);
  for (int id : cons_ids) {
    std::shared_ptr<SAMRAI::hier::Variable> variable{};
    if (vardb.mapIndexToVariable(id, variable)) {
      const std::string& name = variable->getName();
      const int depth =
          static_cast<SAMRAI::pdat::CellVariable<double>*>(variable.get())
              ->getDepth();
      std::shared_ptr flux_variable =
          std::make_shared<SAMRAI::pdat::SideVariable<double>>(
              dim, name + "_On_Face", depth, true);
      aux_desc.flux_ids.push_back(
          vardb.registerVariableAndContext(flux_variable, flux, flux_gcws));
      aux_desc.coarse_fine_ids.push_back(vardb.registerVariableAndContext(
          flux_variable, coarse_fine, SAMRAI::hier::IntVector::getZero(dim)));
    }
  }
  return aux_desc;
}

namespace {

struct BoundaryConditionWrapper : SAMRAI::xfer::RefinePatchStrategy {
  AnyBoundaryCondition* condition_{};
  GriddingAlgorithm* grid_{};
  int level_number_{};

  void setPhysicalBoundaryConditions(
      SAMRAI::hier::Patch& patch, double,
      const SAMRAI::hier::IntVector&) override {
    if (condition_) {
      condition_->FillBoundary(patch, *grid_, level_number_);
    }
  }

  void preprocessRefine(SAMRAI::hier::Patch&, const SAMRAI::hier::Patch&,
                        const SAMRAI::hier::Box&,
                        const SAMRAI::hier::IntVector&) override {}

  void postprocessRefine(SAMRAI::hier::Patch&, const SAMRAI::hier::Patch&,
                         const SAMRAI::hier::Box&,
                         const SAMRAI::hier::IntVector&) override {}

  SAMRAI::hier::IntVector
  getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const override {
    return SAMRAI::hier::IntVector::getZero(dim);
  }
};

}

IntegratorContext::IntegratorContext(std::shared_ptr<GriddingAlgorithm> grid,
                                     HyperbolicMethod method,
                                     AuxialiaryDataDescription aux_desc)
    : ghost_cell_width_(method.flux_method.GetStencilWidth()),
      gridding_(std::move(grid)), method_(std::move(method)),
      aux_desc_(std::move(aux_desc)) {
  const std::size_t nlevel =
      std::size_t(gridding_->GetPatchHierarchy().GetMaxNumberOfLevels());
  time_points_.resize(nlevel);
  regrid_time_points_.resize(nlevel);
  cycles_.resize(nlevel);
  boundaries_.resize(nlevel);
  int level = 0;
  for (auto& boundary : boundaries_) {
    auto wrapper = std::make_shared<BoundaryConditionWrapper>();
    wrapper->condition_ = &gridding_->GetBoundaryCondition();
    wrapper->grid_ = gridding_.get();
    wrapper->level_number_ = level;
    level += 1;
    boundary = std::move(wrapper);
  }

  SAMRAI::tbox::Dimension dim(
      gridding_->GetPatchHierarchy().GetDataDescription().dim);
  fill_scratch_ = std::make_shared<SAMRAI::xfer::RefineAlgorithm>();
  coarsen_scratch_ = std::make_shared<SAMRAI::xfer::CoarsenAlgorithm>(dim);
  coarsen_fluxes_ = std::make_shared<SAMRAI::xfer::CoarsenAlgorithm>(dim);

  const DataDescription& desc =
      gridding_->GetPatchHierarchy().GetDataDescription();
  int k = 0;
  for (int dest : aux_desc_.scratch_ids) {
    fill_scratch_->registerRefine(
        dest, desc.data_ids[k], dest,
        std::make_shared<SAMRAI::pdat::CellDoubleConstantRefine>());
    coarsen_scratch_->registerCoarsen(
        dest, dest,
        std::make_shared<SAMRAI::geom::CartesianCellDoubleWeightedAverage>());
  }
  for (int flux_id : aux_desc_.flux_ids) {
    coarsen_fluxes_->registerCoarsen(
        flux_id, flux_id,
        std::make_shared<SAMRAI::geom::CartesianSideDoubleWeightedAverage>());
  }

  fill_scratch_two_level_schedule_.resize(nlevel);
  fill_scratch_one_level_schedule_.resize(nlevel);
  coarsen_scratch_schedule_.resize(nlevel);
  coarsen_fluxes_schedule_.resize(nlevel);

  ResetHierarchyConfiguration();
}

bool IntegratorContext::LevelExists(int level) const noexcept {
  return level <
         gridding_->GetPatchHierarchy().GetNative()->getNumberOfLevels();
}

void IntegratorContext::ResetHierarchyConfiguration(int level) {
  const PatchHierarchy& hier = gridding_->GetPatchHierarchy();
  const int nlevel = hier.GetNumberOfLevels();
  for (int ilvl = level; ilvl < nlevel; ++ilvl) {
    std::size_t lvl = std::size_t(ilvl);
    std::shared_ptr<SAMRAI::hier::PatchLevel> patch_level =
        hier.GetPatchLevel(ilvl);
    fill_scratch_one_level_schedule_[lvl] =
        fill_scratch_->createSchedule(patch_level, boundaries_[lvl].get());
    if (ilvl > 0) {
      std::shared_ptr<SAMRAI::hier::PatchLevel> prev_level =
        hier.GetPatchLevel(ilvl - 1);
      fill_scratch_two_level_schedule_[lvl] = fill_scratch_->createSchedule(
          patch_level, ilvl - 1, hier.GetNative(), boundaries_[lvl].get());
      coarsen_scratch_schedule_[lvl] =
          coarsen_scratch_->createSchedule(prev_level, patch_level);
      coarsen_fluxes_schedule_[lvl] =
          coarsen_fluxes_->createSchedule(prev_level, patch_level);
    }
  }
}

const SAMRAI::geom::CartesianGridGeometry&
IntegratorContext::GetGeometry(int level) const {
  return gridding_->GetPatchHierarchy().GetGeometry(level);
}

void IntegratorContext::ResetHierarchyConfiguration(
    std::shared_ptr<GriddingAlgorithm> gridding) {
  auto grid = gridding_;
  gridding_ = std::move(gridding);
  ResetHierarchyConfiguration();
}

SAMRAI::hier::IntVector IntegratorContext::GetRatioToCoarserLevel(int) const
    noexcept {
  return gridding_->GetPatchHierarchy().GetOptions().refine_ratio;
}

const PatchHierarchy& IntegratorContext::GetPatchHierarchy() const noexcept {
  return gridding_->GetPatchHierarchy();
}

PatchHierarchy& IntegratorContext::GetPatchHierarchy() noexcept {
  return gridding_->GetPatchHierarchy();
}

MPI_Comm IntegratorContext::GetMpiCommunicator() const noexcept {
  return gridding_->GetPatchHierarchy().GetNative()->getMPI().getCommunicator();
}

Duration IntegratorContext::GetTimePoint(int level) const {
  return time_points_[static_cast<std::size_t>(level)];
}

void IntegratorContext::SetTimePoint(Duration tp, int level) {
  time_points_[static_cast<std::size_t>(level)] = tp;
}

void IntegratorContext::SetCycles(std::ptrdiff_t cycles, int level) {
  cycles_[static_cast<std::size_t>(level)] = cycles;
}

SAMRAI::hier::PatchLevel& IntegratorContext::GetPatchLevel(int level) {
  return *gridding_->GetPatchHierarchy().GetPatchLevel(level);
}

const SAMRAI::hier::PatchLevel&
IntegratorContext::GetPatchLevel(int level) const {
  return *gridding_->GetPatchHierarchy().GetPatchLevel(level);
}

std::ptrdiff_t IntegratorContext::GetCycles(int level) const {
  return cycles_[static_cast<std::size_t>(level)];
}

span<const int> IntegratorContext::GetDataIds() const {
  return gridding_->GetPatchHierarchy().GetDataIds();
}

span<const int> IntegratorContext::GetScratchIds() const {
  return aux_desc_.scratch_ids;
}

span<const int> IntegratorContext::GetFluxIds() const {
  return aux_desc_.flux_ids;
}

int IntegratorContext::PreAdvanceLevel(int level_num,
                                        [[maybe_unused]] Duration dt,
                                        std::pair<int, int> subcycle) {
  if (subcycle.first == 0) {
    gridding_->RegridAllFinerLevels(level_num - 1);
    ResetHierarchyConfiguration(level_num);
    return level_num;
  }
  return 0;
}

Result<void, TimeStepTooLarge>
IntegratorContext::PostAdvanceLevel(int level_num, Duration dt, int subcycle) {
  SetCycles(GetCycles(level_num) + 1, level_num);
  double timepoint = (GetTimePoint(level_num) + dt).count();
  ::MPI_Bcast(&timepoint, 1, MPI_DOUBLE, 0, GetMpiCommunicator());
  SetTimePoint(Duration(timepoint), level_num);
  const std::shared_ptr<SAMRAI::hier::PatchLevel>& level =
      GetPatchHierarchy().GetPatchLevel(level_num);
  level->setTime(GetTimePoint(level_num).count());
  PatchHierarchy& hier = gridding_->GetPatchHierarchy();
  hier.SetTimePoint(GetTimePoint(level_num), level_num);
  hier.SetCycles(GetCycles(level_num), level_num);
  if (subcycle == 0) {
    span<const int> data = GetDataIds();
    span<const int> scratch = GetDataIds();
    SAMRAI::hier::PatchLevel& patches = GetPatchLevel(level_num);
    for (std::shared_ptr<SAMRAI::hier::Patch> patch : patches) {
      for (auto [dest, src] : ranges::views::zip(data, scratch)) {
        patch->getPatchData(dest)->copy(*patch->getPatchData(src));
      }
    }
  }
  return boost::outcome_v2::success();
}

void IntegratorContext::FillGhostLayerTwoLevels(int level,
                                                [[maybe_unused]] int coarse) {
  FUB_ASSERT(coarse == level - 1);
  const double fill_time = time_points_[level].count();
  fill_scratch_two_level_schedule_[level]->fillData(fill_time);
}

void IntegratorContext::FillGhostLayerSingleLevel(int level) {
  const double fill_time = time_points_[level].count();
  fill_scratch_one_level_schedule_[level]->fillData(fill_time);
}

void IntegratorContext::CoarsenConservatively(
    int fine_level, [[maybe_unused]] int coarse_level) {
  coarsen_scratch_schedule_[fine_level]->coarsenData();
}

const AnyBoundaryCondition&
IntegratorContext::GetBoundaryCondition() const {
  return gridding_->GetBoundaryCondition();
}

void IntegratorContext::ResetCoarseFineFluxes([[maybe_unused]] int fine,
                                              int /* coarse */) {
  // SAMRAI::hier::PatchLevel& patch_level = GetPatchLevel(fine);
  // span<const int> coarse_fine_ids = GetCoarseFineIds();
  // for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : fluxes) {
  //   for (int id : coarse_find_ids) {
  //     auto& data = static_cast<SAMRAI::pdat::SideData<double>&>(
  //         *patch->getPatchData(id));
  //     data.fillAll(0.0);
  //   }
  // }
}

void IntegratorContext::AccumulateCoarseFineFluxes(
    [[maybe_unused]] int level, [[maybe_unused]] double scale,
    [[maybe_unused]] Direction dir) {
  // SAMRAI::hier::PatchLevel& fluxes = GetPatchLevel(level);
  // const int ncomps = flux_desc_->getMaxNumberRegisteredComponents() / 2;
  // const int d = static_cast<int>(dir);
  // for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : fluxes) {
  //   for (int component = 0; component < ncomps; ++component) {
  //     auto* pointer_to_src =
  //         dynamic_cast<const SAMRAI::pdat::SideData<double>*>(
  //             patch->getPatchData(component).get());
  //     auto* pointer_to_dest = dynamic_cast<SAMRAI::pdat::SideData<double>*>(
  //         patch->getPatchData(component + ncomps).get());
  //     FUB_ASSERT(pointer_to_src && pointer_to_dest);
  //     span<const double> src =
  //         MakeMdSpan<4>(pointer_to_src->getArrayData(d)).get_span();
  //     span<double> dest =
  //         MakeMdSpan<4>(pointer_to_dest->getArrayData(d)).get_span();
  //     std::transform(src.begin(), src.end(), dest.begin(), dest.begin(),
  //                    [scale](double x, double y) { return y + scale * x; });
  //   }
  // }
}

void IntegratorContext::ApplyFluxCorrection(int fine,
                                            [[maybe_unused]] int coarse,
                                            [[maybe_unused]] fub::Duration dt) {
  coarsen_fluxes_schedule_[fine]->coarsenData();
}

void IntegratorContext::ComputeNumericFluxes(int level, Duration dt,
                                             Direction dir) {
  method_.flux_method.ComputeNumericFluxes(*this, level, dt, dir);
}

Duration IntegratorContext::ComputeStableDt(int level, Direction dir) {
  return method_.flux_method.ComputeStableDt(*this, level, dir);
}

void IntegratorContext::CompleteFromCons(int level, Duration dt) {
  method_.reconstruction.CompleteFromCons(*this, level, dt);
}

void IntegratorContext::UpdateConservatively(int level, Duration dt,
                                             Direction dir) {
  method_.time_integrator.UpdateConservatively(*this, level, dt, dir);
}

AnyBoundaryCondition& IntegratorContext::GetBoundaryCondition() {
  return gridding_->GetBoundaryCondition();
}

const std::shared_ptr<GriddingAlgorithm>&
IntegratorContext::GetGriddingAlgorithm() const noexcept {
  return gridding_;
}

} // namespace fub::samrai
