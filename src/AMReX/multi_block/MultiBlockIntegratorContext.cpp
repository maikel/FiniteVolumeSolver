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

#include "fub/AMReX/multi_block/MultiBlockIntegratorContext.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

namespace fub::amrex {

MultiBlockIntegratorContext::MultiBlockIntegratorContext(
    FlameMasterReactor reactor, std::vector<IntegratorContext> tubes,
    std::vector<cutcell::IntegratorContext> plena,
    std::vector<BlockConnection> connectivity)
    : tubes_{std::move(tubes)}, plena_{std::move(plena)} {
  std::vector<std::shared_ptr<GriddingAlgorithm>> tube_grids(tubes_.size());
  std::transform(
      tubes_.begin(), tubes_.end(), tube_grids.begin(),
      [](IntegratorContext& ctx) { return ctx.GetGriddingAlgorithm(); });
  std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plenum_grids(
      plena_.size());
  std::transform(plena_.begin(), plena_.end(), plenum_grids.begin(),
                 [](cutcell::IntegratorContext& ctx) {
                   return ctx.GetGriddingAlgorithm();
                 });
  gridding_ = std::make_shared<MultiBlockGriddingAlgorithm>(
      std::move(reactor), std::move(tube_grids), std::move(plenum_grids),
      std::move(connectivity));
}

MultiBlockIntegratorContext::MultiBlockIntegratorContext(
    FlameMasterReactor reactor, std::vector<IntegratorContext> tubes,
    std::vector<cutcell::IntegratorContext> plena,
    std::vector<BlockConnection> connectivity,
    std::vector<std::shared_ptr<PressureValve>> valves)
    : tubes_{std::move(tubes)}, plena_{std::move(plena)} {
  std::vector<std::shared_ptr<GriddingAlgorithm>> tube_grids(tubes_.size());
  std::transform(
      tubes_.begin(), tubes_.end(), tube_grids.begin(),
      [](IntegratorContext& ctx) { return ctx.GetGriddingAlgorithm(); });
  std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plenum_grids(
      plena_.size());
  std::transform(plena_.begin(), plena_.end(), plenum_grids.begin(),
                 [](cutcell::IntegratorContext& ctx) {
                   return ctx.GetGriddingAlgorithm();
                 });
  gridding_ = std::make_shared<MultiBlockGriddingAlgorithm>(
      std::move(reactor), std::move(tube_grids), std::move(plenum_grids),
      std::move(connectivity), std::move(valves));
}

void MultiBlockIntegratorContext::CopyDataToScratch(int level) {
  for (IntegratorContext& ctx : tubes_) {
    ctx.CopyDataToScratch(level);
  }
  for (cutcell::IntegratorContext& ctx : plena_) {
    ctx.CopyDataToScratch(level);
  }
}

void MultiBlockIntegratorContext::CopyScratchToData(int level) {
  for (IntegratorContext& ctx : tubes_) {
    ctx.CopyScratchToData(level);
  }
  for (cutcell::IntegratorContext& ctx : plena_) {
    ctx.CopyScratchToData(level);
  }
}

span<IntegratorContext> MultiBlockIntegratorContext::Tubes() noexcept {
  return tubes_;
}
span<const IntegratorContext> MultiBlockIntegratorContext::Tubes() const
    noexcept {
  return tubes_;
}

span<cutcell::IntegratorContext> MultiBlockIntegratorContext::Plena() noexcept {
  return plena_;
}
span<const cutcell::IntegratorContext>
MultiBlockIntegratorContext::Plena() const noexcept {
  return plena_;
}

const std::shared_ptr<MultiBlockGriddingAlgorithm>&
MultiBlockIntegratorContext::GetGriddingAlgorithm() const noexcept {
  return gridding_;
}

Duration MultiBlockIntegratorContext::GetTimePoint(int level) const {
  return plena_[0].GetTimePoint(level);
}

std::ptrdiff_t MultiBlockIntegratorContext::GetCycles(int level) const {
  return plena_[0].GetCycles(level);
}

MPI_Comm MultiBlockIntegratorContext::GetMpiCommunicator() const noexcept {
  return plena_[0].GetMpiCommunicator();
}

const std::shared_ptr<CounterRegistry>&
MultiBlockIntegratorContext::GetCounterRegistry() const noexcept {
  return plena_[0].GetCounterRegistry();
}

bool MultiBlockIntegratorContext::LevelExists(int level) const noexcept {
  return std::any_of(tubes_.begin(), tubes_.end(),
                     [=](const IntegratorContext& ctx) {
                       return ctx.LevelExists(level);
                     }) ||
         std::any_of(plena_.begin(), plena_.end(),
                     [=](const cutcell::IntegratorContext& ctx) {
                       return ctx.LevelExists(level);
                     });
}

int MultiBlockIntegratorContext::GetRatioToCoarserLevel(int level,
                                                        Direction dir) const
    noexcept {
  return plena_[0].GetRatioToCoarserLevel(level, dir);
}

::amrex::IntVect
MultiBlockIntegratorContext::GetRatioToCoarserLevel(int level) const noexcept {
  return plena_[0].GetRatioToCoarserLevel(level);
}

template <typename T, typename S, typename F>
void ForEachBlock(const std::tuple<T, S>& tuple, F function) {
  for (auto&& x : std::get<0>(tuple)) {
    function(x);
  }
  for (auto&& x : std::get<1>(tuple)) {
    function(x);
  }
}

/// \brief Replaces the underlying gridding algorithm with the specified one.
void MultiBlockIntegratorContext::ResetHierarchyConfiguration(
    std::shared_ptr<MultiBlockGriddingAlgorithm> gridding) {
  {
    auto new_tube = gridding->GetTubes().begin();
    for (IntegratorContext& ctx : tubes_) {
      ctx.ResetHierarchyConfiguration(*new_tube);
      ++new_tube;
    }
  }
  {
    auto new_plenum = gridding->GetPlena().begin();
    for (cutcell::IntegratorContext& ctx : plena_) {
      ctx.ResetHierarchyConfiguration(*new_plenum);
      ++new_plenum;
    }
  }
  gridding_ = gridding;
}

/// \brief Whenever the gridding algorithm changes the data hierarchy this
/// function will regrid all distributed helper variables managed by the
/// context.
///
/// \param[in] level  The level number of the coarsest level which changed its
/// shape. Regrid all levels finer than level.
void MultiBlockIntegratorContext::ResetHierarchyConfiguration(int level) {
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}}, [level](auto& block) {
    if (block.LevelExists(level)) {
      block.ResetHierarchyConfiguration(level);
    }
  });
}

/// \brief Sets the cycle count for a specific level number and direction.
void MultiBlockIntegratorContext::SetCycles(std::ptrdiff_t cycle, int level) {
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}},
               [cycle, level](auto& block) {
                 if (block.LevelExists(level)) {
                   block.SetCycles(cycle, level);
                 }
               });
}

/// \brief Sets the time point for a specific level number and direction.
void MultiBlockIntegratorContext::SetTimePoint(Duration t, int level) {
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}}, [t, level](auto& block) {
    if (block.LevelExists(level)) {
      block.SetTimePoint(t, level);
    }
  });
}
/// @}

/// @{
/// \name Member functions relevant for the level integrator algorithm.

/// \brief On each first subcycle this will regrid the data if neccessary.
int MultiBlockIntegratorContext::PreAdvanceLevel(int level_num, Duration dt,
                                                 std::pair<int, int> subcycle) {
  int level_which_changed = std::numeric_limits<int>::max();
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}},
               [&level_which_changed, level_num, dt, subcycle](auto& block) {
                 if (block.LevelExists(level_num)) {
                   level_which_changed =
                       std::min(level_which_changed,
                                block.PreAdvanceLevel(level_num, dt, subcycle));
                 }
               });
  MultiBlockBoundary* boundary = gridding_->GetBoundaries(level_num).begin();
  for (const BlockConnection& connection : gridding_->GetConnectivity()) {
    IntegratorContext& tube = tubes_[connection.tube.id];
    cutcell::IntegratorContext& plenum = plena_[connection.plenum.id];
    boundary->ComputeBoundaryData(plenum.GetPatchHierarchy(),
                                  tube.GetPatchHierarchy());
    ++boundary;
  }
  if (level_which_changed ==
      plena_[0].GetPatchHierarchy().GetMaxNumberOfLevels()) {
    return 1;
  }
}

/// \brief Increases the internal time stamps and cycle counters for the
/// specified level number and direction.
Result<void, TimeStepTooLarge>
MultiBlockIntegratorContext::PostAdvanceLevel(int level_num, Duration dt,
                                              std::pair<int, int> subcycle) {
  for (IntegratorContext& tube : tubes_) {
    if (tube.LevelExists(level_num)) {
      Result<void, TimeStepTooLarge> result =
          tube.PostAdvanceLevel(level_num, dt, subcycle);
      if (!result) {
        return result;
      }
    }
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    if (plenum.LevelExists(level_num)) {
      Result<void, TimeStepTooLarge> result =
          plenum.PostAdvanceLevel(level_num, dt, subcycle);
      if (!result) {
        return result;
      }
    }
  }
  return boost::outcome_v2::success();
}

namespace {
struct WrapBoundaryCondition {
  std::size_t id;
  span<const BlockConnection> connectivity;
  span<MultiBlockBoundary> boundaries;
  AnyBoundaryCondition* base_condition{nullptr};
  cutcell::AnyBoundaryCondition* base_cc_condition{nullptr};

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level) const {
    FUB_ASSERT(base_condition != nullptr);
    base_condition->FillBoundary(mf, gridding, level);
    MultiBlockBoundary* boundary = boundaries.begin();
    for (const BlockConnection& connection : connectivity) {
      if (id == connection.tube.id) {
        boundary->FillBoundary(mf, gridding, level);
      }
      ++boundary;
    }
  }
  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir) const {
    FUB_ASSERT(base_condition != nullptr);
    base_condition->FillBoundary(mf, gridding, level, dir);
    MultiBlockBoundary* boundary = boundaries.begin();
    for (const BlockConnection& connection : connectivity) {
      if (id == connection.tube.id) {
        boundary->FillBoundary(mf, gridding, level, dir);
      }
      ++boundary;
    }
  }

  void FillBoundary(::amrex::MultiFab& mf,
                    const cutcell::GriddingAlgorithm& gridding,
                    int level) const {
    FUB_ASSERT(base_cc_condition != nullptr);
    base_cc_condition->FillBoundary(mf, gridding, level);
    MultiBlockBoundary* boundary = boundaries.begin();
    for (const BlockConnection& connection : connectivity) {
      if (id == connection.plenum.id) {
        boundary->FillBoundary(mf, gridding, level);
      }
      ++boundary;
    }
  }
  void FillBoundary(::amrex::MultiFab& mf,
                    const cutcell::GriddingAlgorithm& gridding, int level,
                    Direction dir) const {
    FUB_ASSERT(base_cc_condition != nullptr);
    base_cc_condition->FillBoundary(mf, gridding, level, dir);
    MultiBlockBoundary* boundary = boundaries.begin();
    for (const BlockConnection& connection : connectivity) {
      if (id == connection.plenum.id) {
        boundary->FillBoundary(mf, gridding, level, dir);
      }
      ++boundary;
    }
  }
};
} // namespace

void MultiBlockIntegratorContext::ApplyBoundaryCondition(int level,
                                                         Direction dir) {
  std::size_t id = 0;
  if (dir == Direction::X) {
    for (IntegratorContext& tube : tubes_) {
      if (tube.LevelExists(level)) {
        AnyBoundaryCondition& bc =
            tube.GetGriddingAlgorithm()->GetBoundaryCondition();
        AnyBoundaryCondition wrapped = WrapBoundaryCondition{
            id, gridding_->GetConnectivity(), gridding_->GetBoundaries(level),
            &bc, nullptr};
        tube.ApplyBoundaryCondition(level, dir, wrapped);
      }
      id += 1;
    }
  }
  id = 0;
  for (cutcell::IntegratorContext& plenum : plena_) {
    if (plenum.LevelExists(level)) {
      cutcell::AnyBoundaryCondition& bc =
          plenum.GetGriddingAlgorithm()->GetBoundaryCondition();
      cutcell::AnyBoundaryCondition wrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(level), nullptr, &bc};
      plenum.ApplyBoundaryCondition(level, dir, wrapped);
    }
    id += 1;
  }
}

/// \brief Fills the ghost layer of the scratch data and interpolates in the
/// coarse fine layer.
void MultiBlockIntegratorContext::FillGhostLayerTwoLevels(int fine,
                                                          int coarse) {
  std::size_t id = 0;
  for (IntegratorContext& tube : tubes_) {
    if (tube.LevelExists(fine)) {
      AnyBoundaryCondition& fbc =
          tube.GetGriddingAlgorithm()->GetBoundaryCondition();
      AnyBoundaryCondition fwrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(fine), &fbc, nullptr};
      AnyBoundaryCondition& cbc =
          tube.GetGriddingAlgorithm()->GetBoundaryCondition();
      AnyBoundaryCondition cwrapped = WrapBoundaryCondition{
          id, gridding_->GetConnectivity(), gridding_->GetBoundaries(coarse),
          &cbc, nullptr};
      tube.FillGhostLayerTwoLevels(fine, fwrapped, coarse, cwrapped);
    }
    id += 1;
  }
  id = 0;
  for (cutcell::IntegratorContext& plenum : plena_) {
    if (plenum.LevelExists(fine)) {
      cutcell::AnyBoundaryCondition& fbc =
          plenum.GetGriddingAlgorithm()->GetBoundaryCondition();
      cutcell::AnyBoundaryCondition fwrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(fine), nullptr, &fbc};
      cutcell::AnyBoundaryCondition& cbc =
          plenum.GetGriddingAlgorithm()->GetBoundaryCondition();
      cutcell::AnyBoundaryCondition cwrapped = WrapBoundaryCondition{
          id, gridding_->GetConnectivity(), gridding_->GetBoundaries(coarse),
          nullptr, &cbc};
      plenum.FillGhostLayerTwoLevels(fine, fwrapped, coarse, cwrapped);
    }
    id += 1;
  }
}

/// \brief Fills the ghost layer of the scratch data and does nothing in the
/// coarse fine layer.
void MultiBlockIntegratorContext::FillGhostLayerSingleLevel(int level) {
  std::size_t id = 0;
  for (IntegratorContext& tube : tubes_) {
    if (tube.LevelExists(level)) {
      AnyBoundaryCondition& bc =
          tube.GetGriddingAlgorithm()->GetBoundaryCondition();
      AnyBoundaryCondition wrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(level), &bc, nullptr};
      tube.FillGhostLayerSingleLevel(level, wrapped);
    }
    id += 1;
  }
  id = 0;
  for (cutcell::IntegratorContext& plenum : plena_) {
    if (plenum.LevelExists(level)) {
      cutcell::AnyBoundaryCondition& bc =
          plenum.GetGriddingAlgorithm()->GetBoundaryCondition();
      cutcell::AnyBoundaryCondition wrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(level), nullptr, &bc};
      plenum.FillGhostLayerSingleLevel(level, wrapped);
    }
    id += 1;
  }
}

/// \brief Returns a estimate for a stable time step size which can be taken
/// for specified level number in direction dir.
Duration MultiBlockIntegratorContext::ComputeStableDt(int level,
                                                      Direction dir) {
  Duration stable_dt{std::numeric_limits<double>::infinity()};
  if (dir == Direction::X) {
    stable_dt = std::accumulate(
        tubes_.begin(), tubes_.end(), stable_dt,
        [=](Duration dt, IntegratorContext& ctx) {
          return ctx.LevelExists(level)
                     ? std::min(dt, ctx.ComputeStableDt(level, dir))
                     : dt;
        });
  }
  stable_dt = std::accumulate(
      plena_.begin(), plena_.end(), stable_dt,
      [=](Duration dt, cutcell::IntegratorContext& ctx) {
        return ctx.LevelExists(level)
                   ? std::min(dt, ctx.ComputeStableDt(level, dir))
                   : dt;
      });
  return stable_dt;
}

void MultiBlockIntegratorContext::PreAdvanceHierarchy() {
  for (cutcell::IntegratorContext& ctx : plena_) {
    ctx.PreAdvanceHierarchy();
  }
}

void MultiBlockIntegratorContext::PostAdvanceHierarchy() {
  for (cutcell::IntegratorContext& ctx : plena_) {
    ctx.PostAdvanceHierarchy();
  }
}

namespace {
::amrex::Box GetBoundary(const ::amrex::Box& box, Direction dir, int side) {
  const int d = static_cast<int>(dir);
  ::amrex::IntVect lo = box.smallEnd();
  ::amrex::IntVect hi = box.bigEnd();
  if (side == 0) {
    hi[d] = lo[d];
  } else {
    lo[d] = hi[d];
  }
  return {lo, hi, ::amrex::IndexType({AMREX_D_DECL(1, 0, 0)})};
}

void AverageFlux_(span<double> average, const ::amrex::MultiFab& fluxes,
                  const ::amrex::EBFArrayBoxFactory& factory,
                  const ::amrex::Box& box, Direction dir, int side,
                  MPI_Comm comm) {
  std::fill(average.begin(), average.end(), 0.0);
  const int d = static_cast<int>(dir);
  double local_total_area = 0.0;
  const ::amrex::MultiCutFab& betas = *factory.getAreaFrac()[d];
  const ::amrex::Box boundary = GetBoundary(box, dir, side);
  ForEachFab(fluxes, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FabType type = factory.getMultiEBCellFlagFab()[mfi].getType();
    if (type == ::amrex::FabType::singlevalued) {
      const ::amrex::CutFab& area = betas[mfi];
      ForEachIndex(mfi.tilebox() & boundary, [&](auto... is) {
        ::amrex::IntVect iv{static_cast<int>(is)...};
        local_total_area += area(iv);
      });
    } else if (type == ::amrex::FabType::regular) {
      ForEachIndex(mfi.tilebox() & boundary,
                   [&](auto...) { local_total_area += 1.0; });
    }
  });
  double total_area = 0.0;
  ::MPI_Allreduce(&local_total_area, &total_area, 1, MPI_DOUBLE, MPI_SUM, comm);

  if (total_area > 0.0) {
    ForEachFab(fluxes, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::MultiCutFab& betas = *factory.getAreaFrac()[d];
      const ::amrex::FArrayBox& fab = fluxes[mfi];
      ::amrex::FabType type = factory.getMultiEBCellFlagFab()[mfi].getType();
      if (type == ::amrex::FabType::singlevalued) {
        const ::amrex::CutFab& area = betas[mfi];
        for (int component = 0; component < average.size(); ++component) {
          ForEachIndex(mfi.tilebox() & boundary, [&](auto... is) {
            ::amrex::IntVect iv{static_cast<int>(is)...};
            average[component] += area(iv) * fab(iv, component);
          });
        }
      } else if (type == ::amrex::FabType::regular) {
        for (int component = 0; component < average.size(); ++component) {
          ForEachIndex(mfi.tilebox() & boundary, [&](auto... is) {
            ::amrex::IntVect iv{static_cast<int>(is)...};
            average[component] += fab(iv, component);
          });
        }
      }
    });
    std::vector<double> buffer(average.size());
    ::MPI_Allreduce(average.data(), buffer.data(), average.size(), MPI_DOUBLE,
                    MPI_SUM, comm);
    std::transform(buffer.begin(), buffer.end(), average.begin(),
                   [=](double flux_sum) { return flux_sum / total_area; });
  }
}

template <typename IntegratorContext>
::amrex::Box RefineCoarseBoxToLevel_(const IntegratorContext& context,
                                     const ::amrex::Box& coarse_box,
                                     int level) {
  ::amrex::Box refined_box = coarse_box;
  for (int ilevel = 1; ilevel <= level; ++ilevel) {
    refined_box.refine(context.GetRatioToCoarserLevel(ilevel));
  }
  return refined_box;
}
} // namespace

/// \brief Fill the flux MultiFab with numeric fluxes based on current states
/// in scratch.
void MultiBlockIntegratorContext::ComputeNumericFluxes(int level, Duration dt,
                                                       Direction dir) {
  if (dir == Direction::X) {
    for (IntegratorContext& tube : tubes_) {
      if (tube.LevelExists(level)) {
        tube.ComputeNumericFluxes(level, dt, dir);
      }
    }
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    if (plenum.LevelExists(level)) {
      plenum.ComputeNumericFluxes(level, dt, dir);
    }
  }
  std::vector<double> average_flux(
      plena_[0].GetFluxes(0, Direction::X).nComp());
  span<const MultiBlockBoundary> boundaries = gridding_->GetBoundaries(level);
  const MultiBlockBoundary* boundary = boundaries.begin();
  for (const BlockConnection& conn : gridding_->GetConnectivity()) {
    if (dir == conn.direction &&
        !(boundary->GetValve() &&
          boundary->GetValve()->state == PressureValveState::closed)) {
      cutcell::IntegratorContext& plenum = plena_[conn.plenum.id];
      IntegratorContext& tube = tubes_[conn.tube.id];
      const ::amrex::EBFArrayBoxFactory& factory =
          plenum.GetEmbeddedBoundary(level);
      const ::amrex::MultiFab& fluxes = plenum.GetFluxes(level, conn.direction);
      ::amrex::Box plenum_mirror_box_on_level =
          RefineCoarseBoxToLevel_(plenum, conn.plenum.mirror_box, level);
      AverageFlux_(average_flux, fluxes, factory, plenum_mirror_box_on_level,
                   conn.direction, conn.side, plenum.GetMpiCommunicator());
      ::amrex::MultiFab& tube_fluxes = tube.GetFluxes(level, Direction::X);
      ::amrex::Box faces = RefineCoarseBoxToLevel_(
          tube, ::amrex::surroundingNodes(conn.tube.mirror_box, 0), level);
      ::amrex::IntVect iv = conn.side == 0 ? faces.bigEnd() : faces.smallEnd();
      ForEachFab(tube_fluxes, [&](const ::amrex::MFIter& mfi) {
        const ::amrex::Box box = mfi.tilebox();
        if (box.contains(iv)) {
          ::amrex::FArrayBox& fab = tube_fluxes[mfi];
          fab(iv, 0) = average_flux[0];
          fab(iv, 1) = average_flux[1];
          int tube_component = 2;
          int plenum_component = 1 + AMREX_SPACEDIM;
          while (tube_component < tube_fluxes.nComp()) {
            fab(iv, tube_component) = average_flux[plenum_component];
            ++tube_component;
            ++plenum_component;
          }
        }
      });
    }
    ++boundary;
  }
}

/// \brief Apply a conservative time update for each conservative variable on
/// the specified level number and direction.
void MultiBlockIntegratorContext::UpdateConservatively(int level, Duration dt,
                                                       Direction dir) {
  if (dir == Direction::X) {
    for (IntegratorContext& tube : tubes_) {
      if (tube.LevelExists(level)) {
        tube.UpdateConservatively(level, dt, dir);
      }
    }
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    if (plenum.LevelExists(level)) {
      plenum.UpdateConservatively(level, dt, dir);
    }
  }
}

/// \brief Reconstruct complete state variables from conservative ones.
void MultiBlockIntegratorContext::CompleteFromCons(int level, Duration dt) {
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}}, [=](auto& block) {
    if (block.LevelExists(level)) {
      block.CompleteFromCons(level, dt);
    }
  });
}

/// \brief Accumulate fluxes on the coarse fine interfaces for a specified
/// fine level number.
void MultiBlockIntegratorContext::AccumulateCoarseFineFluxes(int level,
                                                             double time_scale,
                                                             Direction dir) {
  if (dir == Direction::X) {
    for (IntegratorContext& tube : tubes_) {
      if (tube.LevelExists(level)) {
        tube.AccumulateCoarseFineFluxes(level, time_scale, dir);
      }
    }
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    if (plenum.LevelExists(level)) {
      plenum.AccumulateCoarseFineFluxes(level, time_scale, dir);
    }
  }
}

/// \brief Replace the coarse fluxes by accumulated fine fluxes on the coarse
/// fine interfaces.
void MultiBlockIntegratorContext::ApplyFluxCorrection(int fine, int coarse,
                                                      Duration dt) {
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}}, [=](auto& block) {
    if (block.LevelExists(fine)) {
      block.ApplyFluxCorrection(fine, coarse, dt);
    }
  });
}

/// \brief Resets all accumulates fluxes to zero.
void MultiBlockIntegratorContext::ResetCoarseFineFluxes(int fine, int coarse) {
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}}, [=](auto& block) {
    if (block.LevelExists(fine)) {
      block.ResetCoarseFineFluxes(fine, coarse);
    }
  });
}

/// \brief Coarsen scratch data from a fine level number to a coarse level
/// number.
void MultiBlockIntegratorContext::CoarsenConservatively(int fine, int coarse) {
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}},
               [fine, coarse](auto& block) {
                 if (block.LevelExists(fine)) {
                   block.CoarsenConservatively(fine, coarse);
                 }
               });
}
///@}

} // namespace fub::amrex
