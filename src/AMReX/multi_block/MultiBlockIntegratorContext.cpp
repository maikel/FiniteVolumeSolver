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

template <typename T, typename F> void ForEachBlock(T&& tuple, F&& function) {
  std::apply(
      [&](auto&&... ts) {
        (
            [&](auto&& xs) {
              for (auto&& x : xs) {
                function(x);
              }
            }(ts),
            ...);
      },
      std::forward<T>(tuple));
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
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}}, [=](auto& block) {
    if (block.LevelExists(level)) {
      block.ResetHierarchyConfiguration(level);
    }
  });
}

/// \brief Sets the cycle count for a specific level number and direction.
void MultiBlockIntegratorContext::SetCycles(std::ptrdiff_t cycle, int level) {
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}}, [=](auto& block) {
    if (block.LevelExists(level)) {
      block.SetCycles(cycle, level);
    }
  });
}

/// \brief Sets the time point for a specific level number and direction.
void MultiBlockIntegratorContext::SetTimePoint(Duration t, int level) {
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}}, [=](auto& block) {
    if (block.LevelExists(level)) {
      block.SetTimePoint(t, level);
    }
  });
}
/// @}

/// @{
/// \name Member functions relevant for the level integrator algorithm.

/// \brief On each first subcycle this will regrid the data if neccessary.
void MultiBlockIntegratorContext::PreAdvanceLevel(int level_num, Duration dt,
                                                  int subcycle) {
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}}, [=](auto& block) {
    if (block.LevelExists(level_num)) {
      block.PreAdvanceLevel(level_num, dt, subcycle);
    }
  });
}

/// \brief Increases the internal time stamps and cycle counters for the
/// specified level number and direction.
Result<void, TimeStepTooLarge>
MultiBlockIntegratorContext::PostAdvanceLevel(int level_num, Duration dt,
                                              int subcycle) {
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
  BoundaryCondition* base_condition{nullptr};
  cutcell::BoundaryCondition* base_cc_condition{nullptr};

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration time_point,
                    const GriddingAlgorithm& gridding) const {
    FUB_ASSERT(base_condition != nullptr);
    base_condition->FillBoundary(mf, geom, time_point, gridding);
    MultiBlockBoundary* boundary = boundaries.begin();
    for (const BlockConnection& connection : connectivity) {
      if (id == connection.tube.id) {
        boundary->FillBoundary(mf, geom, time_point, gridding);
      }
      ++boundary;
    }
  }
  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration time_point,
                    const cutcell::GriddingAlgorithm& gridding) const {
    FUB_ASSERT(base_cc_condition != nullptr);
    base_cc_condition->FillBoundary(mf, geom, time_point, gridding);
    MultiBlockBoundary* boundary = boundaries.begin();
    for (const BlockConnection& connection : connectivity) {
      if (id == connection.plenum.id) {
        boundary->FillBoundary(mf, geom, time_point, gridding);
      }
      ++boundary;
    }
  }
};
} // namespace

/// \brief Fills the ghost layer of the scratch data and interpolates in the
/// coarse fine layer.
void MultiBlockIntegratorContext::FillGhostLayerTwoLevels(int fine,
                                                          int coarse) {
  std::size_t id = 0;
  for (IntegratorContext& tube : tubes_) {
    if (tube.LevelExists(fine)) {
      BoundaryCondition& fbc = tube.GetBoundaryCondition(fine);
      BoundaryCondition fwrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(), &fbc, nullptr};
      fwrapped.parent = tube.GetGriddingAlgorithm().get();
      fwrapped.geometry = tube.GetGeometry(fine);
      BoundaryCondition& cbc = tube.GetBoundaryCondition(coarse);
      BoundaryCondition cwrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(), &cbc, nullptr};
      cwrapped.parent = tube.GetGriddingAlgorithm().get();
      cwrapped.geometry = tube.GetGeometry(coarse);
      tube.FillGhostLayerTwoLevels(fine, fwrapped, coarse, cwrapped);
    }
    id += 1;
  }
  id = 0;
  for (cutcell::IntegratorContext& plenum : plena_) {
    if (plenum.LevelExists(fine)) {
      cutcell::BoundaryCondition& fbc = plenum.GetBoundaryCondition(fine);
      cutcell::BoundaryCondition fwrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(), nullptr, &fbc};
      fwrapped.parent = plenum.GetGriddingAlgorithm().get();
      fwrapped.geometry = plenum.GetGeometry(fine);
      cutcell::BoundaryCondition& cbc = plenum.GetBoundaryCondition(coarse);
      cutcell::BoundaryCondition cwrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(), nullptr, &cbc};
      cwrapped.parent = plenum.GetGriddingAlgorithm().get();
      cwrapped.geometry = plenum.GetGeometry(fine);
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
      BoundaryCondition& bc = tube.GetBoundaryCondition(level);
      BoundaryCondition wrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(), &bc, nullptr};
      wrapped.parent = tube.GetGriddingAlgorithm().get();
      wrapped.geometry = tube.GetGeometry(level);
      tube.FillGhostLayerSingleLevel(level, wrapped);
    }
    id += 1;
  }
  id = 0;
  for (cutcell::IntegratorContext& plenum : plena_) {
    if (plenum.LevelExists(level)) {
      cutcell::BoundaryCondition& bc = plenum.GetBoundaryCondition(level);
      cutcell::BoundaryCondition wrapped =
          WrapBoundaryCondition{id, gridding_->GetConnectivity(),
                                gridding_->GetBoundaries(), nullptr, &bc};
      wrapped.parent = plenum.GetGriddingAlgorithm().get();
      wrapped.geometry = plenum.GetGeometry(level);
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
  MultiBlockBoundary* boundary = gridding_->GetBoundaries().begin();
  for (const BlockConnection& connection : gridding_->GetConnectivity()) {
    IntegratorContext& tube = tubes_[connection.tube.id];
    cutcell::IntegratorContext& plenum = plena_[connection.plenum.id];
    boundary->ComputeBoundaryData(plenum.GetPatchHierarchy(),
                                  tube.GetPatchHierarchy());
  }
}

void MultiBlockIntegratorContext::PostAdvanceHierarchy() {
  for (cutcell::IntegratorContext& ctx : plena_) {
    ctx.PostAdvanceHierarchy();
  }
}

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
  ForEachBlock(std::tuple{span{tubes_}, span{plena_}}, [=](auto& block) {
    if (block.LevelExists(fine)) {
      block.CoarsenConservatively(fine, coarse);
    }
  });
}
///@}

} // namespace fub::amrex
