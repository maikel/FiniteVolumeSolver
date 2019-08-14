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
namespace {
Direction LastDirection() noexcept { return Direction{AMREX_SPACEDIM - 1}; }
} // namespace

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
  ForEachBlock(std::tuple{tubes_, plena_},
               [=](auto& block) { block.ResetHierarchyConfiguration(level); });
}

/// \brief Sets the cycle count for a specific level number and direction.
void MultiBlockIntegratorContext::SetCycles(std::ptrdiff_t cycle, int level) {
    for (IntegratorContext& tube : tubes_) {
      tube.SetCycles(cycle, level);
    }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.SetCycles(cycle, level);
  }
}

/// \brief Sets the time point for a specific level number and direction.
void MultiBlockIntegratorContext::SetTimePoint(Duration t, int level) {
  for (IntegratorContext& tube : tubes_) {
    tube.SetTimePoint(t, level);
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.SetTimePoint(t, level);
  }
}
/// @}

/// @{
/// \name Member functions relevant for the level integrator algorithm.

/// \brief On each first subcycle this will regrid the data if neccessary.
void MultiBlockIntegratorContext::PreAdvanceLevel(int level_num,
                                                  Duration dt, int subcycle) {
  for (IntegratorContext& tube : tubes_) {
    tube.PreAdvanceLevel(level_num, dt, subcycle);
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.PreAdvanceLevel(level_num, dt, subcycle);
  }
}

/// \brief Increases the internal time stamps and cycle counters for the
/// specified level number and direction.
Result<void, TimeStepTooLarge>
MultiBlockIntegratorContext::PostAdvanceLevel(int level_num,
                                              Duration dt, int subcycle) {
  for (IntegratorContext& tube : tubes_) {
        Result<void, TimeStepTooLarge> result = tube.PostAdvanceLevel(level_num, dt, subcycle);
	if (!result) {
		return result;
	}
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
      Result<void, TimeStepTooLarge> result = plenum.PostAdvanceLevel(level_num, dt, subcycle);
	if (!result) {
	   return result;
        }
  }
  return boost::outcome_v2::success();
}

/// \brief Fills the ghost layer of the scratch data and interpolates in the
/// coarse fine layer.
void MultiBlockIntegratorContext::FillGhostLayerTwoLevels(int level, int coarse) {
  for (IntegratorContext& tube : tubes_) {
    tube.FillGhostLayerTwoLevels(level, coarse);
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.FillGhostLayerTwoLevels(level, coarse);
  }
}

namespace {
struct WrapBoundaryCondition {
  MultiBlockBoundary* boundary;
  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration time_point,
                    const GriddingAlgorithm& gridding) const {
    boundary->FillBoundary(mf, geom, time_point, gridding);
  }
  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration time_point,
                    const cutcell::GriddingAlgorithm& gridding) const {
    boundary->FillBoundary(mf, geom, time_point, gridding);
  }
};
} // namespace

/// \brief Fills the ghost layer of the scratch data and does nothing in the
/// coarse fine layer.
void MultiBlockIntegratorContext::FillGhostLayerSingleLevel(
    int level) {
  {
    for (IntegratorContext& tube : tubes_) {
      tube.FillGhostLayerSingleLevel(level);
    }
    MultiBlockBoundary* boundary = gridding_->GetBoundaries().begin();
    for (const BlockConnection& connection : gridding_->GetConnectivity()) {
        IntegratorContext& tube = tubes_[connection.tube.id];
        BoundaryCondition bc(WrapBoundaryCondition{boundary});
        bc.parent = tube.GetGriddingAlgorithm().get();
        bc.geometry = tube.GetGeometry(level);
        tube.FillGhostLayerSingleLevel(level, bc);
       ++boundary;
    }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.FillGhostLayerSingleLevel(level);
  }
  }
  {
  MultiBlockBoundary* boundary = gridding_->GetBoundaries().begin();
  for (const BlockConnection& connection : gridding_->GetConnectivity()) {
      cutcell::IntegratorContext& plenum = plena_[connection.plenum.id];
      cutcell::BoundaryCondition bc(WrapBoundaryCondition{boundary});
      bc.parent = plenum.GetGriddingAlgorithm().get();
      bc.geometry = plenum.GetGeometry(level);
      plenum.FillGhostLayerSingleLevel(level, bc);
      ++boundary;
  }
  }
}

/// \brief Returns a estimate for a stable time step size which can be taken
/// for specified level number in direction dir.
Duration MultiBlockIntegratorContext::ComputeStableDt(int level,
                                                      Direction dir) {
  Duration stable_dt{std::numeric_limits<double>::infinity()};
  if (dir == Direction::X) {
    stable_dt =
        std::accumulate(tubes_.begin(), tubes_.end(), stable_dt,
                    [=](Duration dt, IntegratorContext& ctx) {
                      return std::min(dt, ctx.ComputeStableDt(level, dir));
                    });
  }
  stable_dt =
      std::accumulate(plena_.begin(), plena_.end(), stable_dt,
                  [=](Duration dt, cutcell::IntegratorContext& ctx) {
                    return std::min(dt, ctx.ComputeStableDt(level, dir));
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
      tube.ComputeNumericFluxes(level, dt, dir);
    }
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.ComputeNumericFluxes(level, dt, dir);
  }
}

/// \brief Apply a conservative time update for each conservative variable on
/// the specified level number and direction.
void MultiBlockIntegratorContext::UpdateConservatively(int level, Duration dt,
                                                       Direction dir) {
  if (dir == LastDirection()) {
    for (IntegratorContext& tube : tubes_) {
      tube.UpdateConservatively(level, dt, Direction::X);
    }
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.UpdateConservatively(level, dt, dir);
  }
}

/// \brief Reconstruct complete state variables from conservative ones.
void MultiBlockIntegratorContext::CompleteFromCons(int level, Duration dt) {
  for (IntegratorContext& tube : tubes_) {
    tube.CompleteFromCons(level, dt);
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.CompleteFromCons(level, dt);
  }
}

/// \brief Accumulate fluxes on the coarse fine interfaces for a specified
/// fine level number.
void MultiBlockIntegratorContext::AccumulateCoarseFineFluxes(int level,
                                                             double time_scale,
                                                             Direction dir) {
  if (dir == Direction::X) {
    for (IntegratorContext& tube : tubes_) {
      tube.AccumulateCoarseFineFluxes(level, time_scale, dir);
    }
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.AccumulateCoarseFineFluxes(level, time_scale, dir);
  }
}

/// \brief Replace the coarse fluxes by accumulated fine fluxes on the coarse
/// fine interfaces.
void MultiBlockIntegratorContext::ApplyFluxCorrection(int fine, int coarse,
                                                      Duration dt) {
  for (IntegratorContext& tube : tubes_) {
    tube.ApplyFluxCorrection(fine, coarse, dt);
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.ApplyFluxCorrection(fine, coarse, dt);
  }
}

/// \brief Resets all accumulates fluxes to zero.
void MultiBlockIntegratorContext::ResetCoarseFineFluxes(int fine, int coarse) {
    for (IntegratorContext& tube : tubes_) {
      tube.ResetCoarseFineFluxes(fine, coarse);
    }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.ResetCoarseFineFluxes(fine, coarse);
  }
}

/// \brief Coarsen scratch data from a fine level number to a coarse level
/// number.
void MultiBlockIntegratorContext::CoarsenConservatively(int fine, int coarse) {
  for (IntegratorContext& tube : tubes_) {
    tube.CoarsenConservatively(fine, coarse);
  }
  for (cutcell::IntegratorContext& plenum : plena_) {
    plenum.CoarsenConservatively(fine, coarse);
  }
}
///@}

} // namespace fub::amrex
