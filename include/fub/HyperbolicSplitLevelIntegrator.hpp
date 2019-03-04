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

#include "fub/Box.hpp"
#include "fub/CartesianCoordinates.hpp"
#include "fub/CompleteFromCons.hpp"
#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/boundary_condition/NoBoundary.hpp"

#include <mpi.h>

namespace fub {
template <typename IntegratorContext, typename PatchIntegrator,
          typename FluxMethod, typename BoundaryCondition = NoBoundary>
class HyperbolicSplitLevelIntegrator : private IntegratorContext {
public:
  using Equation = typename PatchIntegrator::Equation;
  using Context = IntegratorContext;
  using PatchHandle = typename Context::PatchHandle;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;

  static constexpr int Rank() { return Equation::Rank(); }

  HyperbolicSplitLevelIntegrator(
      IntegratorContext context, PatchIntegrator integrator,
      FluxMethod flux_method,
      BoundaryCondition boundary_condition = BoundaryCondition())
      : Context(std::move(context)), patch_integrator_{std::move(integrator)},
        flux_method_{std::move(flux_method)}, boundary_condition_{
                                                  this, boundary_condition} {}

  HyperbolicSplitLevelIntegrator(HyperbolicSplitLevelIntegrator&&) = default;

  const IntegratorContext& GetIntegratorContext() const noexcept {
    return *this;
  }

  View<Complete> GetData(const PatchHandle& patch) {
    return MakeView(*this, boost::hana::type_c<View<Complete>>,
                    Context::GetData(patch), GetEquation());
  }

  View<Complete> GetScratch(const PatchHandle& patch, Direction dir) {
    return MakeView(*this, boost::hana::type_c<View<Complete>>,
                    Context::GetScratch(patch, dir), GetEquation());
  }

  View<Conservative> GetFluxes(const PatchHandle& patch, Direction dir) {
    return MakeView(*this, boost::hana::type_c<View<Conservative>>,
                    Context::GetFluxes(patch, dir), GetEquation());
  }

  const Equation& GetEquation() const noexcept {
    return patch_integrator_.GetEquation();
  }

  using Context::GetCycles;
  using Context::GetPatchHierarchy;
  using Context::GetTimePoint;
  using Context::ResetHierarchyConfiguration;

  void FillGhostLayerSingleLevel(int level, Direction direction) {
    Context::FillGhostLayerSingleLevel(level, direction, boundary_condition_);
  }

  void FillGhostLayerTwoLevels(int fine, int coarse, Direction direction) {
    Context::FillGhostLayerTwoLevels(fine, coarse, direction,
                                     boundary_condition_);
  }

  /// Returns a stable dt across all levels and in one spatial direction.
  ///
  /// This function takes the refinement level into account.
  /// For stability it is advised to use some additional CFL condition.
  ///
  /// \param[in] dir  The direction into which the time step size is calculated.
  Duration ComputeStableDt(Direction dir) {
    int refine_ratio = 1;
    double coarse_dt = std::numeric_limits<double>::infinity();
    for (int level_num = 0; Context::LevelExists(level_num); ++level_num) {
      double level_dt = std::numeric_limits<double>::infinity();
      if (level_num > 0) {
        FillGhostLayerTwoLevels(level_num, level_num - 1, dir);
      } else {
        FillGhostLayerSingleLevel(level_num, dir);
      }
      Context::ForEachPatch(level_num, [&](const PatchHandle& patch) {
        View<const Complete> scratch = GetScratch(patch, dir);
        const double dx = Context::GetDx(patch, dir);
        const double patch_dt = flux_method_.ComputeStableDt(scratch, dx, dir);
        level_dt = std::min(level_dt, patch_dt);
      });
      refine_ratio *= Context::GetRatioToCoarserLevel(level_num);
      coarse_dt = std::min(refine_ratio * level_dt, coarse_dt);
    }
    const MPI_Comm comm = Context::GetMpiCommunicator();
    double global_min_dt{0};
    MPI_Allreduce(&coarse_dt, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, comm);
    return Duration(global_min_dt);
  }

  /// Advance a specified patch level and all finer levels in the specified
  /// split direction by time dt.
  ///
  /// This method subcycles finer levels.
  ///
  /// \param[in] level_num  An integer denoting the patch level where 0 is the
  /// coarsest level.
  ///
  /// \param[in] direction The dimensional split direction which will be used to
  /// advance.
  ///
  /// \param[in] dt A stable time step size for the level_num-th patch level.
  void AdvanceLevel(int this_level, Direction direction, Duration dt,
                    int subcycle = 0) {
    // PreAdvanceLevel might regrid this and all finer levels.
    // The Context must make sure that scratch data is allocated
    Context::PreAdvanceLevel(this_level, direction, dt, subcycle);

    // Fill the ghost layer which is needed for this split direction
    if (subcycle == 0 && this_level > 0) {
      FillGhostLayerTwoLevels(this_level, this_level - 1, direction);
    } else {
      FillGhostLayerSingleLevel(this_level, direction);
    }

    // Compute fluxes in the specified direction
    Context::ForEachPatch(this_level, [&](const PatchHandle& patch) {
      const double dx = Context::GetDx(patch, direction);
      View<const Complete> scratch = GetScratch(patch, direction);
      View<Conservative> fluxes = GetFluxes(patch, direction);
      flux_method_.ComputeNumericFluxes(fluxes, scratch, dt, dx, direction);
    });

    // Accumulate Fluxes on the coarse-fine interface to the next coarser level.
    if (this_level > 0) {
      Context::AccumulateCoarseFineFluxes(this_level, direction, dt);
    }

    // If a finer level exists in the hierarchy, we subcycle that finer level
    // multiple times and use the fine fluxes on coarse-fine interfaces
    const int next_level = this_level + 1;
    if (Context::LevelExists(next_level)) {
      Context::ResetCoarseFineFluxes(next_level, this_level, direction);
      const int refine_ratio = Context::GetRatioToCoarserLevel(next_level);
      for (int r = 0; r < refine_ratio; ++r) {
        AdvanceLevel(next_level, direction, dt / refine_ratio, r);
      }
      Context::ApplyFluxCorrection(next_level, this_level, dt, direction);
    }

    // Use the updated fluxes to update cons variables at the "SCRATCH" context.
    Context::ForEachPatch(this_level, [&](const PatchHandle& patch) {
      View<Conservative> scratch = GetScratch(patch, direction);
      const int gcw = Context::GetGhostCellWidth(patch, direction);
      const double dx = Context::GetDx(patch, direction);
      View<const Conservative> fluxes = GetFluxes(patch, direction);
      StridedView<Conservative> inner =
          ViewInnerRegion(scratch, direction, gcw - 1);
      patch_integrator_.UpdateConservatively(inner, fluxes, inner, dt, dx,
                                             direction);
    });

    // Coarsen inner regions from next finer level to this level.
    if (Context::LevelExists(next_level)) {
      Context::CoarsenConservatively(next_level, this_level, direction);
    }

    // The conservative update and and the coarsening happened on conservative
    // variables. We have to reconstruct the missing variables in the complete
    // state.
    Context::ForEachPatch(this_level, [&](const PatchHandle& patch) {
      View<Complete> state = GetData(patch);
      View<Conservative> scratch = GetScratch(patch, direction);
      const int gcw = Context::GetGhostCellWidth(patch, direction);
      StridedView<const Conservative> inner =
          ViewInnerRegion(scratch, direction, gcw);
      CompleteFromCons(GetEquation(), state, inner);
    });

    Context::PostAdvanceLevel(this_level, direction, dt, subcycle);
  }

private:
  struct AdaptedBoundaryCondition {
    static constexpr int Rank = Equation::Rank();

    void operator()(PatchHandle patch, Location boundary, Duration time_point) {
      Direction dir = boundary.direction;
      View<Complete> scratch = integrator_->GetScratch(patch, dir);
      const CartesianCoordinates& coordinates =
          integrator_->GetCartesianCoordinates(patch);
      const int gcw = integrator_->GetGhostCellWidth(patch, dir);
      const Box<Rank> fill_box =
          GetBoundaryBox(Extents<0>(scratch), boundary, gcw);
      boundary_condition_.FillBoundary(scratch, fill_box, coordinates, boundary,
                                       time_point);
    }

    HyperbolicSplitLevelIntegrator* integrator_;
    BoundaryCondition boundary_condition_;
  };
  PatchIntegrator patch_integrator_;
  FluxMethod flux_method_;
  AdaptedBoundaryCondition boundary_condition_;
};

} // namespace fub