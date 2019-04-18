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

#ifndef FUB_HYPERBOLIC_SPLIT_LEVEL_INTEGRATOR_HPP
#define FUB_HYPERBOLIC_SPLIT_LEVEL_INTEGRATOR_HPP

#include "fub/CartesianCoordinates.hpp"
#include "fub/CompleteFromCons.hpp"
#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"

#include <mpi.h>

namespace fub {
template <typename IntegratorContext, typename PatchIntegrator,
          typename FluxMethod, typename Reconstruction>
class HyperbolicSplitLevelIntegrator : private IntegratorContext {
public:
  using Equation = typename PatchIntegrator::Equation;
  using Context = IntegratorContext;
  using PatchHandle = typename Context::PatchHandle;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;

  static constexpr int Rank() { return Equation::Rank(); }

  HyperbolicSplitLevelIntegrator(IntegratorContext context,
                                 PatchIntegrator integrator,
                                 FluxMethod flux_method,
                                 Reconstruction reconstruction)
      : Context(std::move(context)), patch_integrator_{std::move(integrator)},
        flux_method_{std::move(flux_method)}, reconstruction_{
                                                  std::move(reconstruction)} {}

  HyperbolicSplitLevelIntegrator(HyperbolicSplitLevelIntegrator&&) = default;

  IntegratorContext& GetIntegratorContext() noexcept { return *this; }

  FluxMethod& GetFluxMethod() noexcept { return flux_method_; }

  const IntegratorContext& GetIntegratorContext() const noexcept {
    return *this;
  }

  BasicView<Complete> GetData(const PatchHandle& patch) {
    typename Context::template MakeView<BasicView<Complete>> make_view{};
    return make_view(*this, Context::GetData(patch), GetEquation());
  }

  BasicView<Complete> GetScratch(const PatchHandle& patch, Direction dir) {
    typename Context::template MakeView<BasicView<Complete>> make_view{};
    return make_view(*this, Context::GetScratch(patch, dir), GetEquation());
  }

  BasicView<Conservative> GetFluxes(const PatchHandle& patch, Direction dir) {
    typename Context::template MakeView<BasicView<Conservative>> make_view{};
    return make_view(*this, Context::GetFluxes(patch, dir), GetEquation());
  }

  const Equation& GetEquation() const noexcept {
    return patch_integrator_.GetEquation();
  }

  using Context::FillGhostLayerSingleLevel;
  using Context::FillGhostLayerTwoLevels;
  using Context::GetCycles;
  using Context::GetPatchHierarchy;
  using Context::GetTimePoint;
  using Context::ResetHierarchyConfiguration;

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
      if (level_num > 0) {
        FillGhostLayerTwoLevels(level_num, level_num - 1, dir);
      } else {
        FillGhostLayerSingleLevel(level_num, dir);
      }
      const double level_dt = Context::Minimum(
          level_num,
          [dir, context = &GetIntegratorContext(),
           fm = flux_method_](const PatchHandle& patch) mutable -> double {
            return fm.ComputeStableDt(*context, patch, dir);
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
    Context::ForEachPatch(
        this_level, [direction, dt, context = &GetIntegratorContext(),
                     fm = flux_method_](const PatchHandle& patch) mutable {
          fm.ComputeNumericFluxes(*context, patch, direction, dt);
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
    Context::ForEachPatch(
        this_level, [direction, dt, context = &GetIntegratorContext(),
                     pi = patch_integrator_](const PatchHandle& patch) mutable {
          pi.UpdateConservatively(*context, patch, direction, dt);
        });

    // Coarsen inner regions from next finer level to this level.
    if (Context::LevelExists(next_level)) {
      Context::CoarsenConservatively(next_level, this_level, direction);
    }

    Context::ForEachPatch(
        this_level, [direction, dt, context = &GetIntegratorContext(),
                     rec = reconstruction_](const PatchHandle& patch) mutable {
          rec.CompleteFromCons(*context, patch, direction, dt);
        });

    // The conservative update and and the coarsening happened on conservative
    // variables. We have to reconstruct the missing variables in the complete
    // state.
    Context::PostAdvanceLevel(this_level, direction, dt, subcycle);
  }

private:
  PatchIntegrator patch_integrator_;
  FluxMethod flux_method_;
  Reconstruction reconstruction_;
};

} // namespace fub

#endif
