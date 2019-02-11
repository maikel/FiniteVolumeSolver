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

#ifndef FUB_SAMRAI_HYPERBOLIC_SPLIT_LEVEL_INTEGRATOR_HPP
#define FUB_SAMRAI_HYPERBOLIC_SPLIT_LEVEL_INTEGRATOR_HPP

#include "fub/core/span.hpp"

#include "fub/SAMRAI/DimensionalSplitBoundaryCondition.hpp"
#include "fub/SAMRAI/DimensionalSplitFluxMethod.hpp"
#include "fub/SAMRAI/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/SAMRAI/RegisterVariables.hpp"

#include <SAMRAI/hier/CoarseFineBoundary.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/OutersideData.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/CoarsenSchedule.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>

#include <boost/container/static_vector.hpp>

#include <array>
#include <vector>

namespace fub {
namespace samrai {
struct BoundaryCondition : public SAMRAI::xfer::RefinePatchStrategy {
  BoundaryCondition(DimensionalSplitBoundaryCondition* base) : base_{base} {}

  void setPhysicalBoundaryConditions(
      SAMRAI::hier::Patch& patch, double fill_time,
      const SAMRAI::hier::IntVector& ghost_width_to_fill) override;

  SAMRAI::hier::IntVector
  getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const override {
    return SAMRAI::hier::IntVector::getZero(dim);
  }

  void preprocessRefine(SAMRAI::hier::Patch&, const SAMRAI::hier::Patch&,
                        const SAMRAI::hier::Box&,
                        const SAMRAI::hier::IntVector&) override {}

  void postprocessRefine(SAMRAI::hier::Patch&, const SAMRAI::hier::Patch&,
                         const SAMRAI::hier::Box&,
                         const SAMRAI::hier::IntVector&) override {}

  DimensionalSplitBoundaryCondition* base_;
};

class HyperbolicSplitLevelIntegrator {
public:
  struct InternalDataIds {
    /// intermediate destination space without ghost cells
    std::vector<int> intermediate;
    /// Scratch space for each direction with ghost cells
    std::array<std::vector<int>, 3> scratch;
    /// SideData for each direction
    std::array<std::vector<int>, 3> fluxes;
    /// OutersideData for each direction
    std::vector<int> outerside_fluxes;
  };

  static constexpr int MaxVariables = 8;

  HyperbolicSplitLevelIntegrator(
      PatchDataIdSet description,
      std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
      HyperbolicSplitPatchIntegrator& patch_integrator,
      DimensionalSplitFluxMethod& flux_method,
      DimensionalSplitBoundaryCondition& boundary_condition);

  double ComputeStableDt(int direction);

  void AdvanceLevel(const SAMRAI::hier::PatchLevel& level, int direction,
                    double dt);

  void FillGhostLayer(const SAMRAI::hier::PatchLevel& level, int direction);

  void ResetHierarchyConfiguration(
      std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy);

private:
  void
  AccumulateFluxesOnOuterside(span<SAMRAI::pdat::OutersideData<double>*> states,
                              span<const SAMRAI::pdat::SideData<double>*> cons,
                              const SAMRAI::hier::Patch& patch, int direction,
                              double dt);

  void ClearOutersideFluxes(const SAMRAI::hier::PatchLevel& level);

  void ConservativelyCoarsenOutersideFluxes(
      const SAMRAI::hier::PatchLevel& next_level,
      const SAMRAI::hier::PatchLevel& coarse_level, int direction);

  void ConservativelyCoarsenInnerRegions(
      const SAMRAI::hier::PatchLevel& fine_level,
      const SAMRAI::hier::PatchLevel& coarse_level, int direction);

  void ApplyCoarsenedFluxesOnLevel(const SAMRAI::hier::PatchLevel& fine_level,
                                   const SAMRAI::hier::PatchLevel& coarse_level,
                                   double dt, int diretion);

  template <typename T>
  using StaticVector = boost::container::static_vector<T, MaxVariables>;

  StaticVector<SAMRAI::pdat::CellData<double>*>
  GetScratch(const SAMRAI::hier::Patch& patch, int direction) const;

  StaticVector<SAMRAI::pdat::CellData<double>*>
  GetScratchCons(const SAMRAI::hier::Patch& patch, int direction) const;

  StaticVector<SAMRAI::pdat::CellData<double>*>
  GetNextState(const SAMRAI::hier::Patch& patch) const;

  StaticVector<SAMRAI::pdat::CellData<double>*>
  GetNextCons(const SAMRAI::hier::Patch& patch) const;

  StaticVector<SAMRAI::pdat::SideData<double>*>
  GetFluxes(const SAMRAI::hier::Patch& patch, int direction) const;

  StaticVector<SAMRAI::pdat::OutersideData<double>*>
  GetOutersideFluxes(const SAMRAI::hier::Patch& patch) const;

  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
  PatchDataIdSet description_;
  HyperbolicSplitPatchIntegrator* patch_integrator_;
  DimensionalSplitFluxMethod* flux_method_;
  BoundaryCondition boundary_condition_;

  InternalDataIds data_ids_;

  std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm> coarsen_outerside_algorithm_;
  std::vector<std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>>
      coarsen_outerside_{};

  std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
      coarsen_inner_region_algorithm_;
  std::vector<std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>>
      coarsen_inner_region_{};

  /// The RefineAlgorithm describes which patch data ids get communicated to
  /// The algorithm is independent of the hierarchy.
  std::array<std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>, 3>
      fill_ghost_layer_algorithm_{};
  /// The RefineSchedule executes an algorithm on specific patch levels.
  /// This vector stores a schedule for each patch level in the hierarchy.
  /// This needs to be rebuilt whenever the hierarchy changes.
  std::array<std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>, 3>
      fill_ghost_layer_{};
};

} // namespace samrai
} // namespace fub

#endif