#ifndef FUB_AMREX_HYPERBOLIC_SPLIT_LEVEL_INTEGRATOR_HPP
#define FUB_AMREX_HYPERBOLIC_SPLIT_LEVEL_INTEGRATOR_HPP

namespace fub {
namespace amrex {

class HyperbolicSplitLevelIntegrator {
  HyperbolicSplitLevelIntegrator(
      PatchHierarchy hierarchy, HyperbolicSplitPatchIntegrator patch_integrator,
      DimensionalSplitFluxMethod flux_method,
      DimensionalSplitBoundaryCondition boundary_condition);

  void AdvanceLevel(PatchLevel& level, int direction, double dt);

  void FillGhostLayer(PatchLevel& level, int direction, double dt);

  void ResetHierarchyConfiguration(PatchHierarchy hierarchy);

private:
  void AddFluxesToOuterface(const Patch& patch, int direction, double dt);

  void CoarsenOuterfaceFluxes(const PatchLevel& next_level,
                              int coarse_level_num, int direction);

  void ConservativeFixUpOnCoarseFineBoundary(const PatchLevel& level,
                                             int diretion);

  void CoarsenInnerRegions(PatchLevel& next_level, int coarse_level_num,
                           int direction);

  PatchHierarchy hierarchy_;
  DataDescription description_;
  HyperbolicSplitPatchIntegrator patch_integrator_;
  DimensionalSplitFluxMethod flux_method_;
  DimensionalSplitBoundaryCondition boundary_condition_;

  struct InternalLevelData {
    /// This corresponds to the "NEW" context which we need to store intermediate
    /// dimensional split steps.
    std::vector<::amrex::MultiFab> new_time_level_;

    /// Scratch space with ghost cell widths
    std::array<std::vector<::amrex::MultiFab>, 3> scratch_;

    /// These arrays will store the fluxes for each patch level which is present
    /// in the patch hierarchy. These will need to be rebuilt if the
    /// PatchHierarchy changes.
    std::array<std::vector<::amrex::MultiFab>, 3> fluxes_;

    /// FluxRegister accumulate fluxes on coarse fine interfaces between
    /// refinement level. These will need to be rebuilt whenever the hierarchy
    /// changes.
    std::vector<::amrex::FluxRegister> coarse_fine_;
  };
  std::vector<InternalLevelData> internal_level_data;
};

} // namespace amrex
} // namespace fub

#endif