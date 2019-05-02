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

#include "fub/CartesianCoordinates.hpp"
#include "fub/Duration.hpp"
#include "fub/core/function_ref.hpp"
#include "fub/core/span.hpp"
#include "fub/grid/SAMRAI/GriddingAlgorithm.hpp"
#include "fub/grid/SAMRAI/RegisterVariables.hpp"

#include <SAMRAI/geom/CartesianGridGeometry.h>
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

#include <mpi.h>

#include <array>
#include <optional>
#include <vector>

namespace fub {
namespace samrai {

/// This integrator context delegates grid related tasks to the SAMRAI library.
class HyperbolicSplitIntegratorContext {
public:
  static constexpr int MaxVariables = 8;
  using PatchHandle = SAMRAI::hier::Patch*;

  template <typename T>
  using vector = boost::container::static_vector<T, MaxVariables>;

  using BoundaryCondition = function_ref<void(PatchHandle, Location, Duration)>;

  HyperbolicSplitIntegratorContext(const GriddingAlgorithm& gridding,
                                   DataDescription description, int gcw);

  /// Invoke a specified function for each patch of a patch level.
  template <typename F> void ForEachPatch(int level_num, F function) {
    const SAMRAI::hier::PatchLevel& level =
        *GetPatchHierarchy()->getPatchLevel(level_num);
    for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
      function(patch.get());
    }
  }

  /// \brief Returns a native handle to the patch hierarchy.
  const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& GetPatchHierarchy() const
      noexcept;

  /// \brief Regrid level_num and all finer levels if it is the first subcycle.
  void PreAdvanceLevel(int level_num, Direction dir, Duration dt, int subcycle);

  /// Do nothing here.
  void PostAdvanceLevel(int level_num, Direction dir, Duration dt,
                        int subcycle);

  /// \brief Returns the geometric cell width of `patch` in direction `dir`.
  double GetDx(PatchHandle patch, Direction dir) const;

  /// \brief Returns a Coordinates object associated with the specified patch.
  CartesianCoordinates GetCartesianCoordinates(PatchHandle patch) const;

  /// \brief Returns the MPI communicator which is associated with this context.
  MPI_Comm GetMpiCommunicator() const noexcept;

  bool LevelExists(int level) const noexcept;

  int GetRatioToCoarserLevel(int level) const noexcept;

  int GetGhostCellWidth(PatchHandle, Direction dir);

  void FillGhostLayerTwoLevels(int fine, int coarse, Direction direction,
                               BoundaryCondition boundary);

  void FillGhostLayerSingleLevel(int level, Direction direction,
                                 BoundaryCondition boundary);

  void AccumulateCoarseFineFluxes(int level, Direction dir, Duration dt);
  void ApplyFluxCorrection(int fine, int coarse, Duration dt, Direction dir);
  void ResetCoarseFineFluxes(int fine, int coarse, Direction dir);
  void CoarsenConservatively(int fine, int coarse, Direction dir);

  vector<SAMRAI::pdat::CellData<double>*> GetData(PatchHandle patch);

  vector<SAMRAI::pdat::CellData<double>*> GetScratch(PatchHandle patch,
                                                     Direction dir);

  vector<SAMRAI::pdat::SideData<double>*> GetFluxes(PatchHandle patch,
                                                    Direction dir);

  vector<SAMRAI::pdat::OutersideData<double>*>
  GetOutersideFluxes(PatchHandle patch) const;

  Duration GetTimePoint(int level, Direction dir) const;
  void SetTimePoint(Duration dt, int level, Direction dir);

  std::ptrdiff_t GetCycles(int level, Direction dir) const;
  void SetCycles(std::ptrdiff_t cycle, int level, Direction dir);

  const SAMRAI::geom::CartesianGridGeometry& GetGeometry(int level) const;

  double ComputeStableDt(int level, Direction direction);

  void ResetHierarchyConfiguration(int level = 0);

  void
  AccumulateFluxesOnOuterside(span<SAMRAI::pdat::OutersideData<double>*> states,
                              span<const SAMRAI::pdat::SideData<double>*> cons,
                              const SAMRAI::hier::Patch& patch, int direction,
                              double dt);

  void ClearOutersideFluxes(int level);

  GriddingAlgorithm gridding_;
  DataDescription description_;
  struct InternalDataIds {
    std::vector<int> intermediate;
    /// Scratch space for each direction with ghost cells
    std::array<std::vector<int>, SAMRAI_MAXIMUM_DIMENSION> scratch;
    /// SideData for each direction
    std::array<std::vector<int>, SAMRAI_MAXIMUM_DIMENSION> fluxes;
    /// OutersideData for each direction
    std::vector<int> outerside_fluxes;
  };
  InternalDataIds data_ids_;
  std::array<std::vector<Duration>, 3> time_points_{};
  std::array<std::vector<Duration>, 3> regrid_time_points_{};
  std::array<std::vector<std::ptrdiff_t>, 3> cycles_{};
  SAMRAI::hier::IntVector ghost_cell_width_;

  /// The CoarsenAlgorithm describes which patch data ids get coarsened.
  /// This does not need to be rebuild if the box level of the hierarchy
  /// changes.
  std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm> coarsen_outerside_algorithm_;
  std::array<std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>,
             SAMRAI_MAXIMUM_DIMENSION>
      coarsen_inner_region_algorithm_;

  /// The schedules execute the algorithm for concrete patch levels.
  std::vector<std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>>
      coarsen_outerside_{};
  std::array<std::vector<std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>>,
             SAMRAI_MAXIMUM_DIMENSION>
      coarsen_inner_region_{};

  /// The RefineAlgorithm describes which patch data ids get communicated to
  /// The algorithm is independent of the specific hierarchy.
  std::array<std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>,
             SAMRAI_MAXIMUM_DIMENSION>
      fill_ghost_layer_algorithm_{};

  /// The RefineSchedule executes an algorithm on specific patch levels.
  /// This vector stores a schedule for each patch level in the hierarchy.
  /// This needs to be rebuilt whenever the hierarchy changes.
  std::array<std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>, 3>
      fill_ghost_two_levels_{};
  std::array<std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>, 3>
      fill_ghost_single_level_{};
};

} // namespace samrai
} // namespace fub

#endif
