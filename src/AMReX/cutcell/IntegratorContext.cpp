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

#include "fub/AMReX/cutcell/IntegratorContext.hpp"

// Polymorphic Strategies
#include "fub/AMReX/cutcell/BoundaryCondition.hpp"
#include "fub/AMReX/cutcell/FluxMethod.hpp"
#include "fub/AMReX/cutcell/Reconstruction.hpp"
#include "fub/AMReX/cutcell/TimeIntegrator.hpp"

#include "fub/AMReX/cutcell/IndexSpace.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/utility.hpp"

#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxReg_C.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>

namespace fub::amrex::cutcell {

////////////////////////////////////////////////////////////////////////////////
//                                                  IntegratorContext::LevelData

////////////////////////////////////////////////////////////////////////////////
//                                                      Move Assignment Operator

IntegratorContext::LevelData& IntegratorContext::LevelData::
operator=(LevelData&& other) noexcept {
  if (other.coarse_fine.fineLevel() > 0) {
    // If we do not invoke clear in beforehand it will throw an error in AMReX
    coarse_fine.clear();
    coarse_fine.define(scratch[0].boxArray(), scratch[0].DistributionMap(),
                       other.coarse_fine.refRatio(),
                       other.coarse_fine.fineLevel(),
                       other.coarse_fine.nComp());
  }
  eb_factory = std::move(other.eb_factory);
  reference_states = std::move(other.reference_states);
  boundary_fluxes = std::move(other.boundary_fluxes);
  scratch = std::move(other.scratch);
  fluxes = std::move(other.fluxes);
  stabilized_fluxes = std::move(other.stabilized_fluxes);
  shielded_left_fluxes = std::move(other.shielded_left_fluxes);
  shielded_right_fluxes = std::move(other.shielded_right_fluxes);
  doubly_shielded_fluxes = std::move(other.doubly_shielded_fluxes);
  time_point = other.time_point;
  regrid_time_point = other.regrid_time_point;
  cycles = other.cycles;
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
//                                                             IntegratorContext

////////////////////////////////////////////////////////////////////////////////
//                                          Constructor and Assignment Operators

IntegratorContext::IntegratorContext(
    std::shared_ptr<GriddingAlgorithm> gridding, HyperbolicMethod nm)
    : ghost_cell_width_{nm.flux_method.GetStencilWidth() + 1},
      gridding_{std::move(gridding)}, data_{}, method_{
                                                   std::move(nm)} {
  data_.reserve(
      static_cast<std::size_t>(GetPatchHierarchy().GetMaxNumberOfLevels()));
  // Allocate auxiliary data arrays for each refinement level in the hierarchy
  ResetHierarchyConfiguration();
}

IntegratorContext::IntegratorContext(const IntegratorContext& other)
    : ghost_cell_width_{other.ghost_cell_width_}, gridding_{other.gridding_},
      data_(static_cast<std::size_t>(GetPatchHierarchy().GetNumberOfLevels())),
      method_{other.method_} {
  // Allocate auxiliary data arrays
  ResetHierarchyConfiguration();
  // Copy time stamps and cycle counters
  std::size_t n_levels =
      static_cast<std::size_t>(GetPatchHierarchy().GetNumberOfLevels());
  for (std::size_t i = 0; i < n_levels; ++i) {
    data_[i].cycles = other.data_[i].cycles;
    data_[i].time_point = other.data_[i].time_point;
    data_[i].regrid_time_point = other.data_[i].regrid_time_point;
  }
}

IntegratorContext& IntegratorContext::IntegratorContext::
operator=(const IntegratorContext& other) {
  // We use the copy and move idiom to provide the strong exception guarantee.
  // If an exception occurs we do not change the original object.
  IntegratorContext tmp{other};
  return (*this = std::move(tmp));
}

///////////////////////////////////////////////////////////////////////////////
//                                                             Member Accessors



}