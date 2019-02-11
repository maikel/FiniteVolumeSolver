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

#ifndef FUB_AMREX_PATCH_HIERARCHY_HPP
#define FUB_AMREX_PATCH_HIERARCHY_HPP

#include "fub/Equation.hpp"
#include "fub/PatchLevel.hpp"

#include <AMReX_FluxRegister.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace fub {
namespace amrex {

CartesianCoordinates GetCartesianCoordinates(const ::amrex::Geometry& geom,
                                             const ::amrex::Box& box) {
  Eigen::Vector3d lower{};
  geom.LoNode(box.smallEnd(), lower.data());
  Eigen::Vector3d upper{};
  geom.HiNode(box.bigEnd(), upper.data());
  std::array<std::ptrdiff_t, 3> extents{1, 1, 1};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[i] = box.length(i);
  }
  return CartesianCoordinates(lower, upper, DynamicExtents<3>(extents));
}

struct PatchHierarchyOptions {
  int n_state_components{1};
  int n_flux_components{1};
  int ghost_layer_width{1};
  int max_refinement_level{0};
};

class PatchHierarchy {
public:
  explicit PatchHierarchy(const CartesianGridGeometry& geometry,
                          const PatchHierarchyOptions& options);

  const ::amrex::Geometry& GetGeometry(int level) const noexcept {
    return patch_level_geometry_[level];
  }
  const PatchHierarchyOptions& GetOptions() const noexcept { return options_; }

  PatchLevel& GetPatchLevel(int level) noexcept {
    return patch_level_[level];
  }

  const PatchLevel& GetPatchLevel(int level) const noexcept {
    return patch_level_[level];
  }

private:
  PatchHierarchyOptions options_;
  std::vector<PatchLevel> patch_level_;
  std::vector<::amrex::Geometry> patch_level_geometry_;
};

PatchLevel::PatchLevel(
    const Equation& equation, int level, double tp,
    const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping,
    int num_ghosts_cells)
    : time_point{tp} {
  constexpr int no_ghost_cells = 0;
  // Get number of components for a complete state.
  // Allocate complete state data for the specified patch level.
  auto state_sizes = GetSizes<State>(equation);
  ForEachMember(
      [&](::amrex::MultiFab& data, ::amrex::MultiFab& scratch,
          int n_components) {
        // Our base data has no ghost cell width
        data.define(box_array, distribution_mapping, n_components,
                    no_ghost_cells);

        // We will fill ghost data on scratch
        scratch.define(box_array, distribution_mapping, n_components,
                       num_ghosts_cells);
      },
      data, scratch, state_sizes);

  // Get number of components for a conservative state.
  // Allocate flux buffers for each direction and for the coarse fine interface
  // regions.
  const auto cons_sizes = GetSizes<Conservative>(equation);
  ForEachMember(
      [&](::amrex::MultiFab& flux,
          std::unique_ptr<::amrex::FluxRegister>& coarse_fine,
          int n_components) {
        flux.define(box_array, distribution_mapping, n_components,
                    no_ghost_cells);
        if (level > 0) {
          const ::amrex::IntVect n_ghosts(AMREX_D_DECL(
              num_ghosts_cells, num_ghosts_cells, num_ghosts_cells));
          coarse_fine = std::make_unique<::amrex::FluxRegister>(
              box_array, distribution_mapping, n_ghosts, level, n_components);
        }
      },
      fluxes, coarse_fine, cons_sizes);
}

PatchHierarchy::PatchHierarchy(const Equation& equation,
                                         const CartesianGridGeometry& geometry,
                                         const PatchHierarchyOptions& options)
    : equation_{equation}, geometry_{geometry}, options_{options},
      patch_level_{}, patch_level_geometry_{} {
  patch_level_.resize(options.max_refinement_level + 1);
  patch_level_geometry_.resize(options.max_refinement_level + 1);
  ::amrex::Box level_box(::amrex::IntVect(),
                         ::amrex::IntVect(geometry.cell_dimensions.data()));
  for (::amrex::Geometry& geom : patch_level_geometry_) {
    geom = ::amrex::Geometry(level_box, &geometry.coordinates);
    level_box.refine(2);
  }
}

void PatchHierarchy<Equation>::MakeLevelFromScratch(
    int level, double time_point, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {
  patch_level_[level] =
      PatchLevel<Equation>(equation_, level, time_point, box_array,
                           distribution_mapping, options_.ghost_layer_width);
}

void PatchHierarchy<Equation>::MakeLevelFromCoarse(
    int lev, double time, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {}

void PatchHierarchy<Equation>::RemakeLevel(
    int lev, double time, const ::amrex::BoxArray& box_array,
    const ::amrex::DistributionMapping& distribution_mapping) {}

void PatchHierarchy<Equation>::ClearLevel(int level) {}

} // namespace amrex
} // namespace fub

#endif