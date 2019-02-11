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

#include "fub/equations/Advection.hpp"
#include "fub/CartesianCoordinates.hpp"

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/ScopeGuard.hpp"
#include "fub/Eigen.hpp"
#include "fub/tagging/GradientDetector.hpp"

#include <AMReX_PlotFileUtil.H>

#include <iostream>

struct WaveData {
  static double wave(double x) {
    return std::exp(-20 * x * x) * std::sin(4.0 * M_PI * x);
  }

  template <typename T>
  void InitializeData(T&& states,
                      const fub::CartesianCoordinates& coords) const {
    fub::ForEachIndex(Mapping(states), [&](auto... is) {
      fub::Advection<2>::State q{{.mass = wave(coords(is...).norm())}};
      states(is...) = q;
    });
  }
};

int main(int argc, char** argv) {
  const fub::amrex::ScopeGuard guard(argc, argv);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM == 2);

  const std::array<int, Dim> n_cells{128, 128};
  const std::array<double, Dim> xlower{-1.0, -1.0};
  const std::array<double, Dim> xupper{+1.0, +1.0};

  const std::array<double, Dim> velocity{1.0, 0.0};
  const fub::Advection equation{velocity};
  const fub::amrex::CartesianGridGeometry geometry{n_cells, {xlower, xupper}};

  fub::amrex::PatchHierarchy hierarchy(
      equation, geometry, {.ghost_layer_width = 1, .max_refinement_level = 1});

  using State = fub::Advection<2>::State;
  fub::GradientDetector tagging(std::pair{&State::mass, 1e-2});
  {
    WaveData initial_data{};
    fub::amrex::GriddingAlgorithm initial_gridding(hierarchy, initial_data,
                                                   tagging);
    initial_gridding.MakeCoarsestLevel(0.0);
  }

  // ...
  // dt = cfl * level_integrator.ComputeStableDt(hierarchy);
  // level_integrator.AdvanceLevel(hierarchy.GetCoarsestLevel(), dt)
  // ...
  // if (regrid_frequency > 0 && cycle % regrid_frequency == 0) {
  //   gridding.Regrid(cycle, time_point);
  /// }

  ::amrex::WriteMultiLevelPlotfile(
      "AMReX_Output", 2,
      {&hierarchy.patch_level(0).data.mass,
       &hierarchy.patch_level(1).data.mass},
      {"mass"}, {hierarchy.geometry(0), hierarchy.geometry(1)}, 0.0,
      {0, 0}, {{1, 1}, {2, 2}});
   ::amrex::WriteMultiLevelPlotfile(
      "AMReX_Output", 2,
      {&hierarchy.patch_level(0).data.mass,
       &hierarchy.patch_level(1).data.mass},
      {"velocity"}, {hierarchy.geometry(0), hierarchy.geometry(1)}, 0.0,
      {0, 0}, {{1, 1}, {2, 2}});
}