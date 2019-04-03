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

#include "fub/CartesianCoordinates.hpp"
#include "fub/equations/PerfectGas.hpp"

#include "fub/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/HyperbolicSplitSystemSolver.hpp"
#include "fub/ext/Eigen.hpp"
#include "fub/flux_method/HllMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"
#include "fub/grid/AMReX/FluxMethod.hpp"
#include "fub/grid/AMReX/GriddingAlgorithm.hpp"
#include "fub/grid/AMReX/HyperbolicSplitIntegratorContext.hpp"
#include "fub/grid/AMReX/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/grid/AMReX/Reconstruction.hpp"
#include "fub/grid/AMReX/ScopeGuard.hpp"
#include "fub/split_method/StrangSplitting.hpp"
#include "fub/tagging/GradientDetector.hpp"
#include "fub/tagging/TagBuffer.hpp"

#include "fub/RunSimulation.hpp"

#include <fmt/format.h>
#include <iostream>

struct ShockData {
  using Equation = fub::PerfectGas<2>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  void InitializeData(const fub::View<Complete>& data,
                      fub::amrex::PatchHandle patch) {
    using namespace fub;
    const ::amrex::Geometry& geom = hierarchy->GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->tilebox();
    CartesianCoordinates x = fub::amrex::GetCartesianCoordinates(geom, box);

    ForEachIndex(Box<0>(data), [&](auto... is) {
      if (x(is...)[0] < -0.04) {
        Store(data, left, {is...});
      } else {
        Store(data, right, {is...});
      }
    });
  }

  std::shared_ptr<fub::amrex::PatchHierarchy> hierarchy;
  Equation equation;
  Complete left;
  Complete right;
};


struct TransmissiveBoundary {
  using Equation = fub::PerfectGas<2>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  void operator()(const fub::PatchDataView<double, 3>& data,
                  fub::amrex::PatchHandle, fub::Location location,
                  int fill_width, fub::Duration) {
    fub::View<Complete> complete =
        fub::amrex::MakeView<fub::View<Complete>>(data, equation);
    const std::size_t dir = static_cast<std::size_t>(location.direction);
    const std::size_t other_dir = (dir + 1) % 2;
    std::array<std::ptrdiff_t, 3> origin = data.Origin();
    if (location.side == 0) {
      for (std::ptrdiff_t j = 0; j < data.Extent(other_dir); ++j) {
        std::array<std::ptrdiff_t, 2> dest_index{origin[0], origin[1]};
        std::array<std::ptrdiff_t, 2> source_index{origin[0], origin[1]};
        dest_index[other_dir] += j;
        source_index[other_dir] += j;
        source_index[dir] += fill_width;
        Load(state, complete, source_index);
        for (int i = 0; i < fill_width; ++i) {
          dest_index[dir] = origin[dir] + i;
          Store(complete, state, dest_index);
        }
      }
    } else if (location.side == 1) {
      const std::ptrdiff_t n = fub::Extents<0>(complete).extent(dir);
      for (std::ptrdiff_t j = 0; j < data.Extent(other_dir); ++j) {
        std::array<std::ptrdiff_t, 2> dest_index{origin[0], origin[1]};
        std::array<std::ptrdiff_t, 2> source_index{origin[0], origin[1]};
        dest_index[other_dir] += j;
        source_index[other_dir] += j;
        source_index[dir] += n - 1 - fill_width;
        Load(state, complete, source_index);
        for (int i = 0; i < fill_width; ++i) {
          dest_index[dir] = origin[dir] + n - fill_width + i;
          Store(complete, state, dest_index);
        }
      }
    }
  }

  std::shared_ptr<fub::amrex::PatchHierarchy> hierarchy;
  Equation equation;
  Complete state{equation};
};

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::amrex::ScopeGuard guard(argc, argv);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM == 2);

  const std::array<int, Dim> n_cells{128, 128};

  const std::array<double, Dim> xlower{-1.0, -1.0};
  const std::array<double, Dim> xupper{+1.0, +1.0};

  fub::PerfectGas<2> equation{};

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = std::array<int, 2>{0, 0};

  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 1;

  auto hierarchy =
      std::make_shared<fub::amrex::PatchHierarchy>(desc, geometry, hier_opts);

  using Complete = fub::PerfectGas<2>::Complete;
  fub::GradientDetector gradient{std::make_pair(&Complete::density, 0.001),
                                 std::make_pair(&Complete::pressure, 0.1)};
  fub::TagBuffer buffer{4};

  auto from_prim = [](Complete& state, const fub::PerfectGas<2>& equation) {
    state.energy = state.pressure * equation.gamma_minus_1_inv;
    state.speed_of_sound =
        std::sqrt(equation.gamma * state.pressure / state.density);
  };

  Complete left;
  left.density = 1.0;
  left.momentum = 0.0;
  left.pressure = 1000.0;
  from_prim(left, equation);

  Complete right;
  right.density = 1.0;
  right.momentum = 0.0;
  right.pressure = 0.01;
  from_prim(right, equation);

  ShockData initial_data{hierarchy, equation, left, right};
  TransmissiveBoundary boundary{hierarchy, equation};
  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      hierarchy, fub::amrex::AdaptInitialData(initial_data, equation),
      fub::amrex::AdaptTagging(equation, hierarchy, gradient, buffer), boundary);
  gridding->InitializeHierarchy(0.0);

  fub::HyperbolicSplitPatchIntegrator patch_integrator{equation};
  fub::MusclHancockMethod flux_method{equation};

  const int gcw = flux_method.GetStencilWidth();
  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      fub::amrex::HyperbolicSplitIntegratorContext(std::move(gridding), gcw),
      fub::amrex::HyperbolicSplitPatchIntegrator(patch_integrator),
      fub::amrex::FluxMethod(flux_method),
      fub::amrex::Reconstruction(equation)));

  std::string base_name = "AMReX/PerfectGas2d_";

  auto output = [&](auto& hierarchy, int cycle, fub::Duration) {
    std::string name = fmt::format("{}{:04}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::WritePlotFile(name, *hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  auto print_msg = [](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(hierarchy, 0, 0.0s);
  fub::RunOptions run_options{};
  
  run_options.final_time = 1.0s;
  run_options.output_frequency = 1;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}
