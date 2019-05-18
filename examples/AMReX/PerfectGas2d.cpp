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

#include "fub/AMReX/FluxMethod.hpp"
#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/HyperbolicSplitIntegratorContext.hpp"
#include "fub/AMReX/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/AMReX/Reconstruction.hpp"
#include "fub/AMReX/ScopeGuard.hpp"
#include "fub/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/HyperbolicSplitSystemSolver.hpp"
#include "fub/boundary_condition/TransmissiveBoundary.hpp"
#include "fub/ext/Eigen.hpp"
#include "fub/flux_method/HllMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"
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
                      const fub::amrex::PatchHierarchy& hierarchy,
                      fub::amrex::PatchHandle patch) {
    using namespace fub;
    const ::amrex::Geometry& geom = hierarchy.GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->tilebox();
    CartesianCoordinates x = fub::amrex::GetCartesianCoordinates(geom, box);

    ForEachIndex(Box<0>(data), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
      if (x(i, j)[0] < -0.04) {
        Store(data, left, {i, j});
      } else {
        Store(data, right, {i, j});
      }
    });
  }

  Equation equation;
  Complete left{};
  Complete right{};
};

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::amrex::ScopeGuard guard(argc, argv);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM >= 2);

  const std::array<int, Dim> n_cells{AMREX_D_DECL(128, 128, 1)};
  const std::array<double, Dim> xlower{AMREX_D_DECL(-1.0, -1.0, -1.0)};
  const std::array<double, Dim> xupper{AMREX_D_DECL(+1.0, +1.0, +1.0)};

  fub::PerfectGas<2> equation{};

  fub::amrex::CartesianGridGeometry geometry{};
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 3;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 2, 1)};

  using Complete = fub::PerfectGas<2>::Complete;
  fub::GradientDetector gradient{equation,
                                 std::make_pair(&Complete::density, 0.001),
                                 std::make_pair(&Complete::pressure, 0.1)};

  auto from_prim = [](Complete& state, const fub::PerfectGas<2>& equation) {
    state.energy = state.pressure * equation.gamma_minus_1_inv;
    state.speed_of_sound =
        std::sqrt(equation.gamma * state.pressure / state.density);
  };

  Complete left;
  left.density = 1.0;
  left.momentum = 0.0;
  left.pressure = 8.0;
  from_prim(left, equation);

  Complete right;
  right.density = 1.0;
  right.momentum = 0.0;
  right.pressure = 1.0;
  from_prim(right, equation);

  ShockData initial_data{equation, left, right};

  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(desc, geometry, hier_opts),
      fub::amrex::AdaptInitialData(initial_data, equation),
      fub::amrex::AdaptTagging(equation, gradient, fub::TagBuffer(4)),
      fub::TransmissiveBoundary(equation));
  gridding->InitializeHierarchy(0.0);

  fub::HyperbolicSplitPatchIntegrator patch_integrator{equation};
  fub::GodunovMethod flux_method{equation};

  const int gcw = flux_method.GetStencilWidth();
  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      fub::amrex::HyperbolicSplitIntegratorContext(std::move(gridding), gcw),
      fub::amrex::HyperbolicSplitPatchIntegrator(patch_integrator),
      fub::amrex::FluxMethod(flux_method),
      fub::amrex::Reconstruction(equation)));

  std::string base_name = "PerfectGas2d/";

  auto output = [&](const fub::amrex::PatchHierarchy& hierarchy,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::WritePlotFile(name, hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  auto print_msg = [](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(solver.GetPatchHierarchy(), 0, 0.0s);
  fub::RunOptions run_options{};
  run_options.final_time = 2ms;
  run_options.output_interval = 0.1ms;
  run_options.cfl = 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}
