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

#include "fub/equations/PerfectGas.hpp"

#include "fub/CartesianCoordinates.hpp"
#include "fub/HyperbolicSplitCutCellPatchIntegrator.hpp"
#include "fub/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/HyperbolicSplitSystemSolver.hpp"

#include "fub/grid/AMReX/FillCutCellData.hpp"
#include "fub/grid/AMReX/GriddingAlgorithm.hpp"
#include "fub/grid/AMReX/ScopeGuard.hpp"
#include "fub/grid/AMReX/cutcell/FluxMethod.hpp"
#include "fub/grid/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/grid/AMReX/cutcell/HyperbolicSplitIntegratorContext.hpp"
#include "fub/grid/AMReX/cutcell/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/grid/AMReX/cutcell/IndexSpace.hpp"
#include "fub/grid/AMReX/cutcell/Reconstruction.hpp"
#include "fub/grid/AMReX/cutcell/Tagging.hpp"

#include "fub/geometry/Halfspace.hpp"
#include "fub/initial_data/RiemannProblem.hpp"

#include "fub/tagging/GradientDetector.hpp"
#include "fub/tagging/TagBuffer.hpp"
#include "fub/tagging/TagCutCells.hpp"

#include "fub/boundary_condition/TransmissiveBoundary.hpp"

#include "fub/cutcell_method/KbnStabilisation.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"

#include "fub/RunSimulation.hpp"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <iostream>

#include <xmmintrin.h>

auto Rectangle(const std::array<double, 2>& lower,
               const std::array<double, 2>& upper) {
  amrex::EB2::PlaneIF lower_x({lower[0], lower[1]}, {0, +1});
  amrex::EB2::PlaneIF lower_y({lower[0], lower[1]}, {+1, 0});
  amrex::EB2::PlaneIF upper_x({upper[0], upper[1]}, {0, -1});
  amrex::EB2::PlaneIF upper_y({upper[0], upper[1]}, {-1, 0});
  return amrex::EB2::makeIntersection(lower_x, lower_y, upper_x, upper_y);
}

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard _(argc, argv);

  const std::array<int, AMREX_SPACEDIM> n_cells{AMREX_D_DECL(256, 256, 1)};
  const std::array<double, AMREX_SPACEDIM> xlower{
      AMREX_D_DECL(-0.10, -0.15, -0.15)};
  const std::array<double, AMREX_SPACEDIM> xupper{
      AMREX_D_DECL(+0.20, +0.15, +0.15)};
  amrex::RealBox xbox(xlower, xupper);
  const std::array<int, AMREX_SPACEDIM> periodicity{};

  const int n_level = 2;

  amrex::Geometry coarse_geom(
      amrex::Box{
          {}, {AMREX_D_DECL(n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1)}},
      &xbox, -1, periodicity.data());
  auto embedded_boundary =
      amrex::EB2::makeUnion(Rectangle({-1.0, +0.015}, {0.0, 1.0}),
                            Rectangle({-1.0, -1.0}, {0.0, -0.015}));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::PerfectGas<2> equation;
  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  fub::amrex::cutcell::PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, coarse_geom, n_level);

  auto hierarchy = std::make_shared<fub::amrex::cutcell::PatchHierarchy>(
      desc, geometry, options);

  //  auto hierarchy = fub::amrex::cutcell::ReadCheckpointFile(
  //      "LinearShock2d/Checkpoint", desc, geometry, options);

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.0;
  cons.momentum << 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> right;
  fub::CompleteFromCons(equation, right, cons);

  cons.energy *= 4;
  fub::Complete<fub::PerfectGas<2>> left;
  fub::CompleteFromCons(equation, left, cons);

  fub::amrex::cutcell::RiemannProblem initial_data(
      hierarchy, equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04), left,
      right);

  using State = fub::Complete<fub::PerfectGas<2>>;
  fub::GradientDetector gradients{std::pair{&State::pressure, 0.05},
                                  std::pair{&State::density, 0.005}};

  fub::TransmissiveBoundary boundary{equation};

  fub::HyperbolicSplitCutCellPatchIntegrator patch_integrator{equation};
  fub::MusclHancockMethod muscl_method{equation};
  fub::KbnCutCellMethod cutcell_method(muscl_method);

  auto gridding = std::make_shared<fub::amrex::cutcell::GriddingAlgorithm>(
      hierarchy, fub::amrex::AdaptInitialData(initial_data, equation),
      fub::amrex::cutcell::AdaptTagging(equation, hierarchy, fub::TagCutCells(),
                                        gradients, fub::TagBuffer(4)),
      boundary);

  gridding->InitializeHierarchy(0.0);

  const int gcw = cutcell_method.GetStencilWidth();
  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      fub::amrex::cutcell::HyperbolicSplitIntegratorContext(gridding, gcw),
      fub::amrex::cutcell::HyperbolicSplitPatchIntegrator(patch_integrator),
      fub::amrex::cutcell::FluxMethod(cutcell_method),
      fub::amrex::cutcell::Reconstruction(equation)));

  std::string base_name = "LinearShock2d";

  auto output = [&](auto& hierarchy, std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}/Plot-{:04}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::cutcell::WritePlotFile(name, *hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };
  auto print_msg = [](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(hierarchy, 0, 0s);
  fub::RunOptions run_options{};
  run_options.final_time = 0.002s;
  run_options.output_interval = 0.0000125s;
  run_options.cfl = 0.5 * 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}
