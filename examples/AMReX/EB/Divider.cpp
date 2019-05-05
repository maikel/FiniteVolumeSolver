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

#include "fub/geometry/ExpandTube.hpp"
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
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <boost/filesystem.hpp>

#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
  static_assert(AMREX_SPACEDIM == 3);

  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  const fub::amrex::ScopeGuard _(argc, argv);

  const int r = 3;

  const std::array<int, 3> n_cells{8 * r * 12, 8 * r * 3, 8 * r * 2};
  const std::array<double, 3> xlower{0.0, -0.01, -0.02};
  const std::array<double, 3> xupper{0.23, 0.05, 0.02};
  const std::array<int, 3> periodicity{0, 0, 0};

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(
      amrex::Box{
          {}, {AMREX_D_DECL(n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1)}},
      &xbox, -1, periodicity.data());

  const int n_level = 1;

  Eigen::Matrix<double, 3, Eigen::Dynamic> centerline =
      fub::ReadPointsFromFile("centerline.txt");

  const double radius = 0.5 * 0.0125;
  auto embedded_boundary = amrex::EB2::makeComplement(amrex::EB2::makeUnion(
      amrex::EB2::CylinderIF(radius, 1.0, 0, {-0.1, radius, 0.0}, false),
      fub::ExpandTube(centerline, radius)));
  //  auto embedded_boundary = amrex::EB2::CylinderIF(radius, 1.0, 0, {-0.1,
  //  radius, 0.0}, false);
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::PerfectGas<3> equation;
  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  fub::amrex::cutcell::PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, coarse_geom, n_level);

  fub::Conservative<fub::PerfectGas<3>> cons;
  cons.density = 1.22;
  cons.momentum << 0.0, 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<3>> right;
  fub::CompleteFromCons(equation, right, cons);

  // cons.density = 3.15736;
  // cons.momentum << 1258.31, 0.0, 0.0;
  // cons.energy = 416595.0 * equation.gamma_minus_1_inv +
  //               0.5 * cons.momentum[0] * cons.momentum[0] / cons.density;
  cons.energy *= 8.0;
  fub::Complete<fub::PerfectGas<3>> left;
  fub::CompleteFromCons(equation, left, cons);

  fub::amrex::cutcell::RiemannProblem initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.04), left, right);

  using State = fub::Complete<fub::PerfectGas<3>>;
  fub::GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                                  std::pair{&State::density, 0.005}};

  fub::HyperbolicSplitCutCellPatchIntegrator patch_integrator{equation};
  fub::MusclHancockMethod flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(std::move(flux_method));

  auto gridding = std::make_shared<fub::amrex::cutcell::GriddingAlgorithm>(
      fub::amrex::cutcell::PatchHierarchy(desc, geometry, options),
      fub::amrex::cutcell::AdaptInitialData(initial_data, equation),
      fub::amrex::cutcell::AdaptTagging(equation, fub::TagCutCells(), gradients,
                                        fub::TagBuffer(4)),
      fub::TransmissiveBoundary{equation});
  gridding->InitializeHierarchy(0.0);

  const int gcw = cutcell_method.GetStencilWidth();
  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      fub::amrex::cutcell::HyperbolicSplitIntegratorContext(std::move(gridding),
                                                            gcw),
      fub::amrex::cutcell::HyperbolicSplitPatchIntegrator(patch_integrator),
      fub::amrex::cutcell::FluxMethod(cutcell_method),
      fub::amrex::cutcell::Reconstruction(equation)));

  std::string base_name = "Divider/";

  auto output = [&](const fub::amrex::cutcell::PatchHierarchy& hierarchy,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::cutcell::WritePlotFile(name, hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  auto print_msg = [&](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(solver.GetPatchHierarchy(), solver.GetCycles(), solver.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.final_time = 0.0005s;
  run_options.output_interval = 0.5 * 0.0000125s;
  run_options.cfl = 0.5 * 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}
