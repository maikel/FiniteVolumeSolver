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

#include "fub/equations/IdealGasMix.hpp"
#include "fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp"

#include "fub/CartesianCoordinates.hpp"
#include "fub/HyperbolicSplitCutCellPatchIntegrator.hpp"
#include "fub/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/HyperbolicSplitSystemSolver.hpp"

#include "fub/AMReX/FillCutCellData.hpp"
#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/ScopeGuard.hpp"
#include "fub/AMReX/cutcell/FluxMethod.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/AMReX/cutcell/HyperbolicSplitIntegratorContext.hpp"
#include "fub/AMReX/cutcell/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/AMReX/cutcell/IndexSpace.hpp"
#include "fub/AMReX/cutcell/Reconstruction.hpp"
#include "fub/AMReX/cutcell/Tagging.hpp"

#include "fub/geometry/ExpandTube.hpp"
#include "fub/geometry/Halfspace.hpp"
#include "fub/initial_data/RiemannProblem.hpp"

#include "fub/tagging/GradientDetector.hpp"
#include "fub/tagging/TagBuffer.hpp"
#include "fub/tagging/TagCutCells.hpp"

#include "fub/boundary_condition/TransmissiveBoundary.hpp"

#include "fub/cutcell_method/KbnStabilisation.hpp"
#include "fub/flux_method/HllMethod.hpp"
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
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard _(argc, argv);

  const std::array<int, 3> n_cells{64, 64, 64};
  const std::array<double, 3> xlower{-0.10, -0.15, -0.15};
  const std::array<double, 3> xupper{+0.20, +0.15, +0.15};
  const std::array<int, 3> periodicity{0, 0, 0};

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(
      amrex::Box{
          {}, {AMREX_D_DECL(n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1)}},
      &xbox, -1, periodicity.data());

  const int n_level = 1;

  auto embedded_boundary = amrex::EB2::makeIntersection(
      amrex::EB2::PlaneIF({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, false),
      amrex::EB2::CylinderIF(0.015, -1.0, 0, {1e6, 0.0, 0.0}, true));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::amrex::cutcell::PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, coarse_geom, n_level);

  fub::Burke2012 mechanism;
  fub::IdealGasMix<3> equation{mechanism};
  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  equation.GetReactor().SetMoleFractions("N2:80,O2:20");
  equation.GetReactor().SetTemperature(300.0);
  equation.GetReactor().SetDensity(1.22);
  fub::Complete<fub::IdealGasMix<3>> right(equation);
  equation.CompleteFromReactor(right);

  equation.GetReactor().SetMoleFractions("N2:80,O2:20");
  equation.GetReactor().SetTemperature(500.0);
  equation.GetReactor().SetDensity(3.15);
  fub::Complete<fub::IdealGasMix<3>> left(equation);
  equation.CompleteFromReactor(left, {400.0, 0.0, 0.0});

  fub::amrex::cutcell::RiemannProblem initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04), left, right);

  using Complete = fub::Complete<fub::IdealGasMix<3>>;
  fub::GradientDetector gradients{equation,
                                  std::pair{&Complete::pressure, 0.05},
                                  std::pair{&Complete::density, 0.005}};

  fub::HyperbolicSplitCutCellPatchIntegrator patch_integrator{equation};

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<3>> signals{};
  fub::Hll riemann_solver{equation, signals};
  fub::MusclHancockMethod flux_method(equation,
                                      fub::HllMethod{equation, signals});
  fub::KbnCutCellMethod cutcell_method{flux_method, riemann_solver};

  auto gridding = std::make_shared<fub::amrex::cutcell::GriddingAlgorithm>(
      fub::amrex::cutcell::PatchHierarchy(desc, geometry, options),
      fub::amrex::cutcell::AdaptInitialData(initial_data, equation),
      fub::amrex::cutcell::AdaptTagging(equation, fub::TagCutCells(), gradients,
                                        fub::TagBuffer(4)),
      fub::TransmissiveBoundary(equation));
  gridding->InitializeHierarchy(0.0);

  const int gcw = cutcell_method.GetStencilWidth();
  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      fub::amrex::cutcell::HyperbolicSplitIntegratorContext(std::move(gridding),
                                                            gcw),
      fub::amrex::cutcell::HyperbolicSplitPatchIntegrator(patch_integrator),
      fub::amrex::cutcell::FluxMethod(std::move(cutcell_method)),
      fub::amrex::cutcell::Reconstruction(equation)));

  std::string base_name = "LinearShock3d/";

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
  run_options.final_time = 0.002s;
  run_options.output_interval = 0.5 * 0.0000125s;
  run_options.cfl = 0.5 * 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}
