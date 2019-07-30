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

#include <xmmintrin.h>

struct RiemannProblem {
  fub::IdealGasMix<1> equation_;
  void
  InitializeData(const fub::View<fub::Complete<fub::IdealGasMix<1>>>& states,
                 const fub::amrex::PatchHierarchy& hierarchy,
                 fub::amrex::PatchHandle patch) {
    const ::amrex::Geometry& geom = hierarchy.GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->tilebox();
    fub::CartesianCoordinates x =
        fub::amrex::GetCartesianCoordinates(geom, box);
    fub::FlameMasterReactor& reactor = equation_.GetReactor();
    reactor.SetMoleFractions("N2:79,O2:21,H2:42");
    const double high_temp = 1150.0;
    const double low_temp = 300.0;
    fub::Complete<fub::IdealGasMix<1>> complete(equation_);
    fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
      const double x0 = x(i)[0];
      const double d = std::clamp(x0 / 0.1, 0.0, 1.0);
      reactor.SetTemperature(d * low_temp + (1.0 - d) * high_temp);
      reactor.SetPressure(101325.0);
      equation_.CompleteFromReactor(complete);
      Store(states, complete, {i});
    });
  }
};

auto MakeTubeSolver(int num_cells, fub::Burke2012& mechanism) {
  const std::array<int, 3> n_cells{num_cells, 1, 1};
  const std::array<double, 3> xlower{0.0, 0.0, 0.0};
  const std::array<double, 3> xupper{+1.0, +1.0, +1.0};

  fub::IdealGasMix<1> equation{fub::FlameMasterReactor(mechanism)};

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 1;
  hier_opts.refine_ratio = amrex::IntVect{2, 1, 1};

  using Complete = fub::IdealGasMix<1>::Complete;
  fub::GradientDetector gradient{equation,
                                 std::make_pair(&Complete::density, 5e-3),
                                 std::make_pair(&Complete::pressure, 5e-2)};

  RiemannProblem initial_data{equation};

  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(desc, geometry, hier_opts),
      fub::amrex::AdaptInitialData(initial_data, equation),
      fub::amrex::AdaptTagging(equation, gradient),
      fub::TransmissiveBoundary(equation));
  gridding->InitializeHierarchy(0.0);

  fub::HyperbolicSplitPatchIntegrator patch_integrator{equation};
  fub::HllMethod hlle(equation,
                      fub::EinfeldtSignalVelocities<fub::IdealGasMix<1>>{});
  fub::MusclHancockMethod flux_method{equation, hlle};

  const int gcw = flux_method.GetStencilWidth();
  fub::HyperbolicSplitSystemSolver system_solver(
      fub::HyperbolicSplitLevelIntegrator(
          fub::amrex::HyperbolicSplitIntegratorContext(gridding, gcw),
          fub::amrex::HyperbolicSplitPatchIntegrator(patch_integrator),
          fub::amrex::FluxMethod(flux_method),
          fub::amrex::Reconstruction(equation)));

  fub::ideal_gas::KineticSourceTerm<1> source_term(equation, gridding);
  fub::SplitSystemSourceSolver solver(system_solver, source_term);

  return solver;
}

auto MakePlenumSolver(int num_cells, fub::Burke2012& mechanism) {
  const std::array<int, 3> n_cells{num_cells, num_cells, num_cells};
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

  fub::IdealGasMix<3> equation{mechanism};
  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  equation.GetReactor().SetMoleFractions("N2:79,O2:21");
  equation.GetReactor().SetTemperature(300.0);
  equation.GetReactor().SetPressure(101325.0);
  fub::Complete<fub::IdealGasMix<3>> right(equation);
  equation.CompleteFromReactor(right);

  fub::amrex::cutcell::RiemannProblem initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04), right, right);

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

  return solver;
}

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard _(argc, argv);
  fub::Burke2012 mechanism{};
  auto plenum = MakePlenumSolver(96, mechanism);
  auto tube = MakeTubeSolver(400, mechanism);

  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO |
                         _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW |
                         _MM_MASK_INVALID);

  fub::amrex::CoupledBoundaryFunction coupled_boundary(
      plenum.GetPatchHierarchy(), tube.GetPatchHierarchy(), 3,
      plenum.GetEquation().GetReactor());

  {
    fub::TransmissiveBoundary transmissive{tube.GetEquation()};
    fub::BoundarySet<fub::amrex::PatchHierarchy> boundaries;
    boundaries.SetBoundaryCondition(fub::Location{0, 0}, transmissive);
    boundaries.SetBoundaryCondition(fub::Location{0, 1}, coupled_boundary);
    tube.GetGriddingAlgorithm()->SetBoundaryCondition(std::move(boundaries));
  }

  {
    fub::TransmissiveBoundary transmissive{plenum.GetEquation()};
    fub::BoundarySet<fub::amrex::cutcell::PatchHierarchy> boundaries;
    boundaries.SetBoundaryCondition(fub::Location{0, 0}, coupled_boundary);
    boundaries.SetBoundaryCondition(fub::Location{0, 1}, transmissive);
    boundaries.SetBoundaryCondition(fub::Location{1, 0}, transmissive);
    boundaries.SetBoundaryCondition(fub::Location{1, 1}, transmissive);
    boundaries.SetBoundaryCondition(fub::Location{2, 0}, transmissive);
    boundaries.SetBoundaryCondition(fub::Location{2, 1}, transmissive);
    plenum.GetGriddingAlgorithm()->SetBoundaryCondition(std::move(boundaries));
  }

  std::shared_ptr<fub::amrex::CoupledBoundary> boundary =
      coupled_boundary.GetSharedState();

  std::string base_name = "CoupledLinearShock3d";
  auto output = [&](const fub::amrex::cutcell::PatchHierarchy& hierarchy,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}_Plenum/{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::cutcell::WritePlotFile(name, hierarchy, plenum.GetEquation());
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
    name = fmt::format("{}_Tube/{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::WritePlotFile(name, tube.GetPatchHierarchy(),
                              tube.GetEquation());
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  auto print_msg = [&](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(plenum.GetPatchHierarchy(), plenum.GetCycles(), plenum.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.final_time = 0.002s;
  run_options.output_interval = 1.25e-5s;
  run_options.cfl = 0.5 * 0.9;

  fub::amrex::RunCoupledSimulation(plenum, tube, *boundary, run_options,
                                   wall_time_reference, output, print_msg);
}
