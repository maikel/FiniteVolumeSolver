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

#include "fub/AMReX.hpp"
#include "fub/Solver.hpp"

#include <boost/program_options.hpp>

#include <fmt/format.h>

#include <iostream>

struct ProgramOptions {
  double final_time{0.20};
  double cfl{0.8};
  int n_cells{200};
  int max_refinement_level{1};
  double domain_length{1.5};
  double output_interval{1.0E-5};
};

void MyMain(const ProgramOptions& opts) {
  // Store a reference timepoint to measure the wall time duration
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  constexpr int Dim = AMREX_SPACEDIM;

  // Setup the domain parameters
  const std::array<int, Dim> n_cells{AMREX_D_DECL(opts.n_cells, 1, 1)};
  const int nlevels = opts.max_refinement_level;
  const std::array<double, Dim> xlower{AMREX_D_DECL(0.0, 0.0, 0.0)};
  const std::array<double, Dim> xupper{
      AMREX_D_DECL(opts.domain_length, +0.03, +0.03)};

  // Define the equation which will be solved
  fub::Burke2012 mechanism{};
  fub::IdealGasMix<1> equation{fub::FlameMasterReactor(mechanism)};

  // Define the GriddingAlgorithm for this simulation and initialize data.
  // {{{
  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  using Complete = fub::IdealGasMix<1>::Complete;
  fub::amrex::GradientDetector gradient{
      equation, std::make_pair(&Complete::pressure, 1e-3),
      std::make_pair(&Complete::density, 1e-3),
      std::make_pair(&Complete::temperature, 1e-1)};

  fub::amrex::BoundarySet boundary;
  using fub::amrex::IsentropicPressureBoundary;
  using fub::amrex::PressureValveBoundary;
  fub::amrex::PressureValveOptions valve_options{};
  boundary.conditions.push_back(PressureValveBoundary{equation, valve_options});
  boundary.conditions.push_back(
      IsentropicPressureBoundary{equation, 101325.0, fub::Direction::X, 1});

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = nlevels;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 1, 1)};
  hier_opts.blocking_factor = amrex::IntVect{AMREX_D_DECL(8, 1, 1)};

  fub::IdealGasMix<1>::Complete state(equation);
  equation.GetReactor().SetMoleFractions("N2:79,O2:21");
  equation.GetReactor().SetTemperature(300.0);
  equation.GetReactor().SetPressure(101325.0);
  equation.CompleteFromReactor(state);

  using fub::amrex::ConstantData;
  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts),
      ConstantData{equation, state}, gradient, boundary);
  gridding->InitializeHierarchy(0.0);
  // }}}

  // Setup the numerical Method used to solve this problem.
  // {{{
  fub::ideal_gas::MusclHancockPrimMethod<1> flux_method(equation);

  fub::amrex::HyperbolicMethod method{fub::amrex::FluxMethodAdapter(flux_method),
                                      fub::amrex::EulerForwardTimeIntegrator(),
                                      fub::amrex::Reconstruction(equation)};

  const int scratch_gcw = 2 * flux_method.GetStencilWidth();
  const int flux_gcw = flux_method.GetStencilWidth();

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<1>,
      fub::amrex::IntegratorContext(gridding, method, scratch_gcw, flux_gcw),
      fub::GodunovSplitting());

  fub::amrex::IgniteDetonationOptions io{};
  io.ignite_interval = fub::Duration(0.01);
  fub::amrex::IgniteDetonation ignite(equation, hier_opts.max_number_of_levels,
                                      io);

  fub::SplitSystemSourceLevelIntegrator ign_solver(
      std::move(system_solver), std::move(ignite), fub::GodunovSplitting{});

  fub::ideal_gas::KineticSourceTerm<1> source_term(equation);

  fub::SplitSystemSourceLevelIntegrator level_integrator(
      std::move(ign_solver), std::move(source_term), fub::StrangSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
  // }}}

  // Run the simulation with given feedback functions

  std::string base_name = "IdealGasMix/";
  int rank = -1;
  MPI_Comm_rank(solver.GetMpiCommunicator(), &rank);

  using namespace std::literals::chrono_literals;
  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> output{};
  output.AddOutput(fub::MakeOutput<fub::amrex::GriddingAlgorithm>(
      {}, {0.00025s}, fub::amrex::PlotfileOutput(equation, base_name)));
  output.AddOutput(
      std::make_unique<fub::CounterOutput<fub::amrex::GriddingAlgorithm>>(
          wall_time_reference, std::vector<std::ptrdiff_t>{},
          std::vector<fub::Duration>{0.001s}));

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.cfl = opts.cfl;
  run_options.final_time = fub::Duration(opts.final_time);
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}

int main() {
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
  {
    fub::amrex::ScopeGuard _{};
    MyMain({});
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}
