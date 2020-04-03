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

void MyMain(const fub::ProgramOptions& opts) {
  // Store a reference timepoint to measure the wall time duration
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  using namespace fub::amrex;
  ScopeGuard scope_guard{};

  constexpr int TubeDim = 1;

  // Define the equation which will be solved
  fub::Burke2012 mechanism{};
  fub::IdealGasMix<TubeDim> equation{fub::FlameMasterReactor(mechanism)};

  fub::SeverityLogger info = fub::GetInfoLogger();

  // Define the GriddingAlgorithm for this simulation and initialize data.

  fub::amrex::CartesianGridGeometry grid_geometry =
      fub::GetOptions(opts, "CartesianGridGeometry");
  BOOST_LOG(info) << "CartesianGridGeometry: ";
  grid_geometry.Print(info);

  PatchHierarchyOptions hierarchy_options =
      fub::GetOptions(opts, "PatchHierarchy");
  BOOST_LOG(info) << "PatchHierarchy: ";
  hierarchy_options.Print(info);

  using Complete = fub::IdealGasMix<1>::Complete;
  fub::amrex::GradientDetector gradient{
      equation, std::make_pair(&Complete::pressure, 1e-3),
      std::make_pair(&Complete::density, 1e-3),
      std::make_pair(&Complete::temperature, 1e-1)};

  BoundarySet boundary;

  PressureValveOptions valve_options = fub::GetOptions(opts, "PressureValveBoundary");
  BOOST_LOG(info) << "PressureValveBoundary:";
  valve_options.Print(info);

  IsentropicPressureBoundaryOptions right_boundary =
      fub::GetOptions(opts, "IsentropicPressureBoundary");
  BOOST_LOG(info) << "IsentropicPressureBoundary: ";
  right_boundary.Print(info);

  boundary.conditions.push_back(
      PressureValveBoundary{equation, valve_options});
  boundary.conditions.push_back(
      IsentropicPressureBoundary{equation, right_boundary});

  fub::IdealGasMix<1>::Complete state(equation);
  equation.GetReactor().SetMoleFractions("N2:79,O2:21");
  equation.GetReactor().SetTemperature(300.0);
  equation.GetReactor().SetPressure(101325.0);
  equation.CompleteFromReactor(state);

  using fub::amrex::ConstantData;
  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, grid_geometry, hierarchy_options),
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
      fub::GodunovSplitting{});

  fub::amrex::IgniteDetonationOptions io(fub::GetOptions(opts, "IgniteDetonation"));
  BOOST_LOG(info) << "IgniteDetonation:";
  io.Print(info);
  fub::amrex::IgniteDetonation ignite(equation, hierarchy_options.max_number_of_levels,
                                      io);

  fub::SplitSystemSourceLevelIntegrator ign_solver(
      std::move(system_solver), std::move(ignite), fub::GodunovSplitting{});

  fub::ideal_gas::KineticSourceTerm<1> source_term(equation);

  fub::SplitSystemSourceLevelIntegrator level_integrator(
      std::move(ign_solver), std::move(source_term), fub::StrangSplittingLumped());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
  // }}}

  // Run the simulation with given feedback functions

  // Run the simulation with given feedback functions

  using namespace std::literals::chrono_literals;
  using Plotfile = PlotfileOutput<fub::IdealGasMix<1>>;
  using CounterOutput =
      fub::CounterOutput<GriddingAlgorithm, std::chrono::milliseconds>;
  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<Plotfile>("Plotfile", equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<WriteHdf5>("HDF5");
  fub::MultipleOutputs<GriddingAlgorithm> output(
      std::move(factory), fub::GetOptions(opts, "Output"));

  fub::RunOptions run_options = fub::GetOptions(opts, "RunOptions");
  BOOST_LOG(info) << "RunOptions: ";
  run_options.Print(info);

  output(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
  pybind11::scoped_interpreter interpreter{};
  std::optional<fub::ProgramOptions> opts = fub::ParseCommandLine(argc, argv);
  if (opts) {
    MyMain(*opts);
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}
