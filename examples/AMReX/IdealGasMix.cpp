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

#include <fub/AMReX/IgniteDetonation.hpp>
#include <fub/AMReX/boundary_condition/PressureValveBoundary.hpp>
#include <xmmintrin.h>

struct ProgramOptions {
  double final_time{0.20};
  double cfl{0.8};
  int n_cells{200};
  int max_refinement_level{1};
  double domain_length{1.5};
  double output_interval{1.0E-5};
};

std::optional<std::pair<ProgramOptions, boost::program_options::variables_map>>
ParseCommandLine(int argc, char** argv) {
  namespace po = boost::program_options;
  ProgramOptions opts;
  po::options_description desc("Program Options");
  desc.add_options()("help", "Print help messages")(
      "cfl", po::value<double>(&opts.cfl)->default_value(opts.cfl),
      "Set the CFL condition")("max_refinement_level",
                               po::value<int>(&opts.max_refinement_level)
                                   ->default_value(opts.max_refinement_level),
                               "Set the maximal refinement level")(
      "n_cells", po::value<int>(&opts.n_cells)->default_value(opts.n_cells),
      "Set number of cells in each direction for the plenum")(
      "final_time",
      po::value<double>(&opts.final_time)->default_value(opts.final_time),
      "Set the final simulation time")(
      "domain_length",
      po::value<double>(&opts.domain_length)->default_value(opts.domain_length),
      "Set the base length for the tube")(
      "output_interval",
      po::value<double>(&opts.output_interval)
          ->default_value(opts.output_interval),
      "Sets the output interval");
  desc.add(fub::amrex::PressureValveOptions::GetCommandLineOptions("valve"));
  desc.add(fub::amrex::IgniteDetonationOptions::GetCommandLineOptions());
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (std::exception& e) {
    amrex::Print()
        << "[Error] An Error occured while reading program options:\n";
    amrex::Print() << e.what() << '\n';
    return {};
  }

  if (vm.count("help")) {
    amrex::Print() << desc << "\n";
    return {};
  }

  opts.max_refinement_level = std::max(1, opts.max_refinement_level);

  amrex::Print() << fmt::format(
      "[Info] Simulation will run with the following options:\n[Info]\n");
  amrex::Print() << fmt::format("[Info] final_time = {}s\n", opts.final_time);
  amrex::Print() << fmt::format("[Info] output_interval = {}s\n",
                                opts.output_interval);
  amrex::Print() << fmt::format("[Info] cfl = {}\n", opts.cfl);
  amrex::Print() << fmt::format("[Info] max_refinement_level = {}\n",
                                opts.max_refinement_level);
  std::array<double, 3> xlo{0.0, 0.0, 0.0};
  std::array<double, 3> xup{opts.domain_length, +0.03, +0.03};
  amrex::Print() << fmt::format(
      "[Info] domain = [{}, {}] x [{}, {}] x [{}, {}]\n", xlo[0], xup[0],
      xlo[1], xup[1], xlo[2], xup[2]);
  amrex::Print() << fmt::format("[Info] n_cells = {}\n[Info] dx = {}\n",
                                opts.n_cells,
                                opts.domain_length / opts.n_cells);

  return std::make_pair(opts, vm);
}

void MyMain(const ProgramOptions& opts,
            const boost::program_options::variables_map& map) {
  // Store a reference timepoint to measure the wall time duration
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  // Enable floating point exceptions.
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO |
                         _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW |
                         _MM_MASK_INVALID);

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
  fub::amrex::PressureValveOptions valve_options(map, "valve");
  boundary.conditions.push_back(PressureValveBoundary{equation, valve_options});
  boundary.conditions.push_back(
      IsentropicPressureBoundary{equation, 101325.0, fub::Direction::X, 1});

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = nlevels;
  hier_opts.refine_ratio = ::amrex::IntVect{AMREX_D_DECL(2, 1, 1)};

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
  //  fub::EinfeldtSignalVelocities<fub::IdealGasMix<1>> signals{};
  //  fub::HllMethod hll_method(equation, signals);

  fub::ideal_gas::MusclHancockPrimMethod<1> flux_method(equation);

  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethod(fub::execution::seq, flux_method),
      fub::amrex::ForwardIntegrator(fub::execution::seq),
      fub::amrex::Reconstruction(fub::execution::seq, equation)};

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<1>, fub::amrex::IntegratorContext(gridding, method),
      fub::GodunovSplitting());

  fub::amrex::IgniteDetonationOptions io(map);
  fub::amrex::IgniteDetonation ignite(equation, gridding, io);
  fub::DimensionalSplitSystemSourceSolver ign_solver(system_solver, ignite,
                                                     fub::GodunovSplitting{});

  fub::ideal_gas::KineticSourceTerm<1> source_term(equation, gridding);

  fub::DimensionalSplitSystemSourceSolver solver(ign_solver, source_term,
                                                 fub::StrangSplitting());
  // }}}

  // Run the simulation with given feedback functions

  std::string base_name = "IdealGasMix/";
  int rank = -1;
  MPI_Comm_rank(solver.GetContext().GetMpiCommunicator(), &rank);
  auto output =
      [&](const std::shared_ptr<fub::amrex::GriddingAlgorithm>& gridding,
          std::ptrdiff_t cycle, fub::Duration timepoint, int = 0) {
        std::string name = fmt::format("{}plt{:05}", base_name, cycle);
        amrex::Print() << "Start output to '" << name << "'.\n";
        fub::amrex::WritePlotFile(name, gridding->GetPatchHierarchy(),
                                  equation);

        fub::amrex::WriteTubeData(
            fmt::format("{}/Matlab/{:07}.dat", base_name, cycle),
            gridding->GetPatchHierarchy(), equation, timepoint, cycle,
            solver.GetContext().GetMpiCommunicator());
        amrex::Print() << "Finished output to '" << name << "'.\n";
      };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.cfl = opts.cfl;
  run_options.final_time = fub::Duration(opts.final_time);
  run_options.output_interval = {fub::Duration(opts.output_interval)};
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
  {
    fub::amrex::ScopeGuard _{};
    std::optional<
        std::pair<ProgramOptions, boost::program_options::variables_map>>
        po = ParseCommandLine(argc, argv);
    if (po) {
      MyMain(std::get<0>(*po), std::get<1>(*po));
    }
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}
