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

#include <xmmintrin.h>

struct TemperatureRamp {
  using Complete = fub::IdealGasMix<1>::Complete;

  fub::IdealGasMix<1> equation_;
  double ignition_pos_{-1.0};
  double air_position_{0.0};
  double equiv_raito_{1.0};

  void InitializeData(::amrex::MultiFab& data, const ::amrex::Geometry& geom) {
    fub::FlameMasterReactor& reactor = equation_.GetReactor();
    const double high_temp = 1400.0;
    const double temp = 300.0;
    const double low_pressure = 101325.0;
    const double high_pressure = 8.0 * low_pressure;
    Complete left(equation_);
    Complete right{equation_};
    Complete right_air{equation_};

    double h2_moles = 42.0 * equiv_raito_;
    //    double o2_moles = std::max(21.0 - 0.5 * h2_moles, 0.0);

    std::string fuel_moles = fmt::format("N2:79,O2:21,H2:{}", h2_moles);
    reactor.SetMoleFractions(fuel_moles);
    reactor.SetTemperature(high_temp);
    reactor.SetPressure(high_pressure);
    equation_.CompleteFromReactor(left);

    reactor.SetMoleFractions("N2:79,O2:21");
    reactor.SetTemperature(temp);
    reactor.SetPressure(low_pressure);
    equation_.CompleteFromReactor(right_air);

    reactor.SetMoleFractions(fuel_moles);
    reactor.SetTemperature(temp);
    reactor.SetPressure(low_pressure);
    equation_.CompleteFromReactor(right);

    fub::amrex::ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
        double x = geom.CellCenter(i, 0);
        if (x < ignition_pos_) {
          fub::Store(states, left, {i});
        } else if (x < air_position_) {
          fub::Store(states, right, {i});
        } else {
          fub::Store(states, right_air, {i});
        }
      });
    });
  }
};

struct ProgramOptions {
  double final_time{0.20};
  double cfl{0.4};
  int n_cells{200};
  int max_refinement_level{1};
  double domain_length{1.5};
  double ignition_position{0.3};
  double air_position{0.5 + 0.35};
  double equiv_ratio{1.0};
  double output_interval{1.0E-5};
  std::ptrdiff_t output_frequency{0};
};

std::optional<ProgramOptions> ParseCommandLine(int argc, char** argv) {
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
      "ignition_pos",
      po::value<double>(&opts.ignition_position)
          ->default_value(opts.ignition_position),
      "Set the position for the ignition of the detonation inside the tube")(
      "air_buffer_start",
      po::value<double>(&opts.air_position)->default_value(opts.air_position),
      "Sets the starting position for an air buffer in the tube.")(
      "equiv_ratio",
      po::value<double>(&opts.equiv_ratio)->default_value(opts.equiv_ratio),
      "Sets the equivalence ratio of the fuel in the tube")(
      "output_interval",
      po::value<double>(&opts.output_interval)
          ->default_value(opts.output_interval),
      "Sets the output interval")(
      "output_frequency",
      po::value<std::ptrdiff_t>(&opts.output_frequency)
          ->default_value(opts.output_frequency),
      "Sets the output frequency");
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
  amrex::Print() << fmt::format("[Info] ignition_pos = {}\n",
                                opts.ignition_position);
  amrex::Print() << fmt::format("[Info] equiv_ratio = {}\n", opts.equiv_ratio);
  amrex::Print() << fmt::format("[Info] air_position = {}\n",
                                opts.air_position);

  return opts;
}

void MyMain(const ProgramOptions& opts) {
  // feenableexcept(FE_DIVBYZERO | FE_INVALID);
  // Store a reference timepoint to measure the wall time duration
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::InitializeLogging(MPI_COMM_WORLD);

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

  TemperatureRamp initial_data{equation};
  initial_data.equiv_raito_ = opts.equiv_ratio;
  initial_data.ignition_pos_ = opts.ignition_position;
  initial_data.air_position_ = opts.air_position;

  using Complete = fub::IdealGasMix<1>::Complete;
  fub::amrex::GradientDetector gradient{
      equation, std::make_pair(&Complete::pressure, 1e-3),
      std::make_pair(&Complete::density, 1e-3),
      std::make_pair(&Complete::temperature, 1e-1)};

  fub::amrex::BoundarySet boundary;
  using fub::amrex::IsentropicPressureBoundary;
  using fub::amrex::ReflectiveBoundary;
  using fub::amrex::TransmissiveBoundary;
  boundary.conditions.push_back(
      ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::X, 0});
  //  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 1});
  boundary.conditions.push_back(
      IsentropicPressureBoundary{equation, 101325.0, fub::Direction::X, 1});

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = nlevels;
  hier_opts.refine_ratio = ::amrex::IntVect{AMREX_D_DECL(2, 1, 1)};
  hier_opts.blocking_factor = amrex::IntVect(AMREX_D_DECL(8, 1, 1));

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts), initial_data,
      gradient, boundary);
  gridding->InitializeHierarchy(0.0);
  // }}}

  // Setup the numerical Method used to solve this problem.
  // {{{
  fub::EinfeldtSignalVelocities<fub::IdealGasMix<1>> signals{};
  fub::HllMethod hll_method(equation, signals);

  fub::amrex::HyperbolicMethod method{fub::amrex::FluxMethod(hll_method),
                                      fub::amrex::ForwardIntegrator(),
                                      fub::amrex::Reconstruction(equation)};

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<1>, fub::amrex::IntegratorContext(gridding, method, 2, 1),
      fub::GodunovSplitting());

  std::shared_ptr<fub::CounterRegistry> registry =
      system_solver.GetCounterRegistry();

  auto diameter = [](double x) -> double {
    if (x < 0.5) {
      return 0.45;
    }
    if (x < 0.6) {
      const double lambda = std::clamp((0.6 - x) / 0.1, 0.0, 1.0);
      return lambda * 0.45 + (1.0 - lambda) * 0.3;
    }
    return 0.3;
  };

  fub::amrex::AxialSourceTerm source_term(equation, diameter,
                                          system_solver.GetGriddingAlgorithm());

  fub::SplitSystemSourceLevelIntegrator axial_solver(std::move(system_solver),
                                                     std::move(source_term),
                                                     fub::GodunovSplitting());

  fub::ideal_gas::KineticSourceTerm<1> kinetic_source(equation);

  fub::SplitSystemSourceLevelIntegrator level_integrator(
      std::move(axial_solver), std::move(kinetic_source),
      fub::GodunovSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
  // }}}

  // Run the simulation with given feedback functions

  std::string base_name = "AxialSourceTerm/";

  int rank = -1;
  MPI_Comm_rank(solver.GetMpiCommunicator(), &rank);
  using namespace std::literals::chrono_literals;
  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> output{};
  output.AddOutput(fub::MakeOutput<fub::amrex::GriddingAlgorithm>(
      {}, {fub::Duration(opts.output_interval)},
      [&](const fub::amrex::GriddingAlgorithm& gridding) {
        std::ptrdiff_t cycle = gridding.GetCycles();
        fub::Duration timepoint = gridding.GetTimePoint();
        std::string name = fmt::format("{}plt{:05}", base_name, cycle);
        amrex::Print() << "Start output to '" << name << "'.\n";
        fub::amrex::WritePlotFile(name, gridding.GetPatchHierarchy(), equation);
        fub::amrex::WriteTubeData(
            base_name + "/Tube.h5", gridding.GetPatchHierarchy(), equation,
            timepoint, cycle, solver.GetMpiCommunicator());
        amrex::Print() << "Finished output to '" << name << "'.\n";
      }));
  output.AddOutput(
      std::make_unique<fub::CounterOutput<fub::amrex::GriddingAlgorithm>>(
          wall_time_reference, std::vector<std::ptrdiff_t>{},
          std::vector<fub::Duration>{0.0001s}));

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.cfl = opts.cfl;
  run_options.final_time = fub::Duration(opts.final_time);
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  {
    fub::amrex::ScopeGuard _{};
    std::optional<ProgramOptions> po = ParseCommandLine(argc, argv);
    if (po) {
      MyMain(*po);
    }
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}
