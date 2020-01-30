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
#include "fub/AMReX_CutCell.hpp"
#include "fub/Solver.hpp"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/algorithm/string/replace.hpp>
#include <boost/core/null_deleter.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/manipulators/add_value.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/file.hpp>

#include <pybind11/stl.h>

#include <cmath>
#include <iostream>

#include <xmmintrin.h>

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = AMREX_SPACEDIM;

static constexpr double r_tube = 0.015;

struct ProgramOptions {
  ProgramOptions() = default;

  ProgramOptions(const fub::ProgramOptions& plenum) {
    checkpoint = fub::GetOptionOr(plenum, "checkpoint", checkpoint);
    plenum_temperature =
        fub::GetOptionOr(plenum, "temperature", plenum_temperature);
    plenum_radius = fub::GetOptionOr(plenum, "radius", plenum_radius);
    plenum_outlet_radius =
        fub::GetOptionOr(plenum, "outlet_radius", plenum_outlet_radius);
    plenum_outlet_exit_radius = fub::GetOptionOr(plenum, "outlet_exit_radius",
                                                 plenum_outlet_exit_radius);
    plenum_outlet_exit_length = fub::GetOptionOr(plenum, "outlet_exit_length",
                                                 plenum_outlet_exit_length);
    plenum_outlet_length =
        fub::GetOptionOr(plenum, "outlet_length", plenum_outlet_length);
    plenum_length = fub::GetOptionOr(plenum, "length", plenum_length);
  }

  template <typename Logger> void Print(Logger& log) const {
    BOOST_LOG(log) << "Problem Options:";
    BOOST_LOG(log) << "  - plenum.radius = " << plenum_radius << " [m]";
    BOOST_LOG(log) << "  - plenum.temperature = " << plenum_temperature
                   << " [K]";
    BOOST_LOG(log) << "  - plenum.outlet_radius= " << plenum_outlet_radius
                   << " [m]";
    BOOST_LOG(log) << "  - plenum.outlet_length = " << plenum_outlet_length
                   << " [m]";
    BOOST_LOG(log) << "  - plenum.outlet_exit_radius = "
                   << plenum_outlet_exit_radius << " [m]";
    BOOST_LOG(log) << "  - plenum.outlet_exit_length = "
                   << plenum_outlet_exit_length << " [m]";
    BOOST_LOG(log) << "  - plenum.length = " << plenum_length << " [m]";

    if (!checkpoint.empty()) {
      BOOST_LOG(log) << "Restart simulation from checkpoint '" << checkpoint
                     << "'!";
    }
  }

  std::string checkpoint{};
  double plenum_radius{0.5};
  double plenum_temperature{300};
  double plenum_outlet_radius{r_tube};
  double plenum_outlet_exit_radius{plenum_outlet_radius};
  double plenum_outlet_exit_length{plenum_outlet_radius};
  double plenum_length{0.25};
  double plenum_outlet_length{0.04 + 0.04 + 0.03};
};

auto MakeTubeSolver(fub::Burke2012& mechanism,
                    const fub::ProgramOptions& options,
                    const std::shared_ptr<fub::CounterRegistry>& counters) {
  using namespace fub::amrex;

  CartesianGridGeometry grid_geometry(fub::GetOptions(options, "GridGeometry"));

  PatchHierarchyOptions hierarchy_options(
      fub::GetOptions(options, "PatchHierarchy"));

  using Complete = fub::IdealGasMix<Tube_Rank>::Complete;
  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 1e-3),
                            std::make_pair(&Complete::pressure, 1e-2),
                            std::make_pair(&Complete::temperature, 1e-1)};

  DataDescription desc = MakeDataDescription(equation);

  ::amrex::Box refine_box{{grid_geometry.cell_dimensions[0] - 5, 0, 0},
                          {grid_geometry.cell_dimensions[0] - 1, 0, 0}};
  ConstantBox constant_box{refine_box};

  Complete state(equation);
  {
    using namespace std::literals;
    const fub::ProgramOptions initial_options =
        fub::GetOptions(options, "InitialCondition");
    const std::string moles =
        fub::GetOptionOr(initial_options, "moles", "N2:79,O2:21"s);
    const double temperature =
        fub::GetOptionOr(initial_options, "temperature", 300.0);
    const double pressure =
        fub::GetOptionOr(initial_options, "pressure", 101325.0);
    equation.GetReactor().SetMoleFractions(moles);
    equation.GetReactor().SetTemperature(temperature);
    equation.GetReactor().SetPressure(pressure);
    equation.CompleteFromReactor(state);
  }
  ConstantData initial_data{equation, state};

  PressureValveOptions valve_options =
      fub::GetOptions(options, "PressureValveBoundary");
  PressureValveBoundary valve{equation, valve_options};
  BoundarySet boundaries{{valve}};

  // If a checkpoint path is specified we will fill the patch hierarchy with
  // data from the checkpoint file, otherwise we will initialize the data by
  // the initial data function.
  std::shared_ptr<GriddingAlgorithm> gridding = [&] {
    std::string checkpoint =
        fub::GetOptionOr(options, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(desc, grid_geometry, hierarchy_options), initial_data,
          TagAllOf(gradient, constant_box), boundaries);
      gridding->GetPatchHierarchy().SetCounterRegistry(counters);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint = fmt::format("{}/Tube", checkpoint);
      PatchHierarchy h = ReadCheckpointFile(checkpoint, desc, grid_geometry,
                                            hierarchy_options);
      std::shared_ptr<GriddingAlgorithm> gridding =
          std::make_shared<GriddingAlgorithm>(std::move(h), initial_data,
                                              TagAllOf(gradient, constant_box),
                                              boundaries);
      return gridding;
    }
  }();

  fub::ideal_gas::MusclHancockPrimMethod<Tube_Rank> flux_method{equation};
  HyperbolicMethod method{FluxMethod(flux_method), EulerForwardTimeIntegrator(),
                          Reconstruction(equation)};

  const int scratch_gcw = 4;
  const int flux_gcw = 2;

  IntegratorContext context(gridding, method, scratch_gcw, flux_gcw);

  return std::pair{std::move(context), valve};
}

auto MakePlenumSolver(fub::Burke2012& mechanism,
                      const fub::ProgramOptions& options) {
  using namespace fub::amrex::cutcell;
  ProgramOptions po = options;
  auto embedded_boundary = amrex::EB2::makeIntersection(
      amrex::EB2::CylinderIF(po.plenum_radius, po.plenum_length, 0,
                             {0.5 * po.plenum_length, 0.0, 0.0}, true),
      amrex::EB2::CylinderIF(r_tube, 0.2, 0, {1e-6, 0.0, 0.0}, true),
      amrex::EB2::CylinderIF(po.plenum_outlet_radius, 0.2, 0,
                             {po.plenum_length, 0.0, 0.0}, true),
      fub::amrex::Geometry(fub::ConeIF({po.plenum_length, 0.0, 0.0},
                                       po.plenum_radius, 0.04, true)),
      amrex::EB2::translate(
          amrex::EB2::rotate(fub::amrex::Geometry(fub::ConeIF(
                                 {0.0, 0.0, 0.0}, po.plenum_outlet_exit_radius,
                                 po.plenum_outlet_exit_length, true)),
                             M_PI, 2),
          {po.plenum_length + po.plenum_outlet_length, 0.0, 0.0}));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::amrex::CartesianGridGeometry grid_geometry(
      fub::GetOptions(options, "GridGeometry"));

  PatchHierarchyOptions hierarchy_options =
      fub::GetOptions(options, "PatchHierarchy");
  hierarchy_options.index_spaces =
      MakeIndexSpaces(shop, grid_geometry, hierarchy_options);

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};
  using Complete = fub::Complete<fub::IdealGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&Complete::pressure, 0.05},
                             std::pair{&Complete::density, 0.01}};

  fub::Complete<fub::IdealGasMix<Plenum_Rank>> state(equation);
  {
    using namespace std::literals;
    const fub::ProgramOptions initial_options =
        fub::GetOptions(options, "InitialCondition");
    const std::string moles =
        fub::GetOptionOr(initial_options, "moles", "N2:79,O2:21"s);
    const double temperature =
        fub::GetOptionOr(initial_options, "temperature", 300.0);
    const double pressure =
        fub::GetOptionOr(initial_options, "pressure", 101325.0);
    equation.GetReactor().SetMoleFractions(moles);
    equation.GetReactor().SetTemperature(temperature);
    equation.GetReactor().SetPressure(pressure);
    equation.CompleteFromReactor(state);
  }

  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04),
                              state, state);

  ::amrex::RealBox xbox = grid_geometry.coordinates;
  ::amrex::Geometry coarse_geom = fub::amrex::GetCoarseGeometry(grid_geometry);

  ::amrex::RealBox inlet{{xbox.lo(0), -r_tube, -r_tube},
                         {0.01, +r_tube, +r_tube}};
  ::amrex::Box refine_box = fub::amrex::BoxWhichContains(inlet, coarse_geom);
  ConstantBox constant_in_box{refine_box};

  ::amrex::RealBox outlet{{0.5, xbox.lo(1), xbox.lo(2)},
                          {xbox.hi(0), xbox.hi(1), xbox.hi(2)}};
  refine_box = fub::amrex::BoxWhichContains(outlet, coarse_geom);
  ConstantBox constant_out_box{refine_box};

  IsentropicPressureBoundaryOptions boundary_options =
      fub::GetOptions(options, "IsentropicPressureBoundary");
  BoundarySet boundary_condition{
      {TransmissiveBoundary{fub::Direction::X, 0},
       IsentropicPressureBoundary{equation, boundary_options}}};

  std::shared_ptr gridding = [&] {
    std::string checkpoint =
        fub::GetOptionOr(options, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(equation, grid_geometry, hierarchy_options),
          initial_data,
          TagAllOf(TagCutCells(), gradients, constant_in_box, constant_out_box,
                   TagBuffer(2)),
          boundary_condition);
      grid->InitializeHierarchy(0.0);
      return grid;
    } else {
      checkpoint += "/Plenum";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), grid_geometry,
          hierarchy_options);
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), gradients, constant_in_box, constant_out_box,
                   TagBuffer(2)),
          boundary_condition);
      return grid;
    }
  }();

  // Make Solver

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::ideal_gas::MusclHancockPrimMethod<Plenum_Rank> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  const int scratch_gcw = 4;
  const int flux_gcw = 2;

  return IntegratorContext(gridding, method, scratch_gcw, flux_gcw);
}

void MyMain(const std::map<std::string, pybind11::object>& vm);

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

void WriteCheckpoint(const std::string& path,
                     const fub::amrex::MultiBlockGriddingAlgorithm& grid,
                     std::shared_ptr<fub::amrex::PressureValve> valve, int rank,
                     const fub::amrex::MultiBlockIgniteDetonation& ignition) {
  auto tubes = grid.GetTubes();
  std::string name = fmt::format("{}/Tube", path);
  fub::amrex::WriteCheckpointFile(name, tubes[0]->GetPatchHierarchy());
  if (rank == 0) {
    name = fmt::format("{}/Valve", path);
    std::ofstream valve_checkpoint(name);
    boost::archive::text_oarchive oa(valve_checkpoint);
    oa << *valve;
  }
  name = fmt::format("{}/Plenum", path);
  fub::amrex::cutcell::WriteCheckpointFile(
      name, grid.GetPlena()[0]->GetPatchHierarchy());
  if (rank == 0) {
    name = fmt::format("{}/Ignition", path);
    std::ofstream ignition_checkpoint(name);
    boost::archive::text_oarchive oa(ignition_checkpoint);
    oa << ignition.GetLastIgnitionTimePoints();
  }
}

void MyMain(const std::map<std::string, pybind11::object>& vm) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};

  fub::Burke2012 mechanism{};

  std::vector<fub::amrex::cutcell::IntegratorContext> plenum{};
  std::vector<fub::amrex::IntegratorContext> tubes{};
  std::vector<fub::amrex::BlockConnection> connectivity{};
  std::vector<std::shared_ptr<fub::amrex::PressureValve>> valves{};

  plenum.push_back(MakePlenumSolver(mechanism, fub::GetOptions(vm, "Plenum")));
  auto counter_database = plenum[0].GetCounterRegistry();

  auto&& [tube, valve] =
      MakeTubeSolver(mechanism, fub::GetOptions(vm, "Tube"), counter_database);
  fub::amrex::BlockConnection connection;
  connection.direction = fub::Direction::X;
  connection.side = 0;
  connection.ghost_cell_width = 4;
  connection.plenum.id = 0;
  connection.tube.id = 0;
  connection.tube.mirror_box =
      tube.GetGriddingAlgorithm()->GetPatchHierarchy().GetGeometry(0).Domain();
  connection.plenum.mirror_box = plenum[0]
                                     .GetGriddingAlgorithm()
                                     ->GetPatchHierarchy()
                                     .GetGeometry(0)
                                     .Domain();
  tubes.push_back(std::move(tube));
  valves.push_back(valve.GetSharedState());
  connectivity.push_back(connection);

  fub::IdealGasMix<Plenum_Rank> plenum_equation{mechanism};
  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};

  fub::amrex::MultiBlockIntegratorContext context(
      fub::FlameMasterReactor(mechanism), std::move(tubes), std::move(plenum),
      std::move(connectivity));

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<Plenum_Rank>, std::move(context), fub::GodunovSplitting{});

  const std::size_t n_tubes = system_solver.GetContext().Tubes().size();
  const int max_number_of_levels = system_solver.GetContext()
                                       .Tubes()[0]
                                       .GetPatchHierarchy()
                                       .GetMaxNumberOfLevels();

  fub::amrex::MultiBlockIgniteDetonation ignition{
      tube_equation, n_tubes, max_number_of_levels,
      fub::amrex::IgniteDetonationOptions(vm, "IgniteDetonation")};

  std::string checkpoint = fub::GetOptionOr(vm, "checkpoint", std::string{});
  if (!checkpoint.empty()) {
    MPI_Comm comm = context.GetMpiCommunicator();
    std::string input =
        fub::ReadAndBroadcastFile(checkpoint + "/Ignition", comm);
    std::istringstream ifs(input);
    {
      boost::archive::text_iarchive ia(ifs);
      std::vector<fub::Duration> last_ignitions;
      ia >> last_ignitions;
      ignition.SetLastIgnitionTimePoints(last_ignitions);
    }
    input = fub::ReadAndBroadcastFile(checkpoint + "/Valve", comm);
    ifs = std::istringstream(input);
    boost::archive::text_iarchive ia(ifs);
    ia >> *valve.GetSharedState();
  }
  fub::SplitSystemSourceLevelIntegrator ign_solver(
      std::move(system_solver), std::move(ignition), fub::GodunovSplitting{});

  fub::amrex::MultiBlockKineticSouceTerm source_term(tube_equation);

  fub::SplitSystemSourceLevelIntegrator level_integrator(
      std::move(ign_solver), std::move(source_term), fub::StrangSplitting{});

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  using namespace fub::amrex;
  struct MakeCheckpoint
      : public fub::OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm> {
    std::string directory_ = "SingleTubePlenum";
    std::shared_ptr<fub::amrex::PressureValve> valve_state_{};
    const fub::amrex::MultiBlockIgniteDetonation* ignition_{};
    MakeCheckpoint(
        const fub::ProgramOptions& options,
        const std::shared_ptr<fub::amrex::PressureValve>& valve_state,
        const fub::amrex::MultiBlockIgniteDetonation* ignite)
        : OutputAtFrequencyOrInterval(options),
          valve_state_(valve_state), ignition_{ignite} {
      directory_ = fub::GetOptionOr(options, "directory", directory_);
    }
    void operator()(const MultiBlockGriddingAlgorithm& grid) override {
      int rank = -1;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::string name =
          fmt::format("{}/Checkpoint/{:05}", directory_, grid.GetCycles());
      amrex::Print() << "Write Checkpoint to '" << name << "'!\n";
      WriteCheckpoint(name, grid, valve_state_, rank, *ignition_);
    }
  };

  fub::OutputFactory<MultiBlockGriddingAlgorithm> factory{};
  factory.RegisterOutput<MultiWriteHdf5>("HDF5");
  factory.RegisterOutput<MultiBlockPlotfileOutput>("Plotfiles");
  factory.RegisterOutput<LogProbesOutput>("LogProbes");
  factory.RegisterOutput<MakeCheckpoint>(
      "Checkpoint", valve.GetSharedState(),
      &solver.GetLevelIntegrator().GetSystem().GetSource());
  using CounterOutput =
      fub::CounterOutput<fub::amrex::MultiBlockGriddingAlgorithm,
                         std::chrono::milliseconds>;
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  fub::MultipleOutputs<MultiBlockGriddingAlgorithm> outputs(
      std::move(factory), fub::GetOptions(vm, "Output"));

  outputs(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, fub::GetOptions(vm, "RunOptions"),
                     wall_time_reference, outputs);
}
