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

#include "fub/AMReX/initial_data/PythonData.hpp"
#include "fub/AMReX/cutcell/initial_data/PythonData.hpp"

#include "fub/AMReX/cutcell/FluxMethodFactory.hpp"
#include "fub/flux_method/MusclHancockMethod2.hpp"

#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>

#include <iostream>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>

#include <boost/container/pmr/vector.hpp>
#include <boost/filesystem.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <xmmintrin.h>

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 2;

auto MakeTubeSolver(fub::Burke2012& mechanism,
                    const fub::ProgramOptions& options,
                    const std::shared_ptr<fub::CounterRegistry>& counters) {
  using namespace fub::amrex;

  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << "==================== Tube =========================";

  CartesianGridGeometry grid_geometry(fub::GetOptions(options, "GridGeometry"));
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  PatchHierarchyOptions hierarchy_options(
      fub::GetOptions(options, "PatchHierarchy"));
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  using Complete = fub::IdealGasMix<Tube_Rank>::Complete;
  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 1e-3),
                            std::make_pair(&Complete::pressure, 1e-2),
                            std::make_pair(&Complete::temperature, 1e-1)};

  DataDescription desc = MakeDataDescription(equation);

  ::amrex::Box refine_box{{grid_geometry.cell_dimensions[0] - 5, 0},
                          {grid_geometry.cell_dimensions[0] - 1, 0}};
  ConstantBox constant_box{refine_box};

  PythonData initial_data{equation, fub::GetOptions(options, "InitialData")};

  PressureValveOptions valve_options =
      fub::GetOptions(options, "PressureValveBoundary");
  PressureValveBoundary valve{equation, valve_options};
  BOOST_LOG(log) << "PressureValveBoundary:";
  valve_options.Print(log);

  BoundarySet boundaries{{valve}};

  // If a checkpoint path is specified we will fill the patch hierarchy with
  // data from the checkpoint file, otherwise we will initialize the data by
  // the initial data function.
  std::shared_ptr<GriddingAlgorithm> gridding = [&] {
    std::string checkpoint =
        fub::GetOptionOr(options, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      BOOST_LOG(log) << "Initialize grid by initial condition...";
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(desc, grid_geometry, hierarchy_options), initial_data,
          TagAllOf(gradient, constant_box), boundaries);
      gridding->GetPatchHierarchy().SetCounterRegistry(counters);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint = fmt::format("{}/Tube", checkpoint);
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
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
  HyperbolicMethod method{FluxMethodAdapter(flux_method),
                          EulerForwardTimeIntegrator(),
                          Reconstruction(equation)};

  IntegratorContext context(gridding, method,
                            fub::GetOptions(options, "IntegratorContext"));
  BOOST_LOG(log) << "IntegratorContext:";
  context.GetOptions().Print(log);

  BOOST_LOG(log) << "==================== End Tube =========================";

  return std::pair{context, valve};
}

fub::Polygon ReadPolygonData(std::istream& input) {
  std::string line{};
  namespace pmr = boost::container::pmr;
  pmr::vector<double> xs{};
  pmr::vector<double> ys{};
  while (std::getline(input, line)) {
    double x{}, y{};
    std::istringstream linestream(line);
    linestream >> x >> y;
    if (linestream) {
      xs.push_back(x);
      ys.push_back(y);
    }
  }
  if (xs.empty() || ys.empty()) {
    throw std::invalid_argument{"Invalid Input File:: No data found!"};
  }
  const double x0 = xs.front();
  const double x1 = xs.back();
  const double y0 = ys.front();
  const double y1 = ys.back();
  if (x0 != x1 || y0 != y1) {
    throw std::invalid_argument{
        "Invalid Input File: First and last entries are not the same point."};
  }
  return fub::Polygon(std::move(xs), std::move(ys));
}

auto MakePlenumSolver(fub::Burke2012& mechanism,
                      const std::map<std::string, pybind11::object>& options) {
  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << "==================== Plenum =========================";
  using namespace fub::amrex::cutcell;

  fub::amrex::CartesianGridGeometry grid_geometry =
      fub::GetOptions(options, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  PatchHierarchyOptions hierarchy_options =
      fub::GetOptions(options, "PatchHierarchy");
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  BOOST_LOG(log) << "Read Wall Data Files...";
  std::vector<std::string> wall_filenames{};
  wall_filenames = fub::GetOptionOr(options, "wall_filenames", wall_filenames);

  std::vector<fub::PolymorphicGeometry<2>> geometries;
  std::transform(wall_filenames.begin(), wall_filenames.end(),
                 std::back_inserter(geometries),
                 [](const std::string& filename) {
                   std::ifstream ifs(filename);
                   return fub::PolymorphicGeometry<2>(ReadPolygonData(ifs));
                 });
  BOOST_LOG(log) << "Compute EB level set data...";
  auto embedded_boundary =
      fub::amrex::Geometry(fub::PolymorphicUnion(geometries));
  auto shop = amrex::EB2::makeShop(embedded_boundary);
  hierarchy_options.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, grid_geometry, hierarchy_options);
  std::shared_ptr<fub::CounterRegistry> registry =
      std::make_shared<fub::CounterRegistry>();
  {
    fub::Timer timer = registry->get_timer("MakeIndexSpaces");
    hierarchy_options.index_spaces =
        MakeIndexSpaces(shop, grid_geometry, hierarchy_options);
  }

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};

  PythonData initial_data{equation, fub::GetOptions(options, "InitialData")};

  //  using Complete = fub::Complete<fub::PerfectGas<Plenum_Rank>>;
  using Complete = fub::Complete<fub::IdealGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&Complete::pressure, 0.05},
                             std::pair{&Complete::density, 0.05}};

  auto seq = fub::execution::seq;
  BoundarySet boundary_condition{
      {TransmissiveBoundary{fub::Direction::X, 0},
       TransmissiveBoundary{fub::Direction::X, 1},
       ReflectiveBoundary{seq, equation, fub::Direction::Y, 0},
       TransmissiveBoundary{fub::Direction::Y, 1}}};

  ::amrex::RealBox xbox = grid_geometry.coordinates;
  ::amrex::Geometry coarse_geom = fub::amrex::GetCoarseGeometry(grid_geometry);

  ::amrex::RealBox inlet{{xbox.lo(0), -0.15}, {0.18, +0.15}};
  ::amrex::Box refine_box = fub::amrex::BoxWhichContains(inlet, coarse_geom);
  ConstantBox constant_refinebox{refine_box};

  std::shared_ptr gridding = [&] {
    std::string checkpoint =
        fub::GetOptionOr(options, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      BOOST_LOG(log) << "Initialize grid by initial condition...";
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(equation, grid_geometry, hierarchy_options),
          initial_data, TagAllOf(TagCutCells(), constant_refinebox),
          boundary_condition);
      grid->InitializeHierarchy(0.0);
      return grid;
    } else {
      checkpoint += "/Plenum";
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), grid_geometry,
          hierarchy_options);
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), constant_refinebox), boundary_condition);
      return grid;
    }
  }();

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::ideal_gas::MusclHancockPrimMethod<Plenum_Rank> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  IntegratorContext context(gridding, method,
                            fub::GetOptions(options, "IntegratorContext"));
  BOOST_LOG(log) << "IntegratorContext:";
  context.GetOptions().Print(log);

  BOOST_LOG(log) << "==================== End Plenum =========================";

  return context;
}

void MyMain(const fub::ProgramOptions& options);

int main(int argc, char** argv) {
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  int provided{-1};
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
  if (provided < MPI_THREAD_FUNNELED) {
    fmt::print(
        "Aborting execution. MPI could not provide a thread-safe instance.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
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
                     const fub::amrex::MultiBlockGriddingAlgorithm2& grid,
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
    oa << ignition.GetNextIgnitionTimePoints();
  }
}

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};

  fub::Burke2012 mechanism{};

  std::vector<fub::amrex::cutcell::IntegratorContext> plenum{};
  std::vector<fub::amrex::IntegratorContext> tubes{};
  std::vector<fub::amrex::BlockConnection> connectivity{};
  std::vector<std::shared_ptr<fub::amrex::PressureValve>> valves{};

  plenum.push_back(
      MakePlenumSolver(mechanism, fub::GetOptions(options, "Plenum")));
  auto counter_database = plenum[0].GetCounterRegistry();

  auto&& [tube, valve] = MakeTubeSolver(
      mechanism, fub::GetOptions(options, "Tube"), counter_database);
  fub::amrex::BlockConnection connection;
  connection.direction = fub::Direction::X;
  connection.side = 0;
  connection.ghost_cell_width = std::max(tube.GetOptions().scratch_gcw,
                                         plenum[0].GetOptions().scratch_gcw);
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

  fub::amrex::MultiBlockIntegratorContext2 context(
      tube_equation, plenum_equation, std::move(tubes), std::move(plenum),
      std::move(connectivity));

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<Plenum_Rank>, std::move(context), fub::GodunovSplitting{});

  const std::size_t n_tubes = system_solver.GetContext().Tubes().size();
  const int max_number_of_levels = system_solver.GetContext()
                                       .Tubes()[0]
                                       .GetPatchHierarchy()
                                       .GetMaxNumberOfLevels();

  fub::amrex::IgniteDetonationOptions ignite_options(
      fub::GetOptions(options, "IgniteDetonation"));
  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << "IgniteDetonation:";
  ignite_options.Print(log);
  fub::amrex::MultiBlockIgniteDetonation ignition{
      tube_equation, n_tubes, max_number_of_levels, ignite_options};

  std::string checkpoint =
      fub::GetOptionOr(options, "checkpoint", std::string{});
  if (!checkpoint.empty()) {
    MPI_Comm comm = context.GetMpiCommunicator();
    std::string input =
        fub::ReadAndBroadcastFile(checkpoint + "/Ignition", comm);
    std::istringstream ifs(input);
    {
      boost::archive::text_iarchive ia(ifs);
      std::vector<fub::Duration> last_ignitions;
      ia >> last_ignitions;
      ignition.SetNextIgnitionTimePoints(last_ignitions);
    }
    input = fub::ReadAndBroadcastFile(checkpoint + "/Valve", comm);
    ifs = std::istringstream(input);
    boost::archive::text_iarchive ia(ifs);
    ia >> *valve.GetSharedState();
  }
  fub::SplitSystemSourceLevelIntegrator ign_solver(
      std::move(system_solver), std::move(ignition), fub::GodunovSplitting{});

  fub::amrex::MultiBlockKineticSouceTerm source_term(tube_equation);

  fub::amrex::AxiSymmetricSourceTerm symmetry_source_term(plenum_equation);

  auto advance_level =
      [](fub::amrex::AxiSymmetricSourceTerm& source_term,
         fub::amrex::MultiBlockIntegratorContext2& multi_context, int level,
         fub::Duration dt, const amrex::IntVect& ngrow) {
        fub::amrex::cutcell::IntegratorContext& plenum =
            multi_context.Plena()[0];
        return source_term.AdvanceLevel(plenum, level, dt, ngrow);
      };
  auto compute_stable_dt =
      [](fub::amrex::AxiSymmetricSourceTerm& source_term,
         const fub::amrex::MultiBlockIntegratorContext2& multi_context,
         int level) {
        const fub::amrex::cutcell::IntegratorContext& plenum =
            multi_context.Plena()[0];
        return source_term.ComputeStableDt(plenum, level);
      };
  fub::amrex::MultiBlockLevelIntegrator multi_symmetric_source_term(
      advance_level, compute_stable_dt, symmetry_source_term);

  fub::SplitSystemSourceLevelIntegrator symmetric_level_integrator(
      std::move(ign_solver), std::move(multi_symmetric_source_term),
      fub::GodunovSplitting());

  fub::SplitSystemSourceLevelIntegrator level_integrator(
      std::move(symmetric_level_integrator), std::move(source_term),
      fub::StrangSplittingLumped{});

  // fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
  fub::NoSubcycleSolver solver(std::move(level_integrator));

  using namespace fub::amrex;
  struct MakeCheckpoint
      : public fub::OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm2> {
    std::string directory_ = "ConvergentNozzle";
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
    void operator()(const MultiBlockGriddingAlgorithm2& grid) override {
      int rank = -1;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::string name = fmt::format("{}/{:09}", directory_, grid.GetCycles());
      fub::SeverityLogger log = fub::GetInfoLogger();
      BOOST_LOG(log) << fmt::format("Write Checkpoint to '{}'!", name);
      WriteCheckpoint(name, grid, valve_state_, rank, *ignition_);
    }
  };

  using MultiBlockPlotfileOutput =
      MultiBlockPlotfileOutput2<fub::IdealGasMix<1>, fub::IdealGasMix<2>>;
  fub::OutputFactory<MultiBlockGriddingAlgorithm2> factory{};
  factory.RegisterOutput<MultiWriteHdf5WithNames>("HDF5", plenum_equation,
                                                  tube_equation);
  factory.RegisterOutput<MultiBlockPlotfileOutput>("Plotfiles", tube_equation,
                                                   plenum_equation);
  factory.RegisterOutput<LogProbesOutput>("LogProbes");
  factory.RegisterOutput<MakeCheckpoint>(
      "Checkpoint", valve.GetSharedState(),
      &solver.GetLevelIntegrator().GetSystem().GetSystem().GetSource());
  using CounterOutput =
      fub::CounterOutput<fub::amrex::MultiBlockGriddingAlgorithm2,
                         std::chrono::milliseconds>;
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  fub::MultipleOutputs<MultiBlockGriddingAlgorithm2> outputs(
      std::move(factory), fub::GetOptions(options, "Output"));

  outputs(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);
  fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
}
