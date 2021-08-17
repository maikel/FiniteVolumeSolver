// Copyright (c) 2021 Christian Zenker
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

#include "fub/AMReX/cutcell/FluxMethodFactory.hpp"
#include "fub/flux_method/MusclHancockMethod2.hpp"

#include "fub/AMReX/boundary_condition/MassflowBoundary_PerfectGas.hpp"
#include "fub/AMReX/boundary_condition/ShockValveBoundary.hpp"
#include "fub/AMReX/cutcell/AxiSymmetricSourceTerm_PerfectGas.hpp"
#include "fub/equations/perfect_gas/InitializeShock_CutCell_Multiblock.hpp"

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

auto MakeTubeSolver(const fub::ProgramOptions& options,
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

  using Complete = fub::PerfectGas<Tube_Rank>::Complete;
  fub::PerfectGas<Tube_Rank> tube_equation;

  fub::ProgramOptions equation_options = fub::GetOptions(options, "Equation");
  tube_equation.gamma =
      fub::GetOptionOr(equation_options, "gamma", tube_equation.gamma);
  tube_equation.Rspec =
      fub::GetOptionOr(equation_options, "Rpsec", tube_equation.Rspec);
  tube_equation.gamma_minus_1_inv = 1.0 / (tube_equation.gamma - 1.0);
  tube_equation.gamma_array_ = fub::Array1d::Constant(tube_equation.gamma);
  tube_equation.gamma_minus_1_inv_array_ =
      fub::Array1d::Constant(tube_equation.gamma_minus_1_inv);

  BOOST_LOG(log) << "Equation:";
  BOOST_LOG(log) << fmt::format(" - gamma = {}", tube_equation.gamma);
  BOOST_LOG(log) << fmt::format(" - Rspec = {}", tube_equation.Rspec);

  GradientDetector gradient{tube_equation,
                            std::make_pair(&Complete::density, 1e-3),
                            std::make_pair(&Complete::pressure, 1e-2)};

  DataDescription desc = MakeDataDescription(tube_equation);

  ::amrex::Box refine_box{{grid_geometry.cell_dimensions[0] - 5, 0},
                          {grid_geometry.cell_dimensions[0] - 1, 0}};
  ConstantBox constant_box{refine_box};

  fub::Complete<fub::PerfectGas<Tube_Rank>> initial_state;
  {
    using namespace std::literals;
    fub::Primitive<fub::PerfectGas<Tube_Rank>> prim;
    const fub::ProgramOptions initial_options =
        fub::GetOptions(options, "InitialCondition");
    const double density = fub::GetOptionOr(initial_options, "density", 1.22);
    const double u_velocity =
        fub::GetOptionOr(initial_options, "u_velocity", 0.0);
    [[maybe_unused]] const double v_velocity =
        fub::GetOptionOr(initial_options, "v_velocity", 0.0);
    [[maybe_unused]] const double w_velocity =
        fub::GetOptionOr(initial_options, "w_velocity", 0.0);
    const double pressure =
        fub::GetOptionOr(initial_options, "pressure", 101325.0);
    prim.density = density;
    prim.velocity << u_velocity; //AMREX_D_DECL(u_velocity, v_velocity, w_velocity);
    prim.pressure = pressure;
    fub::CompleteFromPrim(tube_equation, initial_state, prim);
  }

  ConstantData<fub::PerfectGas<Tube_Rank>> initial_data{tube_equation,
                                                        initial_state};

  ShockValveOptions valve_options =
      fub::GetOptions(options, "ShockValveBoundary");
  ShockValveBoundary valve{tube_equation, valve_options};
  BOOST_LOG(log) << "ShockValveBoundary:";
  valve_options.Print(log);

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
          TagAllOf(gradient, constant_box), valve);
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
                                              valve);
      return gridding;
    }
  }();

  auto [flux_method, time_integrator] =
      fub::amrex::GetFluxMethod(fub::GetOptions(options, "FluxMethod"),
                                gridding->GetPatchHierarchy(), tube_equation);
  HyperbolicMethod method{flux_method, time_integrator,
                          Reconstruction(tube_equation)};

  IntegratorContext context(gridding, method,
                            fub::GetOptions(options, "IntegratorContext"));
  BOOST_LOG(log) << "IntegratorContext:";
  context.GetOptions().Print(log);

  BOOST_LOG(log) << "==================== End Tube =========================";

  return context;
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

auto MakePlenumSolver(const std::map<std::string, pybind11::object>& options) {
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
                 [&log](const std::string& filename) {
                   BOOST_LOG(log) << "Reading... " << filename;
                   std::ifstream ifs(filename);
                   return fub::PolymorphicGeometry<2>(ReadPolygonData(ifs));
                 });
  BOOST_LOG(log) << "Compute EB level set data...";
  auto embedded_boundary =
      fub::amrex::Geometry(fub::PolymorphicUnion(geometries));
  auto shop = amrex::EB2::makeShop(embedded_boundary);
  hierarchy_options.index_spaces = fub::amrex::cutcell::MakeIndexSpaces(
      shop, grid_geometry, hierarchy_options);
  std::shared_ptr<fub::CounterRegistry> registry =
      std::make_shared<fub::CounterRegistry>();
  {
    fub::Timer timer = registry->get_timer("MakeIndexSpaces");
    hierarchy_options.index_spaces =
        MakeIndexSpaces(shop, grid_geometry, hierarchy_options);
  }

  fub::PerfectGas<Plenum_Rank> plenum_equation;

  fub::ProgramOptions equation_options = fub::GetOptions(options, "Equation");
  plenum_equation.gamma =
      fub::GetOptionOr(equation_options, "gamma", plenum_equation.gamma);
  plenum_equation.Rspec =
      fub::GetOptionOr(equation_options, "Rpsec", plenum_equation.Rspec);
  plenum_equation.gamma_minus_1_inv = 1.0 / (plenum_equation.gamma - 1.0);
  plenum_equation.gamma_array_ = fub::Array1d::Constant(plenum_equation.gamma);
  plenum_equation.gamma_minus_1_inv_array_ =
      fub::Array1d::Constant(plenum_equation.gamma_minus_1_inv);

  fub::Complete<fub::PerfectGas<2>> initial_state;
  {
    using namespace std::literals;
    fub::Primitive<fub::PerfectGas<2>> prim;
    const fub::ProgramOptions initial_options =
        fub::GetOptions(options, "InitialCondition");
    const double density = fub::GetOptionOr(initial_options, "density", 1.22);
    const double u_velocity =
        fub::GetOptionOr(initial_options, "u_velocity", 0.0);
    [[maybe_unused]] const double v_velocity =
        fub::GetOptionOr(initial_options, "v_velocity", 0.0);
    [[maybe_unused]] const double w_velocity =
        fub::GetOptionOr(initial_options, "w_velocity", 0.0);
    const double pressure =
        fub::GetOptionOr(initial_options, "pressure", 101325.0);
    prim.density = density;
    prim.velocity << AMREX_D_DECL(u_velocity, v_velocity, w_velocity);
    prim.pressure = pressure;
    fub::CompleteFromPrim(plenum_equation, initial_state, prim);
  }

  RiemannProblem<fub::PerfectGas<2>, fub::Halfspace> initial_data{
      plenum_equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.0), initial_state,
      initial_state};

  using Complete = fub::Complete<fub::PerfectGas<Plenum_Rank>>;
  GradientDetector gradients{plenum_equation,
                             std::pair{&Complete::pressure, 0.05},
                             std::pair{&Complete::density, 0.05}};

  auto seq = fub::execution::seq;
  BoundarySet boundary_condition{
      {TransmissiveBoundary{fub::Direction::X, 0},
       TransmissiveBoundary{fub::Direction::X, 1},
       ReflectiveBoundary{seq, plenum_equation, fub::Direction::Y, 0},
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
          PatchHierarchy(plenum_equation, grid_geometry, hierarchy_options),
          initial_data, TagAllOf(TagCutCells(), constant_refinebox),
          boundary_condition);
      grid->InitializeHierarchy(0.0);
      return grid;
    } else {
      checkpoint += "/Plenum";
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(plenum_equation),
          grid_geometry, hierarchy_options);
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), constant_refinebox), boundary_condition);
      return grid;
    }
  }();

  using namespace std::literals;

  auto flux_method = fub::amrex::cutcell::GetCutCellMethod(
      fub::GetOptions(options, "FluxMethod"), plenum_equation);

  HyperbolicMethod method{flux_method, TimeIntegrator{},
                          Reconstruction{plenum_equation}};

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
  // fub::InitializeLogging(MPI_COMM_WORLD);
  pybind11::scoped_interpreter interpreter{};
  std::optional<fub::ProgramOptions> opts = fub::ParseCommandLine(argc, argv);
  if (opts) {
    fub::InitializeLogging(MPI_COMM_WORLD,
                           fub::GetOptions(*opts, "LogOptions"));
    MyMain(*opts);
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}

const fub::amrex::ShockValveBoundary&
GetValveBoundary(const fub::amrex::MultiBlockGriddingAlgorithm2& grid) {
  auto tubes = grid.GetTubes();
  auto boundaries =
      tubes[0]->GetBoundaryCondition().Cast<const fub::amrex::BoundarySet>();
  FUB_ASSERT(boundaries);
  auto valve =
      boundaries->conditions[0].Cast<const fub::amrex::ShockValveBoundary>();
  FUB_ASSERT(valve);
  return *valve;
}

fub::amrex::ShockValveBoundary&
GetValveBoundary(fub::amrex::MultiBlockGriddingAlgorithm2& grid) {
  auto tubes = grid.GetTubes();
  auto boundaries =
      tubes[0]->GetBoundaryCondition().Cast<fub::amrex::BoundarySet>();
  FUB_ASSERT(boundaries);
  auto valve = boundaries->conditions[0].Cast<fub::amrex::ShockValveBoundary>();
  FUB_ASSERT(valve);
  return *valve;
}

void WriteCheckpoint(const std::string& path,
                     const fub::amrex::MultiBlockGriddingAlgorithm2& grid,
                     int rank) {
  auto tubes = grid.GetTubes();
  std::string name = fmt::format("{}/Tube", path);
  fub::amrex::WriteCheckpointFile(name, tubes[0]->GetPatchHierarchy());
  if (rank == 0) {
    const fub::amrex::ShockValveBoundary& valve = GetValveBoundary(grid);
    name = fmt::format("{}/Valve", path);
    std::ofstream valve_checkpoint(name);
    boost::archive::text_oarchive oa(valve_checkpoint);
    oa << valve;
  }
  name = fmt::format("{}/Plenum", path);
  fub::amrex::cutcell::WriteCheckpointFile(
      name, grid.GetPlena()[0]->GetPatchHierarchy());
}

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};

  std::vector<fub::amrex::cutcell::IntegratorContext> plenum{};
  std::vector<fub::amrex::IntegratorContext> tubes{};
  std::vector<fub::amrex::BlockConnection> connectivity{};

  plenum.push_back(MakePlenumSolver(fub::GetOptions(options, "Plenum")));
  auto counter_database = plenum[0].GetCounterRegistry();

  fub::amrex::IntegratorContext tube =
      MakeTubeSolver(fub::GetOptions(options, "Tube"), counter_database);
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
  connectivity.push_back(connection);

  fub::PerfectGas<Plenum_Rank> plenum_equation;
  fub::PerfectGas<Tube_Rank> tube_equation;

  fub::amrex::MultiBlockIntegratorContext2 context(
      tube_equation, plenum_equation, std::move(tubes), std::move(plenum),
      std::move(connectivity));

  fub::amrex::cutcell::feedback_functions::ShockOptions shock_options =
      fub::GetOptions(options, "schock_feedback");
  fub::amrex::cutcell::feedback_functions::ShockFeedback shock_feedback(
      tube_equation, plenum_equation, shock_options);
  context.SetPostAdvanceHierarchyFeedback(shock_feedback);

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<Plenum_Rank>, std::move(context), fub::GodunovSplitting{});

  std::string checkpoint =
      fub::GetOptionOr(options, "checkpoint", std::string{});
  if (!checkpoint.empty()) {
    MPI_Comm comm = context.GetMpiCommunicator();
    fub::amrex::ShockValveBoundary& valve =
        GetValveBoundary(*system_solver.GetGriddingAlgorithm());
    std::string input = fub::ReadAndBroadcastFile(checkpoint + "/Valve", comm);
    std::istringstream ifs(input);
    boost::archive::text_iarchive ia(ifs);
    ia >> valve;
  }

  using namespace fub::amrex;
  struct MakeCheckpoint
      : public fub::OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm2> {
    std::string directory_ = "Divider2D/Checkpoint/";
    MakeCheckpoint(const fub::ProgramOptions& options)
        : OutputAtFrequencyOrInterval(options) {
      directory_ = fub::GetOptionOr(options, "directory", directory_);
      fub::SeverityLogger log = fub::GetInfoLogger();
      BOOST_LOG(log) << "Checkpoint Output configured:";
      BOOST_LOG(log) << fmt::format("  - directory = {}", directory_);
      OutputAtFrequencyOrInterval::Print(log);
    }
    void operator()(const MultiBlockGriddingAlgorithm2& grid) override {
      int rank = -1;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::string name = fmt::format("{}/{:09}", directory_, grid.GetCycles());
      fub::SeverityLogger log = fub::GetInfoLogger();
      BOOST_LOG(log) << fmt::format("Write Checkpoint to '{}'!", name);
      WriteCheckpoint(name, grid, rank);
    }
  };

  using MultiBlockPlotfileOutput =
      MultiBlockPlotfileOutput2<fub::PerfectGas<1>, fub::PerfectGas<2>>;
  fub::OutputFactory<MultiBlockGriddingAlgorithm2> factory{};
  factory.RegisterOutput<MultiWriteHdf5WithNames>("HDF5", plenum_equation,
                                                  tube_equation);
  factory.RegisterOutput<MultiBlockPlotfileOutput>("Plotfiles", tube_equation,
                                                   plenum_equation);
  factory.RegisterOutput<LogProbesOutput>("LogProbes");
  factory.RegisterOutput<MakeCheckpoint>("Checkpoint");
  using CounterOutput =
      fub::CounterOutput<fub::amrex::MultiBlockGriddingAlgorithm2,
                         std::chrono::milliseconds>;
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  fub::MultipleOutputs<MultiBlockGriddingAlgorithm2> outputs(
      std::move(factory), fub::GetOptions(options, "Output"));

  fub::SeverityLogger log = fub::GetInfoLogger();
  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);

  fub::ProgramOptions prog_run_options = fub::GetOptions(options, "RunOptions");
  int AxiSymmetric = fub::GetOptionOr(prog_run_options, "AxiSymmetric", 0);

  if (AxiSymmetric) {
    fub::amrex::AxiSymmetricSourceTerm_PerfectGas symmetry_source_term(
        plenum_equation);

    auto advance_level =
        [](fub::amrex::AxiSymmetricSourceTerm_PerfectGas& source_term,
           fub::amrex::MultiBlockIntegratorContext2& multi_context, int level,
           fub::Duration dt, const amrex::IntVect& ngrow) {
          fub::amrex::cutcell::IntegratorContext& plenum =
              multi_context.Plena()[0];
          return source_term.AdvanceLevel(plenum, level, dt, ngrow);
        };
    auto compute_stable_dt =
        [](fub::amrex::AxiSymmetricSourceTerm_PerfectGas& source_term,
           const fub::amrex::MultiBlockIntegratorContext2& multi_context,
           int level) {
          const fub::amrex::cutcell::IntegratorContext& plenum =
              multi_context.Plena()[0];
          return source_term.ComputeStableDt(plenum, level);
        };
    fub::amrex::MultiBlockLevelIntegrator multi_symmetric_source_term(
        advance_level, compute_stable_dt, symmetry_source_term);

    fub::SplitSystemSourceLevelIntegrator symmetric_level_integrator(
        std::move(system_solver), std::move(multi_symmetric_source_term),
        fub::GodunovSplitting());

    fub::SubcycleFineFirstSolver solver(std::move(symmetric_level_integrator));
    // fub::NoSubcycleSolver solver(std::move(level_integrator));

    outputs(*solver.GetGriddingAlgorithm());

    fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
  } else {
    fub::SubcycleFineFirstSolver solver(std::move(system_solver));
    // fub::NoSubcycleSolver solver(std::move(level_integrator));

    outputs(*solver.GetGriddingAlgorithm());

    fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
  }
}
