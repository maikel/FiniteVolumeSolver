// Copyright (c) 2021 Christian Zenker
// Copyright (c) 2021 Maikel Nadolski
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

#include "fub/AMReX/cutcell/AxiSymmetricSourceTerm_PerfectGas.hpp"
#include "fub/AMReX/cutcell/boundary_condition/MassflowBoundary_PerfectGas.hpp"
#include "fub/AMReX/cutcell/boundary_condition/ShockValveBoundary.hpp"
#include "fub/equations/perfect_gas/InitializeShock.hpp"

#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>

#include <iostream>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>

#include <boost/container/pmr/vector.hpp>
#include <boost/filesystem.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <xmmintrin.h>

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

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};
  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::PerfectGas<2> equation;

  fub::ProgramOptions equation_options = fub::GetOptions(options, "Equation");
  equation.gamma = fub::GetOptionOr(equation_options, "gamma", equation.gamma);
  equation.Rspec = fub::GetOptionOr(equation_options, "Rpsec", equation.Rspec);
  equation.gamma_minus_1_inv = 1.0 / (equation.gamma - 1.0);
  equation.gamma_array_ = fub::Array1d::Constant(equation.gamma);
  equation.gamma_minus_1_inv_array_ =
      fub::Array1d::Constant(equation.gamma_minus_1_inv);

  BOOST_LOG(log) << "Equation:";
  BOOST_LOG(log) << fmt::format(" - gamma = {}", equation.gamma);
  BOOST_LOG(log) << fmt::format(" - Rspec = {}", equation.Rspec);

  fub::amrex::CartesianGridGeometry geometry =
      fub::GetOptions(options, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  geometry.Print(log);

  const int scratch_gcw = fub::GetOptionOr(options, "scratch_gcw", 4);
  const int flux_gcw = fub::GetOptionOr(options, "flux_gcw", 2);
  fub::amrex::cutcell::PatchHierarchyOptions hier_opts =
      fub::GetOptions(options, "PatchHierarchy");
  hier_opts.ngrow_eb_level_set =
      std::max(scratch_gcw + 1, hier_opts.ngrow_eb_level_set);
  BOOST_LOG(log) << "PatchHierarchy:";
  hier_opts.Print(log);

  BOOST_LOG(log) << "Read Wall Data...";
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
  hier_opts.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, geometry, hier_opts);

  using namespace fub::amrex::cutcell;
  using State = fub::Complete<fub::PerfectGas<2>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

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
    fub::CompleteFromPrim(equation, initial_state, prim);
  }

  RiemannProblem<fub::PerfectGas<2>, fub::Halfspace> initial_data{
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.0), initial_state,
      initial_state};

  auto seq = fub::execution::seq;
  BoundarySet boundary_condition{
      {TransmissiveBoundary{fub::Direction::X, 0},
       TransmissiveBoundary{fub::Direction::X, 1},
       ReflectiveBoundary{seq, equation, fub::Direction::Y, 0},
       TransmissiveBoundary{fub::Direction::Y, 1}}};

  ShockValveOptions valve_options =
      fub::GetOptions(options, "ShockValveBoundary");
  ShockValveBoundary shock_valve(equation, valve_options);
  BOOST_LOG(log) << "ShockValveBoundary:";
  valve_options.Print(log);
  boundary_condition.conditions.push_back(std::move(shock_valve));

  std::shared_ptr gridding = [&] {
    fub::SeverityLogger log = fub::GetInfoLogger();
    std::string checkpoint =
        fub::GetOptionOr(options, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      BOOST_LOG(log) << "Initialize grid by initial condition...";
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(equation, geometry, hier_opts), initial_data,
          TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
      grid->InitializeHierarchy(0.0);
      return grid;
    } else {
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), geometry,
          hier_opts);
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
      return grid;
    }
  }();

  using namespace std::literals;

  auto flux_method = fub::amrex::cutcell::GetCutCellMethod(
      fub::GetOptions(options, "FluxMethod"), equation);

  HyperbolicMethod method{flux_method, TimeIntegrator{},
                          Reconstruction{equation}};

  BOOST_LOG(log) << fmt::format("scratch_gcw = {}", scratch_gcw);
  BOOST_LOG(log) << fmt::format("flux_gcw = {}", flux_gcw);

  IntegratorContext context(gridding, method, scratch_gcw, flux_gcw);

  fub::amrex::cutcell::feedback_functions::ShockOptions shock_options =
      fub::GetOptions(options, "schock_feedback");
  fub::amrex::cutcell::feedback_functions::ShockFeedback shock_feedback(
      equation, shock_options);
  context.SetFeedbackFunction(shock_feedback);

  using CounterOutput =
      fub::CounterOutput<GriddingAlgorithm, std::chrono::milliseconds>;
  using Plotfile = PlotfileOutput<fub::PerfectGas<2>>;
  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  auto field_names = fub::VarNames<State, std::vector<std::string>>(equation);
  factory.RegisterOutput<WriteHdf5>("HDF5", field_names);
  factory.RegisterOutput<Plotfile>("Plotfiles", equation);

  struct MakeCheckpoint
      : public fub::OutputAtFrequencyOrInterval<GriddingAlgorithm> {
    std::string directory_ = "Divider2D/Checkpoint/";
    MakeCheckpoint(const fub::ProgramOptions& options)
        : OutputAtFrequencyOrInterval(options) {
      directory_ = fub::GetOptionOr(options, "directory", directory_);
      fub::SeverityLogger log = fub::GetInfoLogger();
      BOOST_LOG(log) << "Checkpoint Output configured:";
      BOOST_LOG(log) << fmt::format("  - directory = {}", directory_);
      OutputAtFrequencyOrInterval::Print(log);
    }
    void operator()(const GriddingAlgorithm& grid) override {
      std::string name = fmt::format("{}/{:09}", directory_, grid.GetCycles());
      fub::SeverityLogger log = fub::GetInfoLogger();
      BOOST_LOG(log) << fmt::format("Write Checkpoint to '{}'!", name);
      fub::amrex::cutcell::WriteCheckpointFile(name, grid.GetPatchHierarchy());
    }
  };
  factory.RegisterOutput<MakeCheckpoint>("Checkpoint");
  fub::MultipleOutputs<GriddingAlgorithm> outputs(
      std::move(factory), fub::GetOptions(options, "Output"));

  fub::DimensionalSplitLevelIntegrator level_integrator(fub::int_c<2>,
                                                        std::move(context));

  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);

  fub::ProgramOptions prog_run_options = fub::GetOptions(options, "RunOptions");
  int AxiSymmetric = fub::GetOptionOr(prog_run_options, "AxiSymmetric", 0);

  if (AxiSymmetric) {
    fub::amrex::AxiSymmetricSourceTerm_PerfectGas symmetry_source_term(
        equation);

    fub::SplitSystemSourceLevelIntegrator symmetric_level_integrator(
        std::move(level_integrator), std::move(symmetry_source_term),
        fub::GodunovSplitting());

    fub::NoSubcycleSolver solver(std::move(symmetric_level_integrator));

    outputs(*solver.GetGriddingAlgorithm());
    fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
  } else {
    fub::NoSubcycleSolver solver(std::move(level_integrator));

    outputs(*solver.GetGriddingAlgorithm());
    fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
  }

  // fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
}

int main(int argc, char** argv) {
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  MPI_Init(nullptr, nullptr);
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
