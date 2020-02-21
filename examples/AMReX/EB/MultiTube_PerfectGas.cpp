// Copyright (c) 2020 Maikel Nadolski
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

#include "fub/AMReX/cutcell/output/PerfectGasProbesOutput.hpp"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

#include <cmath>
#include <iostream>

#include <pybind11/numpy.h>

static constexpr int Plenum_Rank = 3;

static constexpr double r_tube = 0.015;
static constexpr double r_inner = 0.5 * 0.130;
static constexpr double r_outer = 0.5 * 0.385;
static constexpr double r_tube_center = 0.5 * r_inner + 0.5 * r_outer;
static constexpr double alpha = 2. * M_PI / 6.;

auto Center(double x, double phi) -> ::amrex::RealArray {
  using std::cos;
  using std::sin;
  return {x, r_tube_center * sin(phi), r_tube_center * cos(phi)};
}

auto MakePlenumSolver(fub::PerfectGas<3>& equation,
                      const fub::ProgramOptions& options) {
  using namespace fub::amrex::cutcell;
  fub::SeverityLogger log = fub::GetInfoLogger();
  fub::amrex::CartesianGridGeometry grid_geometry(
      fub::GetOptions(options, "GridGeometry"));
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  auto Cylinder = [&](double radius, double height,
                      const std::array<double, 3>& center) {
    return amrex::EB2::CylinderIF(radius, height, 0, center, true);
  };

  auto embedded_boundary = amrex::EB2::makeUnion(
      amrex::EB2::makeIntersection(
          Cylinder(r_outer, 1.0, {0.5, 0.0, 0.0}),
          Cylinder(r_tube, 2.1, Center(-1.0, 0.0 * alpha)),
          Cylinder(r_tube, 2.1, Center(-1.0, 1.0 * alpha)),
          Cylinder(r_tube, 2.1, Center(-1.0, 2.0 * alpha)),
          Cylinder(r_tube, 2.1, Center(-1.0, 3.0 * alpha)),
          Cylinder(r_tube, 2.1, Center(-1.0, 4.0 * alpha)),
          Cylinder(r_tube, 2.1, Center(-1.0, 5.0 * alpha))),
      amrex::EB2::CylinderIF(r_inner, 1.0, 0, {0.25, 0.0, 0.0}, false));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  PatchHierarchyOptions hierarchy_options(
      fub::GetOptions(options, "PatchHierarchy"));

  BOOST_LOG(log) << "Build EB level set.";

  hierarchy_options.index_spaces =
      MakeIndexSpaces(shop, grid_geometry, hierarchy_options);
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  using namespace std::literals;
  const fub::ProgramOptions initial_options =
      fub::GetOptions(options, "InitialCondition");
  std::string source{"initial_data.h5"};
  source = fub::GetOptionOr(initial_options, "data", source);

  InterpolateFrom1d initial_data(equation, std::move(source));

  using State = fub::Complete<fub::PerfectGas<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.01}};

  // ::amrex::RealBox outlet{{0.5, -0.5, -0.5}, {0.54, +0.5, +0.5}};
  // const ::amrex::Box outlet_box =
  //     fub::amrex::BoxWhichContains(outlet, coarse_geometry);

  PressureOutflowOptions boundary_options =
      fub::GetOptions(options, "PressureBoundary");
  BOOST_LOG(log) << "PressureBoundary:";
  boundary_options.Print(log);
  auto seq = fub::execution::seq;

  BoundarySet boundary_condition{
      {ReflectiveBoundary{seq, equation, fub::Direction::X, 0},
       PressureOutflowBoundary{equation, boundary_options},
       ReflectiveBoundary{seq, equation, fub::Direction::Z, 0},
       ReflectiveBoundary{seq, equation, fub::Direction::Z, 1},
       ReflectiveBoundary{seq, equation, fub::Direction::Y, 0},
       ReflectiveBoundary{seq, equation, fub::Direction::Y, 1}}};
  //  perfect_gas::IsentropicPressureBoundary{equation, boundary_options}}};

  // If a checkpoint path is specified we will fill the patch hierarchy with
  // data from the checkpoint file, otherwise we will initialize the data by
  // the initial data function.
  std::shared_ptr gridding = [&] {
    std::string checkpoint{};
    if (options.count("checkpoint")) {
      checkpoint = options.at("checkpoint").cast<std::string>();
    }
    if (checkpoint.empty()) {
      BOOST_LOG(log) << "Initialize hierarchy.";
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(equation, grid_geometry, hierarchy_options),
          initial_data, TagAllOf(TagCutCells(), gradients, TagBuffer(2)),
          boundary_condition);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      BOOST_LOG(log) << "Initialize from checkpoint.";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), grid_geometry,
          hierarchy_options);
      return std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), gradients, TagBuffer(2)), boundary_condition);
    }
  }();

  // Make Solver

  // fub::EinfeldtSignalVelocities<fub::PerfectGas<3>> signals;
  // fub::HllMethod hll_method{equation, signals};v
//  fub::perfect_gas::HllemMethod<3> hllem_method{equation};
  fub::FluxMethod<fub::perfect_gas::MusclHancockPrim<3>> muscl_method{equation};
  fub::KbnCutCellMethod cutcell_method(muscl_method);

  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  const int scratch_gcw = 4;
  const int flux_gcw = 2;

  return IntegratorContext(gridding, method, scratch_gcw, flux_gcw);
}

void WriteCheckpoint(const std::string& path,
                     const fub::amrex::cutcell::GriddingAlgorithm& grid) {
  std::ptrdiff_t cycles = grid.GetCycles();
  std::string name = fmt::format("{}/{:09}", path, cycles);
  fub::amrex::cutcell::WriteCheckpointFile(name, grid.GetPatchHierarchy());
}

struct CheckpointOutput
    : fub::OutputAtFrequencyOrInterval<fub::amrex::cutcell::GriddingAlgorithm> {
  CheckpointOutput(const fub::ProgramOptions& options)
      : OutputAtFrequencyOrInterval(options), log{fub::GetInfoLogger()} {
    directory_ = fub::GetOptionOr(options, "directory", directory_);
    BOOST_LOG(log) << "CheckpointOutput configured:";
    BOOST_LOG(log) << " - directory = " << directory_;
  }

  void operator()(const fub::amrex::cutcell::GriddingAlgorithm& grid) override {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", grid.GetTimePoint().count());
    BOOST_LOG(log) << fmt::format("Write checkpoint to '{}'.", directory_);
    WriteCheckpoint(directory_, grid);
  }

  std::string directory_{"./Checkpoint/"};
  fub::SeverityLogger log;
};

void MyMain(const fub::ProgramOptions& options) {
  using namespace fub::amrex::cutcell;
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};

  fub::PerfectGas<Plenum_Rank> equation{};

  IntegratorContext plenum = MakePlenumSolver(equation, options);

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<Plenum_Rank>, std::move(plenum), fub::StrangSplitting{});

  // fub::NoSubcycleSolver solver(std::move(level_integrator));
  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<WriteHdf5>("HDF5");
  factory.RegisterOutput<PlotfileOutput<fub::PerfectGas<3>>>("Plotfile",
                                                             equation);
  factory.RegisterOutput<CheckpointOutput>("Checkpoint");
  using CounterOutput =
      fub::CounterOutput<GriddingAlgorithm, std::chrono::milliseconds>;
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<DebugOutput>("DebugOutput");
  factory.RegisterOutput<PerfectGasProbesOutput>("ProbesOutput");
  fub::MultipleOutputs<GriddingAlgorithm> outputs(
      std::move(factory), fub::GetOptions(options, "Output"));

  fub::SeverityLogger log = fub::GetInfoLogger();
  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);
  outputs(*solver.GetGriddingAlgorithm());
  BOOST_LOG(log) << "Start simulation.";
  fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
  pybind11::scoped_interpreter interpreter{};
  std::optional<fub::ProgramOptions> opts = fub::ParseCommandLine(argc, argv);
  if (opts) {
    // fub::EnableFloatingPointExceptions();
    MyMain(*opts);
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}
