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
// OUT OF OR ILN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "fub/AMReX.hpp"
#include "fub/AMReX_CutCell.hpp"
#include "fub/Solver.hpp"

#include "fub/equations/PerfectGasMix.hpp"
#include "fub/flux_method/MusclHancockMethod2.hpp"

#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

auto Rectangle(const std::array<double, 2>& lower,
               const std::array<double, 2>& upper) {
  amrex::EB2::PlaneIF lower_x({lower[0], lower[1]}, {0, +1});
  amrex::EB2::PlaneIF lower_y({lower[0], lower[1]}, {+1, 0});
  amrex::EB2::PlaneIF upper_x({upper[0], upper[1]}, {0, -1});
  amrex::EB2::PlaneIF upper_y({upper[0], upper[1]}, {-1, 0});
  return amrex::EB2::makeIntersection(lower_x, lower_y, upper_x, upper_y);
}

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};
  fub::SeverityLogger log = fub::GetInfoLogger();

  auto embedded_boundary =
      amrex::EB2::makeUnion(Rectangle({-1.0, +0.015}, {0.0, 1.0}),
                            Rectangle({-1.0, -1.0}, {0.0, -0.015}));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::PerfectGasMix<2> equation(mech);

  fub::amrex::CartesianGridGeometry grid_geometry(
      fub::GetOptions(options, "GridGeometry"));

  using namespace fub::amrex::cutcell;
  PatchHierarchyOptions hierarchy_options(
      fub::GetOptions(options, "PatchHierarchy"));

  BOOST_LOG(log) << "Compute EB level set data...";
  std::shared_ptr<fub::CounterRegistry> registry =
      std::make_shared<fub::CounterRegistry>();
  {
    fub::Timer timer = registry->get_timer("MakeIndexSpaces");
    hierarchy_options.index_spaces =
        MakeIndexSpaces(shop, grid_geometry, hierarchy_options);
  }

  fub::KineticState<fub::PerfectGasMix<2>> kin;
  kin.temperature = 1.0;
  kin.density = 1.0;
  kin.mole_fractions[0] = 1.0;
  kin.mole_fractions[1] = 0.0;
  kin.mole_fractions[2] = 0.0;

  fub::Complete<fub::PerfectGasMix<2>> right{equation};
  fub::euler::CompleteFromKineticState(equation, right, kin);

  kin.density = 1.0 / 4.0;
  fub::Complete<fub::PerfectGasMix<2>> left{equation};  
  fub::euler::CompleteFromKineticState(equation, left, kin);

  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04),
                              left, right);

  using State = fub::Complete<fub::PerfectGasMix<2>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.05}};

  auto seq = fub::execution::seq;
  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  ReflectiveBoundary{seq, equation, fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, grid_geometry, hierarchy_options), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  using HLLEM = fub::perfect_gas::HllemMethod<fub::PerfectGasMix<2>>;

  using CharacteristicsReconstruction = fub::FluxMethod<fub::MusclHancock2<
      fub::PerfectGasMix<2>,
      fub::CharacteristicsGradient<
          fub::PerfectGasMix<2>,
          fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
      fub::CharacteristicsReconstruction<fub::PerfectGasMix<1>>, HLLEM>>;

  HLLEM hllem_method{equation};
  CharacteristicsReconstruction flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hllem_method);

  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method, 4, 2),
      fub::GodunovSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  using namespace std::literals::chrono_literals;
  using Plotfile = PlotfileOutput<fub::PerfectGasMix<2>>;
  using CounterOutput =
      fub::CounterOutput<GriddingAlgorithm, std::chrono::milliseconds>;
  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<Plotfile>("Plotfile", equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<WriteHdf5>("HDF5");
  fub::MultipleOutputs<GriddingAlgorithm> output(
      std::move(factory), fub::GetOptions(options, "Output"));
  output(*gridding);
  fub::RunOptions run_options(fub::GetOptions(options, "RunOptions"));
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
