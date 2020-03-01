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

#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Translation.H>

#include <iostream>

void MyMain(const fub::ProgramOptions& options) {
  using namespace fub::amrex::cutcell;
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard _{};

  auto MakePolygon = [](auto... points) {
    fub::Polygon::Vector xs{};
    fub::Polygon::Vector ys{};

    (xs.push_back(std::get<0>(points)), ...);
    (ys.push_back(std::get<1>(points)), ...);

    return fub::Polygon(std::move(xs), std::move(ys));
  };

  auto DivergentInlet = [&](double height, const std::array<double, 3>& center) {
    const double xlo = center[0] - height;
    const double xhi = center[0];
    const double xdiv = xhi - 0.075;
    const double r = 0.015;
    const double r2 = 0.015;
    auto polygon = MakePolygon(std::pair{xlo, r}, std::pair{xdiv, r}, std::pair{xhi, r2}, 
                               std::pair{xhi, -r2}, std::pair{xdiv, -r}, std::pair{xlo, -r}, std::pair{xlo, r});
    auto tube_in_zero = fub::amrex::Geometry(fub::Invert(fub::RotateAxis(polygon)));
    amrex::RealArray real_center{center[0], center[1], center[2]};
    auto tube_in_center = amrex::EB2::TranslationIF(fub::amrex::Geometry(tube_in_zero), real_center);
    return tube_in_center;
  };

  auto embedded_boundary = amrex::EB2::makeIntersection(
      DivergentInlet(1.5, {0.0, 0.0, 0.0}),
      amrex::EB2::SphereIF(0.150, {0.0, -0.15 + 0.01501, 0.0}, true));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::amrex::CartesianGridGeometry grid_geometry =
      fub::GetOptions(options, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  PatchHierarchyOptions hierarchy_options =
      fub::GetOptions(options, "PatchHierarchy");
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);
    
  BOOST_LOG(log) << "Build EB Level set...";
  hierarchy_options.index_spaces =
      MakeIndexSpaces(shop, grid_geometry, hierarchy_options);

//   fub::Burke2012 mechanism;
//   fub::IdealGasMix<3> equation{mechanism};
   fub::PerfectGas<3> equation{};

//   fub::Complete<fub::IdealGasMix<3>> left(equation);
  fub::Conservative<fub::PerfectGas<3>> cons;
  cons.density = 1.0;
  cons.momentum << 0.0, 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<3>> right;
  fub::CompleteFromCons(equation, right, cons);

  cons.energy *= 4;
  fub::Complete<fub::PerfectGas<3>> left;
  fub::CompleteFromCons(equation, left, cons);

//   fub::Complete<fub::IdealGasMix<3>> right(equation);
//   {
//     using namespace std::literals;
//     const fub::ProgramOptions initial_options =
//         fub::GetOptions(options, "InitialCondition");
//     const fub::ProgramOptions left_options =
//         fub::GetOptions(initial_options, "left");
//     std::string moles = fub::GetOptionOr(left_options, "moles", "N2:79,O2:21"s);
//     double temperature = fub::GetOptionOr(left_options, "temperature", 300.0);
//     double pressure = fub::GetOptionOr(left_options, "pressure", 101325.0);
//     std::array<double, 3> velocity{};
//     velocity = fub::GetOptionOr(left_options, "velocity", velocity);
//     equation.GetReactor().SetMoleFractions(moles);
//     equation.GetReactor().SetTemperature(temperature);
//     equation.GetReactor().SetPressure(pressure);
//     equation.CompleteFromReactor(left);
//     fub::CompleteFromCons(equation, left, left);

//     const fub::ProgramOptions right_options =
//         fub::GetOptions(initial_options, "right");
//     moles = fub::GetOptionOr(right_options, "moles", "N2:79,O2:21"s);
//     temperature = fub::GetOptionOr(right_options, "temperature", 300.0);
//     pressure = fub::GetOptionOr(right_options, "pressure", 101325.0);
//     velocity = std::array<double, 3>{};
//     velocity = fub::GetOptionOr(right_options, "velocity", velocity);
//     equation.GetReactor().SetMoleFractions(moles);
//     equation.GetReactor().SetTemperature(temperature);
//     equation.GetReactor().SetPressure(pressure);
//     equation.CompleteFromReactor(right);
//     fub::CompleteFromCons(equation, right, right);
//   }
  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.6),
                              left, right);

  using Complete = fub::Complete<fub::PerfectGas<3>>;
//   using Complete = fub::Complete<fub::IdealGasMix<3>>;
  GradientDetector gradients{equation, std::pair{&Complete::pressure, 0.05},
                             std::pair{&Complete::density, 0.05}};

  auto seq = fub::execution::seq;
  BoundarySet boundary_condition{{ReflectiveBoundary{seq, equation, fub::Direction::X, 0},
                                  ReflectiveBoundary{seq, equation, fub::Direction::Z, 0},
                                  ReflectiveBoundary{seq, equation, fub::Direction::Z, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, grid_geometry, hierarchy_options), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(2)), boundary_condition);

  BOOST_LOG(log) << "Initialize Hierarchy...";
  gridding->InitializeHierarchy(0.0);

//   fub::EinfeldtSignalVelocities<fub::IdealGasMix<3>> signals{};
  fub::EinfeldtSignalVelocities<fub::PerfectGas<3>> signals{};
  fub::HllMethod hll_method{equation, signals};
  // fub::MusclHancockMethod flux_method(equation, hll_method);
  // fub::ideal_gas::MusclHancockPrimMethod<3> flux_method(equation);
  fub::FluxMethod<fub::perfect_gas::MusclHancockPrim<3>> flux_method{equation};
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  const int scratch_gcw = 4;
  const int flux_gcw = 2;

  IntegratorContext context(gridding, method, scratch_gcw, flux_gcw);

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<3>, std::move(context), fub::StrangSplitting{});

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  using namespace std::literals::chrono_literals;
//   using Plotfile = PlotfileOutput<fub::IdealGasMix<3>>;
  using Plotfile = PlotfileOutput<fub::PerfectGas<3>>;
  using CounterOutput =
      fub::CounterOutput<GriddingAlgorithm, std::chrono::milliseconds>;
  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<Plotfile>("Plotfile", equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);

  fub::MultipleOutputs<GriddingAlgorithm> outputs(
      std::move(factory), fub::GetOptions(options, "Output"));

  outputs(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);
  fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
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
