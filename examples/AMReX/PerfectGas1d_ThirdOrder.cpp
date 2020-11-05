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

#include <cmath>
#include <fmt/format.h>
#include <iostream>

#include <fenv.h>

#include <boost/log/utility/manipulators/add_value.hpp>

struct SinusProblem {
  using Equation = fub::PerfectGas<1>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  Equation equation_;

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration /*time*/) const {
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& data = patch_level.data;
    fub::amrex::ForEachFab(
        fub::execution::openmp, data, [&](amrex::MFIter& mfi) {
          fub::View<Complete> state = fub::amrex::MakeView<Complete>(
              data[mfi], equation_, mfi.tilebox());
          fub::ForEachIndex(fub::Box<0>(state), [this, &state,
                                                 &geom](std::ptrdiff_t i) {
            const double x = geom.CellCenter(int(i), 0);
            double temperature = (0.95 < x && x < 1.05) ? 310.0 : 300.0;
            const double pressure = 101325.0;
            const double density = pressure / (equation_.Rspec * temperature);
            fub::Array<double, 1, 1> velocity{1.0};
            Complete complete =
                equation_.CompleteFromPrim(density, velocity, pressure);
            fub::Store(state, complete, {i});
          });
        });
  }
};

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard scope_guard{};

  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::ProgramOptions equation_options = fub::GetOptions(options, "Equation");
  
  fub::PerfectGas<1> equation{};
  equation.Rspec = fub::GetOptionOr(equation_options, "R_specific", equation.Rspec);
  equation.gamma = fub::GetOptionOr(equation_options, "gamma", equation.gamma);
  equation.gamma_minus_1_inv = 1.0 / (equation.gamma - 1.0);
  equation.gamma_array_ = fub::Array1d::Constant(equation.gamma);
  equation.gamma_minus_1_inv_array_ = fub::Array1d::Constant(equation.gamma_minus_1_inv);

  fub::perfect_gas::ThirdOrderRungeKuttaMethod<1> runge_kutta{equation};
  runge_kutta.k_0 = 0.0; // 0.017746 / std::sqrt(300.0) * 5000.0;
  runge_kutta.eta_0 = 0.0; // 2.23e-5 / std::sqrt(293.0) * 5000.0;

  runge_kutta.k_0 = fub::GetOptionOr(equation_options, "k_0", runge_kutta.k_0);
  runge_kutta.eta_0 = fub::GetOptionOr(equation_options, "eta_0", runge_kutta.eta_0);

  BOOST_LOG(log) << "Equation:";
  BOOST_LOG(log) << fmt::format(" - R_specific = {}", equation.Rspec);
  BOOST_LOG(log) << fmt::format(" - gamma = {}", equation.gamma);
  BOOST_LOG(log) << fmt::format(" - k_0 = {}", runge_kutta.k_0);
  BOOST_LOG(log) << fmt::format(" - eta_0 = {}", runge_kutta.eta_0);

  fub::amrex::CartesianGridGeometry grid_geometry(
      fub::GetOptions(options, "GridGeometry"));
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  using Complete = fub::PerfectGas<1>::Complete;
  fub::amrex::GradientDetector gradient{
      equation, std::make_pair(&Complete::density, 5e-3),
      std::make_pair(&Complete::pressure, 5e-2)};

  SinusProblem initial_data{equation};

  fub::amrex::PatchHierarchyOptions hierarchy_options(
      fub::GetOptions(options, "PatchHierarchy"));
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, grid_geometry, hierarchy_options), initial_data,
      gradient);
  gridding->InitializeHierarchy(0.0);

  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethodAdapter(fub::execution::seq, runge_kutta),
      fub::amrex::EulerForwardTimeIntegrator(),
      fub::amrex::Reconstruction(equation)};

  const int scratch_ghost_cell_width = 6;
  const int flux_ghost_cell_width = 0;

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<1>,
      fub::amrex::IntegratorContext(gridding, method, scratch_ghost_cell_width,
                                    flux_ghost_cell_width),
      fub::GodunovSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  using namespace fub::amrex;
  using namespace std::literals::chrono_literals;

  fub::OutputFactory<fub::amrex::GriddingAlgorithm> factory{};
  factory.RegisterOutput<fub::amrex::WriteHdf5>("HDF5");
  using CounterOutput =
      fub::CounterOutput<fub::amrex::GriddingAlgorithm,
                         std::chrono::milliseconds>;
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> outputs(
      std::move(factory), fub::GetOptions(options, "Output"));

  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);

  outputs(*solver.GetGriddingAlgorithm());
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