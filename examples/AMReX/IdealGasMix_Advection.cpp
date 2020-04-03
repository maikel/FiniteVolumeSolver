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

#include <fmt/format.h>
#include <iostream>

#include <boost/log/utility/manipulators/add_value.hpp>

struct RiemannProblem {
  using Equation = fub::IdealGasMix<1>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  Equation equation_;
  Complete left_{equation_};
  Complete right_{equation_};

  void InitializeData(amrex::MultiFab& data, const amrex::Geometry& geom) {
    fub::amrex::ForEachFab(
        fub::execution::openmp, data, [&](amrex::MFIter& mfi) {
          fub::View<Complete> state = fub::amrex::MakeView<Complete>(
              data[mfi], equation_, mfi.tilebox());
          fub::ForEachIndex(fub::Box<0>(state),
                            [this, &state, &geom](std::ptrdiff_t i) {
                              const double x = geom.CellCenter(int(i), 0);
                              if (x < 0.0) {
                                Store(state, left_, {i});
                              } else {
                                Store(state, right_, {i});
                              }
                            });
        });
  }
};

int main() {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  // fub::EnableFloatingPointExceptions();
  const fub::amrex::ScopeGuard guard{};
  fub::InitializeLogging(MPI_COMM_WORLD);

  constexpr int Dim = AMREX_SPACEDIM;

  const std::array<int, Dim> n_cells{AMREX_D_DECL(200, 1, 1)};
  const std::array<double, Dim> xlower{AMREX_D_DECL(-1.0, -1.0, -1.0)};
  const std::array<double, Dim> xupper{AMREX_D_DECL(+1.0, +1.0, +1.0)};

  fub::Burke2012 burke{};
  fub::IdealGasMix<1> equation{burke};

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  using Complete = fub::IdealGasMix<1>::Complete;
  fub::amrex::GradientDetector gradient{
      equation, std::make_pair(&Complete::density, 5e-3),
      std::make_pair(&Complete::pressure, 5e-2)};

  fub::FlameMasterReactor& reactor = equation.GetReactor();
  reactor.SetTemperature(300.);
  reactor.SetPressure(101325.0);
  Complete left(equation);
  equation.CompleteFromReactor(left);

  reactor.SetTemperature(600.);
  reactor.SetPressure(0.1 * 101325.0);
  Complete right(equation);
  equation.CompleteFromReactor(right);

  RiemannProblem initial_data{equation, left, right};

  fub::amrex::BoundarySet boundary;
  using fub::amrex::ReflectiveBoundary;
  auto seq = fub::execution::seq;
  boundary.conditions.push_back(
      ReflectiveBoundary{seq, equation, fub::Direction::X, 0});
  boundary.conditions.push_back(
      ReflectiveBoundary{seq, equation, fub::Direction::X, 1});

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 1;
  hier_opts.refine_ratio = ::amrex::IntVect{AMREX_D_DECL(2, 1, 1)};
  hier_opts.blocking_factor = ::amrex::IntVect{AMREX_D_DECL(8, 1, 1)};

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts), initial_data,
      gradient, boundary);
  gridding->InitializeHierarchy(0.0);


  fub::ideal_gas::MusclHancockPrimMethod<1> muscl_prim{equation};
  fub::amrex::HyperbolicMethod method{fub::amrex::FluxMethod(muscl_prim),
                                      fub::amrex::EulerForwardTimeIntegrator(),
                                      fub::amrex::Reconstruction(equation)};

  const int scratch_ghost_cell_width = 2;
  const int flux_ghost_cell_width = 0;

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<1>,
      fub::amrex::IntegratorContext(gridding, method, scratch_ghost_cell_width,
                                    flux_ghost_cell_width),
      fub::GodunovSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  std::string base_name = "IdealGasMix1d_muscl_new/";

  using namespace fub::amrex;
  using namespace std::literals::chrono_literals;
  namespace log = boost::log;

  log::sources::severity_logger<log::trivial::severity_level> lg(
      log::keywords::severity = log::trivial::info);

  double mass0 = 0.0;
  auto conservation_error = [&lg, &mass0](const GriddingAlgorithm& grid) {
    const ::amrex::MultiFab& data =
        grid.GetPatchHierarchy().GetPatchLevel(0).data;
    const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(0);
    const double volume_per_cell = geom.CellSize(0) * geom.CellSize(1);
    const double density_sum = data.sum(0);
    const double mass = density_sum * volume_per_cell;
    if (mass0 == 0.0) {
      mass0 = mass;
    }
    const double mass_error = mass - mass0;
    const double time_point = grid.GetTimePoint().count();
    BOOST_LOG(lg) << log::add_value("Time", time_point)
                  << fmt::format("Conservation Error in Mass: {:.6e}",
                                 mass_error);
  };

  fub::MultipleOutputs<GriddingAlgorithm> output{};

  output.AddOutput(
      fub::MakeOutput<GriddingAlgorithm>({1}, {}, conservation_error));

  output.AddOutput(std::make_unique<PlotfileOutput<fub::IdealGasMix<1>>>(std::vector<std::ptrdiff_t>{},
          std::vector<fub::Duration>{0.00001s}, equation, base_name));

  output.AddOutput(
      std::make_unique<fub::CounterOutput<GriddingAlgorithm>>(
          wall_time_reference, std::vector<std::ptrdiff_t>{},
          std::vector<fub::Duration>{0.5s}));

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.cfl = 0.5;
  run_options.final_time = 0.001s;
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}
