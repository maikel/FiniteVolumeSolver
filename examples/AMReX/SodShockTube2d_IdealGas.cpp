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

struct ShockTubeData {
  using Equation = fub::IdealGasMix<2>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  void InitializeData(amrex::MultiFab& data, const amrex::Geometry& geom) {
    fub::amrex::ForEachFab(data, [&](const amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());

      fub::FlameMasterReactor& reactor = equation_.GetReactor();
      reactor.SetDensity(1.0);
      reactor.SetMoleFractions("N2:79,O2:21");

      auto SetPressure = [](fub::FlameMasterReactor& reactor, double pressure) {
        reactor.SetTemperature(pressure / (reactor.GetDensity() * reactor.GetUniversalGasConstant() / reactor.GetMeanMolarMass()));
      };
      Complete state(equation_);
      ForEachIndex(fub::Box<0>(states),
                   [&](std::ptrdiff_t i, std::ptrdiff_t j) {
                     double xy[AMREX_SPACEDIM];
                     geom.CellCenter({AMREX_D_DECL(int(i), int(j), 0)}, xy);
                     const double x = xy[0];
                     const double y = xy[1];

                     // "Left" states of Sod Shock Tube.
                     if (x + y < 0.0) {
                       reactor.SetDensity(1.0);
                       SetPressure(reactor, 1.01325e5);
                     }
                     // "Right" states.
                     else {
                       reactor.SetDensity(1.0e-1);
                       SetPressure(reactor, 1.25e-1 * 1.01325e5);
                     }
                     equation_.CompleteFromReactor(state);
                     equation_.CompleteFromCons(state, state);
                     Store(states, state, {i, j});
                   });
    });
  }

  Equation equation_;
};

int main() {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::amrex::ScopeGuard guard{};
  fub::InitializeLogging(MPI_COMM_WORLD);

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM >= 2);

  fub::Burke2012 burke{};
  fub::IdealGasMix<2> equation{burke};

  fub::amrex::CartesianGridGeometry geometry{};
  geometry.cell_dimensions = std::array<int, Dim>{AMREX_D_DECL(128, 128, 1)};
  geometry.coordinates = amrex::RealBox({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                                        {AMREX_D_DECL(+1.0, +1.0, +1.0)});
//  geometry.periodicity = std::array<int, Dim>{1, 1};

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 1;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 2, 1)};
  hier_opts.blocking_factor = amrex::IntVect{AMREX_D_DECL(8, 8, 1)};

  using Complete = fub::IdealGasMix<2>::Complete;
  fub::amrex::GradientDetector gradient(equation,
                                        std::pair{&Complete::density, 1.0e-2},
                                        std::pair{&Complete::pressure, 1.0e-2});

  using fub::amrex::ReflectiveBoundary;
  fub::amrex::BoundarySet boundaries{
      {ReflectiveBoundary(equation, fub::Direction::X, 0),
       ReflectiveBoundary(equation, fub::Direction::X, 1),
       ReflectiveBoundary(equation, fub::Direction::Y, 0),
       ReflectiveBoundary(equation, fub::Direction::Y, 1)}};

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts),
      ShockTubeData{equation},
      fub::amrex::TagAllOf(gradient, fub::amrex::TagBuffer(4)), boundaries);
  gridding->InitializeHierarchy(0.0);

//  fub::EinfeldtSignalVelocities<fub::IdealGasMix<2>> signals{};
//  fub::HllMethod hll_method(equation, signals);
  fub::ideal_gas::MusclHancockPrimMethod<2> muscl_method{equation};
  // fub::amrex::HyperbolicMethod method{fub::amrex::FluxMethod(hll_method),
  fub::amrex::HyperbolicMethod method{fub::amrex::FluxMethod(muscl_method),
                                      fub::amrex::EulerForwardTimeIntegrator(),
                                      fub::amrex::Reconstruction(equation)};

//  const int base_gcw = hll_method.GetStencilWidth();
  const int base_gcw = muscl_method.GetStencilWidth();
//  const int scratch_ghost_cell_width = 1 * base_gcw;
//  const int flux_ghost_cell_width = 0 * base_gcw;
  const int scratch_ghost_cell_width = 2 * base_gcw;
  const int flux_ghost_cell_width = 1 * base_gcw;

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>,
      fub::amrex::IntegratorContext(gridding, method, scratch_ghost_cell_width,
                                    flux_ghost_cell_width),
    // fub::GodunovSplitting());
       fub::StrangSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  std::string base_name = "SodShockTube_IdealGas_muscl_prim/";

  using namespace fub::amrex;
  using namespace std::literals::chrono_literals;

  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);

  // Define a function to compute the mass in the domain and output the
  // difference to the intial mass.
  double mass0 = 0.0;
  auto compute_mass = [&log, &mass0](const GriddingAlgorithm& grid) {
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
    BOOST_LOG(log) << boost::log::add_value("Time", time_point)
                   << fmt::format("Conservation Error in Mass: {:.6e}",
                                  mass_error);
  };

  fub::MultipleOutputs<GriddingAlgorithm> output{};

  // Add output to show the conservation error in mass after each time step
  output.AddOutput(fub::MakeOutput<GriddingAlgorithm>({1}, {}, compute_mass));

  // Add output to write AMReX plotfiles in a set time interval
  output.AddOutput(std::make_unique<fub::amrex::PlotfileOutput<fub::IdealGasMix<2>>>(std::vector<std::ptrdiff_t>{},
          std::vector<fub::Duration>{0.0001s}, equation, base_name));

  // Add output for the timer database after each 25 cycles
  output.AddOutput(
      std::make_unique<
          fub::CounterOutput<GriddingAlgorithm, std::chrono::milliseconds>>(
          wall_time_reference, std::vector<std::ptrdiff_t>{25},
          std::vector<fub::Duration>{}));

  using namespace std::literals::chrono_literals;
  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.final_time = 0.01s;
  run_options.cfl = 0.4;
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}
