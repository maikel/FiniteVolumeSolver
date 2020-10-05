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

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration) {
    amrex::MultiFab& data = patch_level.data;
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    InitializeData(data, geom);
  }

  void InitializeData(amrex::MultiFab& data, const amrex::Geometry& geom) {
    fub::amrex::ForEachFab(data, [&](const amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());

      fub::FlameMasterReactor& reactor = equation_.GetReactor();
      reactor.SetDensity(1.0);
      reactor.SetMoleFractions("N2:79,O2:21");

      auto SetPressure = [](fub::FlameMasterReactor& reactor, double pressure) {
        reactor.SetTemperature(pressure / (reactor.GetDensity() *
                                           reactor.GetUniversalGasConstant() /
                                           reactor.GetMeanMolarMass()));
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
                       reactor.SetDensity(1.25e-1);
                       SetPressure(reactor, 1.e-1 * 1.01325e5);
                     }
                     equation_.CompleteFromReactor(state);
                     equation_.CompleteFromCons(state, state);
                     Store(states, state, {i, j});
                   });
    });
  }

  Equation equation_;
};

using FactoryFunction =
    std::function<fub::AnyFluxMethod<fub::amrex::IntegratorContext>(
        const fub::IdealGasMix<2>&)>;

template <typename... Pairs> auto GetFluxMethodFactory(Pairs... ps) {
  std::map<std::string, FactoryFunction> factory;
  ((factory[ps.first] = ps.second), ...);
  return factory;
}

template <typename FluxMethod> struct MakeFlux {
  fub::AnyFluxMethod<fub::amrex::IntegratorContext>
  operator()(const fub::IdealGasMix<2>& eq) const {
    FluxMethod flux_method{eq};
    fub::amrex::FluxMethodAdapter adapter(std::move(flux_method));
    return adapter;
  }
};


void MyMain(const fub::ProgramOptions& opts) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard guard{};
  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::Burke2012 burke{};
  fub::IdealGasMix<2> equation{burke};

  fub::amrex::CartesianGridGeometry geometry =
      fub::GetOptions(opts, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  geometry.Print(log);

  fub::amrex::PatchHierarchyOptions hier_opts =
      fub::GetOptions(opts, "PatchHierarchy");
  BOOST_LOG(log) << "PatchHierarchy:";
  hier_opts.Print(log);

  using Complete = fub::IdealGasMix<2>::Complete;
  fub::amrex::GradientDetector gradient(equation,
                                        std::pair{&Complete::density, 1.0e-2});

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

  using namespace std::literals;
  using HLLE = fub::HllMethod<fub::IdealGasMix<2>, fub::EinfeldtSignalVelocities<fub::IdealGasMix<2>>>;
  using ConservativeReconstruction = fub::MusclHancockMethod<fub::IdealGasMix<2>, HLLE, fub::VanLeer>;
  using PrimitiveReconstruction = fub::ideal_gas::MusclHancockPrimMethod<2>;
  using CharacteristicReconstruction = fub::ideal_gas::MusclHancockCharMethod<2>;

  auto flux_method_factory = GetFluxMethodFactory(
      std::pair{"HLLE"s, MakeFlux<HLLE>()},
      std::pair{"Conservative"s, MakeFlux<ConservativeReconstruction>()},
      std::pair{"Primitive"s, MakeFlux<PrimitiveReconstruction>()},
      std::pair{"Characteristics"s, MakeFlux<CharacteristicReconstruction>()});

  std::string reconstruction =
      fub::GetOptionOr(opts, "reconstruction", "Primitive"s);
  auto flux_method = flux_method_factory.at(reconstruction)(equation);

  fub::amrex::HyperbolicMethod method{flux_method,
  // fub::ideal_gas::MusclHancockPrimMethod<2> muscl_method{equation};
  //fub::amrex::HyperbolicMethod method{
    //  fub::amrex::FluxMethodAdapter(muscl_method),
      fub::amrex::EulerForwardTimeIntegrator(),
      fub::amrex::Reconstruction(equation)};

  const int base_gcw = flux_method.GetStencilWidth();
  const int scratch_ghost_cell_width = 4 * base_gcw;
  const int flux_ghost_cell_width = 3 * base_gcw;

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>,
      fub::amrex::IntegratorContext(gridding, method, scratch_ghost_cell_width,
                                    flux_ghost_cell_width),
      // fub::GodunovSplitting());
      fub::StrangSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  using namespace fub::amrex;
  using namespace std::literals::chrono_literals;

  // Define a function to compute the mass in the domain and output the
  // difference to the intial mass.
  auto compute_mass = [&log, mass0 = 0.0](const GriddingAlgorithm& grid) mutable {
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

  using namespace std::literals::chrono_literals;
  using Plotfile = fub::amrex::PlotfileOutput<fub::IdealGasMix<2>>;
  using CounterOutput = fub::CounterOutput<fub::amrex::GriddingAlgorithm,
                                           std::chrono::milliseconds>;
  fub::OutputFactory<fub::amrex::GriddingAlgorithm> factory{};
  factory.RegisterOutput<Plotfile>("Plotfile", equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<fub::amrex::WriteHdf5>("HDF5");
  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> output(
      std::move(factory), fub::GetOptions(opts, "Output"));

  // Add output to show the conservation error in mass after each time step
  output.AddOutput(fub::MakeOutput<GriddingAlgorithm>({1}, {}, compute_mass));

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options = fub::GetOptions(opts, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);
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
