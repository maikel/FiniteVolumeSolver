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

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration) {
    amrex::MultiFab& data = patch_level.data;
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    InitializeData(data, geom);
  }

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

using FactoryFunction =
    std::function<fub::AnyFluxMethod<fub::amrex::IntegratorContext>(
        const fub::IdealGasMix<1>&)>;

template <typename... Pairs> auto GetFluxMethodFactory(Pairs... ps) {
  std::map<std::string, FactoryFunction> factory;
  ((factory[ps.first] = ps.second), ...);
  return factory;
}

template <typename FluxMethod> struct MakeFlux {
  fub::AnyFluxMethod<fub::amrex::IntegratorContext>
  operator()(const fub::IdealGasMix<1>& eq) const {
    FluxMethod flux_method{eq};
    fub::amrex::FluxMethodAdapter adapter(std::move(flux_method));
    return adapter;
  }
};

// double MinMod(const double a, const double b) {
// 	return (a * b > 0.0) ? ((fabs(a) < fabs(b)) ? a : b) : 0.0;
// }

// template <>
// struct MakeFlux<fub::ideal_gas::MusclHancockCharMethod<1>> {
//   fub::AnyFluxMethod<fub::amrex::IntegratorContext>
//   operator()(const fub::IdealGasMix<1>& eq) const {
//     fub::ideal_gas::MusclHancockCharMethod<1> flux_method{eq, &MinMod};
//     fub::amrex::FluxMethodAdapter adapter(fub::execution::seq, std::move(flux_method));
//     return adapter;
//   }
// };

void MyMain(const fub::ProgramOptions& opts) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard guard{};
  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::Burke2012 burke{};
  fub::IdealGasMix<1> equation{burke};

  fub::amrex::CartesianGridGeometry geometry =
      fub::GetOptions(opts, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  geometry.Print(log);

  fub::amrex::PatchHierarchyOptions hier_opts =
      fub::GetOptions(opts, "PatchHierarchy");
  BOOST_LOG(log) << "PatchHierarchy:";
  hier_opts.Print(log);

  using Complete = fub::IdealGasMix<1>::Complete;
  fub::amrex::GradientDetector gradient{
      equation, std::make_pair(&Complete::density, 5e-3),
      std::make_pair(&Complete::pressure, 5e-2)};

  fub::FlameMasterReactor& reactor = equation.GetReactor();
  const double M = reactor.GetMeanMolarMass();
  const double R = reactor.GetUniversalGasConstant();
  const double Rspec = R / M;
  const double atm = 101325.0;
  const double p_left = 1.0 * atm;
  const double rho_left = 1.0;
  const double T_left = p_left / (rho_left * Rspec);
  reactor.SetDensity(rho_left);
  reactor.SetTemperature(T_left);
  Complete left(equation);
  equation.CompleteFromReactor(left);

  const double p_right = 0.1 * atm;
  const double rho_right = 0.125;
  const double T_right = p_right / (rho_right * Rspec);
  reactor.SetDensity(rho_right);
  reactor.SetTemperature(T_right);
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

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts), initial_data,
      gradient, boundary);
  gridding->InitializeHierarchy(0.0);

  using namespace std::literals;
  using ConservativeReconstruction = fub::MusclHancockMethod<
      fub::IdealGasMix<1>,
      fub::HllMethod<fub::IdealGasMix<1>, fub::EinfeldtSignalVelocities<fub::IdealGasMix<1>>>, fub::VanLeer>;
  using PrimitiveReconstruction = fub::ideal_gas::MusclHancockPrimMethod<1>;
  using CharacteristicReconstruction = fub::ideal_gas::MusclHancockCharMethod<1>;

  auto flux_method_factory = GetFluxMethodFactory(
      std::pair{"Conservative"s, MakeFlux<ConservativeReconstruction>()},
      std::pair{"Primitive"s, MakeFlux<PrimitiveReconstruction>()},
      std::pair{"Characteristics"s, MakeFlux<CharacteristicReconstruction>()});

  std::string reconstruction =
      fub::GetOptionOr(opts, "reconstruction", "Primitive"s);
  auto flux_method = flux_method_factory.at(reconstruction)(equation);

  fub::amrex::HyperbolicMethod method{std::move(flux_method),
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

  using namespace std::literals::chrono_literals;
  using Plotfile = fub::amrex::PlotfileOutput<fub::IdealGasMix<1>>;
  using CounterOutput = fub::CounterOutput<fub::amrex::GriddingAlgorithm,
                                           std::chrono::milliseconds>;
  fub::OutputFactory<fub::amrex::GriddingAlgorithm> factory{};
  factory.RegisterOutput<Plotfile>("Plotfile", equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<fub::amrex::WriteHdf5>("HDF5");
  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> output(
      std::move(factory), fub::GetOptions(opts, "Output"));

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