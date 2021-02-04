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

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>

static_assert(AMREX_SPACEDIM == 2);

using Coord = Eigen::Vector2d;

Coord OrthogonalTo(const Coord& x) { return Coord{x[1], -x[0]}; }

auto Triangle(const Coord& p1, const Coord& p2, const Coord& p3) {
  Coord norm1 = OrthogonalTo(p2 - p1).normalized();
  Coord norm2 = OrthogonalTo(p3 - p2).normalized();
  Coord norm3 = OrthogonalTo(p1 - p3).normalized();
  amrex::EB2::PlaneIF plane1({p1[0], p1[1]}, {norm1[0], norm1[1]});
  amrex::EB2::PlaneIF plane2({p2[0], p2[1]}, {norm2[0], norm2[1]});
  amrex::EB2::PlaneIF plane3({p3[0], p3[1]}, {norm3[0], norm3[1]});
  return amrex::EB2::makeIntersection(plane1, plane2, plane3);
}

using FactoryFunction =
    std::function<fub::AnyFluxMethod<fub::amrex::cutcell::IntegratorContext>(
        const fub::PerfectGas<2>&)>;

template <typename... Pairs> auto GetFluxMethodFactory(Pairs... ps) {
  std::map<std::string, FactoryFunction> factory;
  ((factory[ps.first] = ps.second), ...);
  return factory;
}

template <typename FluxMethod> struct MakeFlux {
  fub::AnyFluxMethod<fub::amrex::cutcell::IntegratorContext>
  operator()(const fub::PerfectGas<2>& eq) const {
    fub::EinfeldtSignalVelocities<fub::PerfectGas<2>> signals{};
    fub::HllMethod hll_method{eq, signals};
    FluxMethod flux_method{eq};
    fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);
    fub::amrex::cutcell::FluxMethod adapter(std::move(cutcell_method));
    return adapter;
  }
};

void MyMain(const fub::ProgramOptions& opts) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard guard{};
  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::PerfectGas<2> equation{};

  fub::amrex::CartesianGridGeometry geometry =
      fub::GetOptions(opts, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  geometry.Print(log);

  fub::amrex::cutcell::PatchHierarchyOptions hier_opts =
      fub::GetOptions(opts, "PatchHierarchy");
  BOOST_LOG(log) << "PatchHierarchy:";
  hier_opts.Print(log);

  BOOST_LOG(log) << "Compute EB level set data...";
  auto embedded_boundary =
      Triangle({0.02, 0.05}, {0.05, 0.0655}, {0.05, 0.0345});
  auto shop = amrex::EB2::makeShop(embedded_boundary);
  hier_opts.index_spaces =
      MakeIndexSpaces(shop, geometry, hier_opts);

  using Complete = fub::PerfectGas<2>::Complete;
  fub::amrex::cutcell::GradientDetector gradient(equation,
                                        std::pair{&Complete::density, 1.0e-2});

  using namespace std::literals;
  using HLLE = fub::HllMethod<fub::PerfectGas<2>, fub::EinfeldtSignalVelocities<fub::PerfectGas<2>>>;
  using HLLEM = fub::perfect_gas::HllemMethod<fub::PerfectGas<2>>;
  using ConservativeReconstruction = fub::MusclHancockMethod<fub::PerfectGas<2>, HLLE, fub::VanLeer>;
  using ConservativeReconstructionM = fub::MusclHancockMethod<fub::PerfectGas<2>, HLLEM, fub::VanLeer>;
  using PrimitiveReconstruction = fub::FluxMethod<fub::perfect_gas::MusclHancockPrim<2>>;
  using CharacteristicReconstruction = fub::perfect_gas::MusclHancockCharMethod<2>;

  auto flux_method_factory = GetFluxMethodFactory(
      std::pair{"HLLE"s, MakeFlux<HLLE>()},
      std::pair{"HLLEM"s, MakeFlux<HLLEM>()},
      std::pair{"Primitive"s, MakeFlux<PrimitiveReconstruction>()},
      std::pair{"Conservative"s, MakeFlux<ConservativeReconstruction>()},
      std::pair{"ConservativeM"s, MakeFlux<ConservativeReconstructionM>()},
      std::pair{"Characteristics"s, MakeFlux<CharacteristicReconstruction>()});

  std::string reconstruction =
      fub::GetOptionOr(opts, "reconstruction", "Characteristics"s);
  BOOST_LOG(log) << "Reconstruction: " << reconstruction;
  auto flux_method = flux_method_factory.at(reconstruction)(equation);

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.0;
  cons.momentum << 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> right;
  fub::CompleteFromCons(equation, right, cons);

  cons.energy *= 5;
  fub::Complete<fub::PerfectGas<2>> left;
  fub::CompleteFromCons(equation, left, cons);

  fub::amrex::cutcell::RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.015),
                              left, right);

  using State = fub::Complete<fub::PerfectGas<2>>;
  fub::amrex::cutcell::GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

  using fub::amrex::cutcell::TransmissiveBoundary;
  fub::amrex::cutcell::BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1}}};

  using fub::amrex::cutcell::GriddingAlgorithm;
  using fub::amrex::cutcell::PatchHierarchy;
  using fub::amrex::cutcell::TagCutCells;
  using fub::amrex::cutcell::TagAllOf;
  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, hier_opts), initial_data,
      TagAllOf(TagCutCells(), gradients), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  using fub::amrex::cutcell::TimeIntegrator;
  using fub::amrex::cutcell::Reconstruction;
  fub::amrex::cutcell::HyperbolicMethod method{flux_method,
                          TimeIntegrator{},
                          Reconstruction{equation}};

  const int base_gcw = flux_method.GetStencilWidth();
  const int scratch_gcw = base_gcw * 4;
  const int flux_gcw = base_gcw * 3;
  using fub::amrex::cutcell::IntegratorContext;
  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method, scratch_gcw, flux_gcw),
      fub::StrangSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
//   fub::NoSubcycleSolver solver(std::move(level_integrator));

  std::string base_name = "Schardin/";
  using namespace std::literals::chrono_literals;
  using Plotfile = fub::amrex::cutcell::PlotfileOutput<fub::PerfectGas<2>>;
  using CounterOutput = fub::CounterOutput<GriddingAlgorithm,
                                           std::chrono::milliseconds>;
  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<Plotfile>("Plotfile", equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<fub::amrex::cutcell::WriteHdf5>("HDF5");
  fub::MultipleOutputs<GriddingAlgorithm> output(
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