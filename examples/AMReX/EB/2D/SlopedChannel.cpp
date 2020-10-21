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

#include "fub/AMReX/cutcell/boundary_condition/ConstantBoundary.hpp"

#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Plane.H>


static_assert(AMREX_SPACEDIM == 2);

using Coord = Eigen::Vector2d;

struct WaveFunction {
  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::cutcell::GriddingAlgorithm& grid, int level,
                      fub::Duration /*time*/) const {
    ::amrex::MultiFab& data = patch_level.data;
    const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    const std::shared_ptr<::amrex::EBFArrayBoxFactory>& factory =
        grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
    const ::amrex::MultiFab& volfrac = factory->getVolFrac();
    fub::amrex::ForEachFab(fub::execution::openmp, data, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& fab = data[mfi];
      const ::amrex::FArrayBox& alpha = volfrac[mfi];
      ::amrex::Box box = mfi.tilebox();
      auto states = fub::amrex::MakeView<fub::Complete<fub::PerfectGas<2>>>(fab, equation_, box);
      fub::amrex::ForEachIndex(box, [&](std::ptrdiff_t i, std::ptrdiff_t j) {
        Coord x(geom.CellCenter(i, 0), geom.CellCenter(j, 1));
        const ::amrex::IntVect iv{int(i), int(j)};
        if (alpha(iv) > 0.0) {
          const double relative_x = (x - origin_).dot(direction_);
          const double exponent = 2.0 * std::abs(relative_x) / width_;
          const double exponent2 = exponent * exponent;
          double rho = rho_0_ + std::exp(-exponent2);
          fub::Array<double, 2, 1> u = u_0_ * direction_;
          double p = p_0_;
          fub::Complete<fub::PerfectGas<2>> state = equation_.CompleteFromPrim(rho, u, p);
          fub::Store(states, state, {i, j});
        } else {
          for (int comp = 0; comp < fab.nComp(); ++comp) {
            fab(iv, comp) = 0.0;
          }
        }
      });
    });
  }

  fub::PerfectGas<2> equation_;
  Coord origin_;
  Coord direction_;
  double rho_0_{1.2};
  double u_0_{30.0};
  double p_0_{101325.0};
  double width_{0.0141};
};

Coord OrthogonalTo(const Coord& x) { return Coord{x[1], -x[0]}; }

auto Plane(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1) {
  Eigen::Vector2d norm1 = OrthogonalTo(p1 - p0).normalized();
  amrex::EB2::PlaneIF plane1({p0[0], p0[1]}, {norm1[0], norm1[1]}, false);
  return amrex::EB2::makeComplement(plane1);
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
  const double theta = M_PI * 30.0 / 180.0;
  const double W = 0.0141;

  using std::sin;
  using std::cos;
  Eigen::Vector2d p0{0.0, 0.0};
  Eigen::Vector2d p1{cos(theta), sin(theta)};
  Eigen::Vector2d q0{p0[0] - W * sin(theta), p0[1] + W * cos(theta)};
  Eigen::Vector2d q1{p1[0] - W * sin(theta), p1[1] + W * cos(theta)};

  auto embedded_boundary = ::amrex::EB2::makeUnion(Plane(p0, p1), Plane(q1, q0));
  auto shop = amrex::EB2::makeShop(embedded_boundary);
  hier_opts.index_spaces =
      MakeIndexSpaces(shop, geometry, hier_opts);

  using Complete = fub::PerfectGas<2>::Complete;
  fub::amrex::cutcell::GradientDetector gradient(equation,
                                        std::pair{&Complete::density, 1.0e-2});

  using namespace std::literals;
  using HLLE = fub::HllMethod<fub::PerfectGas<2>, fub::EinfeldtSignalVelocities<fub::PerfectGas<2>>>;
  using HLLEM = fub::perfect_gas::HllemMethod<2>;
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

  const double relative_origin = 0.035;
  const Coord origin{relative_origin * cos(theta), relative_origin * sin(theta)};
  const Coord direction{cos(theta), sin(theta)};
  WaveFunction initial_data{equation, origin, direction};

  fub::Array<double, 2, 1> u = initial_data.u_0_ * initial_data.direction_;
  fub::Complete<fub::PerfectGas<2>> state = equation.CompleteFromPrim(initial_data.rho_0_, u, initial_data.p_0_);

  using State = fub::Complete<fub::PerfectGas<2>>;
  fub::amrex::cutcell::GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

  using fub::amrex::cutcell::ConstantBoundary;
  fub::amrex::cutcell::BoundarySet boundary_condition{
                                 {  ConstantBoundary<fub::PerfectGas<2>>{fub::Direction::X, 0, equation, state},
                                    ConstantBoundary<fub::PerfectGas<2>>{fub::Direction::X, 1, equation, state},
                                    ConstantBoundary<fub::PerfectGas<2>>{fub::Direction::Y, 0, equation, state},
                                    ConstantBoundary<fub::PerfectGas<2>>{fub::Direction::Y, 1, equation, state}}};

  using fub::amrex::cutcell::GriddingAlgorithm;
  using fub::amrex::cutcell::PatchHierarchy;
  using fub::amrex::cutcell::TagCutCells;
  using fub::amrex::cutcell::TagAllOf;
  using fub::amrex::cutcell::TagBuffer;
  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, hier_opts), initial_data,
      TagAllOf(TagCutCells()), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  using fub::amrex::cutcell::TimeIntegrator;
  using fub::amrex::cutcell::Reconstruction;
  fub::amrex::cutcell::HyperbolicMethod method{flux_method,
                          TimeIntegrator{},
                          Reconstruction{equation}};

  const int base_gcw = flux_method.GetStencilWidth();
  const int scratch_gcw = base_gcw * 2;
  const int flux_gcw = base_gcw * 1;
  using fub::amrex::cutcell::IntegratorContext;
  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method, scratch_gcw, flux_gcw),
      fub::StrangSplitting());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
//   fub::NoSubcycleSolver solver(std::move(level_integrator));

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