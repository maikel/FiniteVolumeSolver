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

#include "fub/cutcell_method/MyStabilisation.hpp"

#include "fub/AMReX/cutcell/boundary_condition/ConstantBoundary.hpp"

#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

static_assert(AMREX_SPACEDIM == 2);

using Coord = Eigen::Vector2d;

static double constant_function(double, double rho_0, double) noexcept {
  return rho_0;
}

static double initial_function1(double rel_x, double rho_0,
                                double width) noexcept {
  return std::min(std::max(rho_0 + rel_x, rho_0), rho_0 + width);
}

static double initial_function2(double rel_x, double rho_0,
                                double width) noexcept {
  const double exponent = 2.0 * std::abs(rel_x) / width;
  const double exponent2 = exponent * exponent;
  const double rho = rho_0 + std::exp(-exponent2);
  return rho;
}

struct WaveFunction {
  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::cutcell::GriddingAlgorithm& grid,
                      int level, fub::Duration /*time*/) const {
    ::amrex::MultiFab& data = patch_level.data;
    const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    const std::shared_ptr<::amrex::EBFArrayBoxFactory>& factory =
        grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
    const ::amrex::MultiFab& volfrac = factory->getVolFrac();
    const auto& facecent = factory->getFaceCent();
    const auto& bdrycent = factory->getBndryCent();
    const auto& bdrynorm = factory->getBndryNormal();
    FUB_ASSERT(facecent[0] && facecent[1]);
    fub::amrex::ForEachFab(
        fub::execution::openmp, data, [&](const ::amrex::MFIter& mfi) {
          ::amrex::FArrayBox& fab = data[mfi];
          const ::amrex::FArrayBox& alpha = volfrac[mfi];
          const auto& fx = (*facecent[0])[mfi];
          const auto& fy = (*facecent[1])[mfi];
          const auto& fB = bdrycent[mfi];
          const auto& bn = bdrynorm[mfi];
          ::amrex::Box box = mfi.tilebox();
          auto states = fub::amrex::MakeView<fub::Complete<fub::PerfectGas<2>>>(
              fab, equation_, box);
          fub::amrex::ForEachIndex(box, [&](std::ptrdiff_t i,
                                            std::ptrdiff_t j) {
            double xhi[2];
            double xlo[2];
            double yhi[2];
            double ylo[2];
            const amrex::IntVect iv{int(i), int(j)};
            geom.LoFace(iv, 0, xlo);
            geom.HiFace(iv, 0, xhi);
            geom.LoFace(iv, 1, ylo);
            geom.HiFace(iv, 1, yhi);
            // Coord xhi(geom.CellCenter(i, 0), geom.CellCenter(j, 1));
            fub::Array<double, 2, 1> u = u_0_ * direction_;
            double p = p_0_;
            double rho = rho_0_;
            if (alpha(iv) == 1.0) {
              Coord x1{xlo[0], xlo[1]};
              Coord x2{xhi[0], xhi[1]};
              Coord x3{ylo[0], ylo[1]};
              Coord x4{yhi[0], yhi[1]};
              const double relative_x1 = (x1 - origin_).dot(direction_);
              const double relative_x2 = (x2 - origin_).dot(direction_);
              const double relative_x3 = (x3 - origin_).dot(direction_);
              const double relative_x4 = (x4 - origin_).dot(direction_);
              const double rho1 = initial_function(relative_x1, rho_0_, width_);
              const double rho2 = initial_function(relative_x2, rho_0_, width_);
              const double rho3 = initial_function(relative_x3, rho_0_, width_);
              const double rho4 = initial_function(relative_x4, rho_0_, width_);
              rho = (rho1 + rho2 + rho3 + rho4) / 4.0;
              fub::Complete<fub::PerfectGas<2>> state =
                  equation_.CompleteFromPrim(rho, u, p);
              fub::Store(states, state, {i, j});
            } else if (alpha(iv) > 0.0) {
              geom.LoFace(iv, 0, xlo);
              geom.HiFace(iv, 0, xhi);
              geom.LoFace(iv, 1, ylo);
              geom.HiFace(iv, 1, yhi);
              amrex::IntVect ivxR = iv;
              ivxR.shift({1, 0});
              amrex::IntVect ivyR = iv;
              ivyR.shift({0, 1});
              double xlo_offset = fx(iv);
              double xhi_offset = fx(ivxR);
              double ylo_offset = fy(iv);
              double yhi_offset = fy(ivyR);
              double xb_offset = fB(iv, 0);
              double yb_offset = fB(iv, 1);
              const double dx = geom.CellSize(0);
              const double dy = geom.CellSize(1);
              Coord x1{xlo[0], xlo[1] + xlo_offset * dy};
              Coord x2{xhi[0], xhi[1] + xhi_offset * dy};
              Coord x3{ylo[0] + ylo_offset * dx, ylo[1]};
              Coord x4{yhi[0] + yhi_offset * dx, yhi[1]};
              Coord x5{geom.CellCenter(i, 0) + xb_offset * dx,
                       geom.CellCenter(j, 1) + yb_offset * dy};
              Coord n{bn(iv, 0), bn(iv, 1)};

              std::array<Coord, 5> xs{x1, x2, x3, x4, x5};
              std::array<double, 5> rel_x{};
              int rel_n = 0;
              for (int i = 0; i < 5; ++i) {
                if ((xs[i] - x5).dot(n) <= 0.0) {
                  rel_x[rel_n] = (xs[i] - origin_).dot(direction_);
                  rel_n += 1;
                }
              }
              rho = 0.0;
              for (int i = 0; i < rel_n; ++i) {
                rho +=
                    initial_function(rel_x[i], rho_0_, width_) / double(rel_n);
              }
              fub::Complete<fub::PerfectGas<2>> state =
                  equation_.CompleteFromPrim(rho, u, p);
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
  std::function<double(double, double, double)> initial_function{
      &initial_function2};
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
    fub::MyCutCellMethod<fub::PerfectGas<2>, FluxMethod> cutcell_method(eq);
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
  const double theta = fub::GetOptionOr(opts, "theta", M_PI * 30.0 / 180.0);
  const double W = 0.0141;

  using std::cos;
  using std::sin;
  Eigen::Vector2d p0{0.0, 0.0};
  Eigen::Vector2d p1{cos(theta), sin(theta)};
  Eigen::Vector2d q0{p0[0] - W * sin(theta), p0[1] + W * cos(theta)};
  Eigen::Vector2d q1{p1[0] - W * sin(theta), p1[1] + W * cos(theta)};

  auto embedded_boundary =
      ::amrex::EB2::makeUnion(Plane(p0, p1), Plane(q1, q0));
  auto shop = amrex::EB2::makeShop(embedded_boundary);
  hier_opts.index_spaces = MakeIndexSpaces(shop, geometry, hier_opts);

  using namespace std::literals;
  // using HLLE =
  //     fub::HllMethod<fub::PerfectGas<2>,
  //                    fub::EinfeldtSignalVelocities<fub::PerfectGas<2>>>;
  using HLLEM = fub::perfect_gas::HllemMethod<fub::PerfectGas<2>>;
  // using ConservativeReconstruction =
  //     fub::MusclHancockMethod<fub::PerfectGas<2>, HLLE, fub::VanLeer>;
  using ConservativeReconstructionNoGradient =
      fub::MusclHancockMethod<fub::PerfectGas<2>, HLLEM, fub::NoGradient>;
  using ConservativeReconstructionNoLimiter =
      fub::MusclHancockMethod<fub::PerfectGas<2>, HLLEM, fub::NoLimiter>;
  using ConservativeReconstructionMinMod =
      fub::MusclHancockMethod<fub::PerfectGas<2>, HLLEM, fub::MinMod>;
  using ConservativeReconstructionVanLeer =
      fub::MusclHancockMethod<fub::PerfectGas<2>, HLLEM, fub::VanLeer>;
  // using PrimitiveReconstruction =
  //     fub::FluxMethod<fub::perfect_gas::MusclHancockPrim<2>>;
  // using CharacteristicReconstruction =
  //     fub::perfect_gas::MusclHancockCharMethod<2>;

  auto flux_method_factory = GetFluxMethodFactory(
      // std::pair{"HLLE"s, MakeFlux<HLLE>()},
      // std::pair{"HLLEM"s, MakeFlux<HLLEM>()},
      // std::pair{"Primitive"s, MakeFlux<PrimitiveReconstruction>()},
      // std::pair{"Conservative"s, MakeFlux<ConservativeReconstruction>()},
      std::pair{"ConservativeNoGradient"s,
                MakeFlux<ConservativeReconstructionNoGradient>()},
      std::pair{"ConservativeNoLimiter"s,
                MakeFlux<ConservativeReconstructionNoLimiter>()},
      std::pair{"ConservativeVanLeer"s,
                MakeFlux<ConservativeReconstructionVanLeer>()},
      std::pair{"ConservativeMinMod"s,
                MakeFlux<ConservativeReconstructionMinMod>()});
  // std::pair{"Characteristics"s, MakeFlux<CharacteristicReconstruction>()});

  std::string reconstruction =
      fub::GetOptionOr(opts, "reconstruction", "Characteristics"s);
  BOOST_LOG(log) << "Reconstruction: " << reconstruction;
  auto flux_method = flux_method_factory.at(reconstruction)(equation);

  const double relative_origin = fub::GetOptionOr(opts, "origin", 0.035);
  const Coord origin{relative_origin * cos(theta),
                     relative_origin * sin(theta)};
  const Coord direction{cos(theta), sin(theta)};
  WaveFunction initial_data{equation, origin, direction};
  std::string initial_function =
      fub::GetOptionOr(opts, "initial_function", "Linear"s);
  if (initial_function == "Linear") {
    initial_data.initial_function = &initial_function1;
  } else if (initial_function == "Constant") {
    initial_data.initial_function = &constant_function;
  }
  double jump = fub::GetOptionOr(opts, "initial_data_jump", 0.0);
  BOOST_LOG(log) << "Initial Function: " << initial_function;
  BOOST_LOG(log) << "Initial Data Jump: " << jump;

  fub::Array<double, 2, 1> u = initial_data.u_0_ * initial_data.direction_;
  fub::Complete<fub::PerfectGas<2>> stateL =
      equation.CompleteFromPrim(initial_data.rho_0_, u, initial_data.p_0_);

  fub::Complete<fub::PerfectGas<2>> stateR = equation.CompleteFromPrim(
      initial_data.rho_0_ + jump, u, initial_data.p_0_);

  using fub::amrex::cutcell::ConstantBoundary;
  using fub::amrex::cutcell::ReflectiveBoundary;
  fub::amrex::cutcell::BoundarySet boundary_condition{
      {ConstantBoundary<fub::PerfectGas<2>>{fub::Direction::X, 0, equation,
                                            stateL},
       ConstantBoundary<fub::PerfectGas<2>>{fub::Direction::X, 1, equation,
                                            stateR},
       ConstantBoundary<fub::PerfectGas<2>>{fub::Direction::Y, 0, equation,
                                            stateL},
       ConstantBoundary<fub::PerfectGas<2>>{fub::Direction::Y, 1, equation,
                                            stateR}}};
  // auto seq = fub::execution::seq;
  // fub::amrex::cutcell::BoundarySet boundary_condition{
  //     {ReflectiveBoundary{seq, equation, fub::Direction::X, 0},
  //      ReflectiveBoundary{seq, equation, fub::Direction::X, 1},
  //      ReflectiveBoundary{seq, equation, fub::Direction::Y, 0},
  //      ReflectiveBoundary{seq, equation, fub::Direction::Y, 1}}};

  using fub::amrex::cutcell::GriddingAlgorithm;
  using fub::amrex::cutcell::PatchHierarchy;
  using fub::amrex::cutcell::TagAllOf;
  using fub::amrex::cutcell::TagBuffer;
  using fub::amrex::cutcell::TagCutCells;
  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, hier_opts), initial_data,
      TagAllOf(TagCutCells()), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  using fub::amrex::cutcell::Reconstruction;
  using fub::amrex::cutcell::TimeIntegrator2;
  fub::amrex::cutcell::HyperbolicMethod method{flux_method, TimeIntegrator2{},
                                               Reconstruction{equation}};

  const int base_gcw = flux_method.GetStencilWidth();
  const int scratch_gcw = base_gcw + 1;
  const int flux_gcw = 0;
  using fub::amrex::cutcell::IntegratorContext;
  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>, IntegratorContext(gridding, method, scratch_gcw, flux_gcw),
      fub::GodunovSplitting());

  // fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
  fub::NoSubcycleSolver solver(std::move(level_integrator));

  using namespace std::literals::chrono_literals;
  using Plotfile = fub::amrex::cutcell::PlotfileOutput<fub::PerfectGas<2>>;
  using CounterOutput =
      fub::CounterOutput<GriddingAlgorithm, std::chrono::milliseconds>;
  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<Plotfile>("Plotfile", equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<fub::amrex::cutcell::WriteHdf5>("HDF5");
  factory.RegisterOutput<fub::amrex::cutcell::DebugOutput>(
      "DebugOutput",
      solver.GetGriddingAlgorithm()->GetPatchHierarchy().GetDebugStorage());
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