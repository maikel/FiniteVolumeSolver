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

#include "fub/flux_method/MusclHancockMethod2.hpp"

#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>

#include <iostream>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>

#include <boost/container/pmr/vector.hpp>
#include <boost/filesystem.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using HLLEM = fub::perfect_gas::HllemMethod<fub::PerfectGas<2>>;

using ConservativeReconstruction = fub::FluxMethod<fub::MusclHancock2<
    fub::PerfectGas<2>,
    fub::ConservativeGradient<
        fub::PerfectGas<2>,
        fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
    fub::ConservativeReconstruction<fub::PerfectGas<2>>, HLLEM>>;

using PrimitiveReconstruction = fub::FluxMethod<fub::MusclHancock2<
    fub::PerfectGas<2>,
    fub::PrimitiveGradient<fub::PerfectGas<2>, fub::CentralDifferenceGradient<
                                                    fub::VanLeerLimiter>>,
    fub::PrimitiveReconstruction<fub::PerfectGas<2>>, HLLEM>>;

using CharacteristicsReconstruction = fub::FluxMethod<fub::MusclHancock2<
    fub::PerfectGas<2>,
    fub::CharacteristicsGradient<
        fub::PerfectGas<2>,
        fub::CentralDifferenceGradient<fub::VanLeerLimiter>>,
    fub::CharacteristicsReconstruction<fub::PerfectGas<2>>, HLLEM>>;

fub::Polygon ReadPolygonData(std::istream& input) {
  std::string line{};
  namespace pmr = boost::container::pmr;
  pmr::vector<double> xs{};
  pmr::vector<double> ys{};
  while (std::getline(input, line)) {
    double x{}, y{};
    std::istringstream linestream(line);
    linestream >> x >> y;
    if (linestream) {
      xs.push_back(x);
      ys.push_back(y);
    }
  }
  if (xs.empty() || ys.empty()) {
    throw std::invalid_argument{"Invalid Input File:: No data found!"};
  }
  const double x0 = xs.front();
  const double x1 = xs.back();
  const double y0 = ys.front();
  const double y1 = ys.back();
  if (x0 != x1 || y0 != y1) {
    throw std::invalid_argument{
        "Invalid Input File: First and last entries are not the same point."};
  }
  return fub::Polygon(std::move(xs), std::move(ys));
}

template <typename Geometry>
struct ShockMachnumber
    : fub::amrex::cutcell::RiemannProblem<fub::PerfectGas<2>, Geometry> {
  static fub::PerfectGas<2>::Complete
  ComputePreShockState(fub::PerfectGas<2>& equation,
                       const fub::PerfectGas<2>::Complete& post_shock,
                       double M_S, const fub::Array<double, 2, 1>& normal) {
    const double g = equation.gamma;
    const double gp1 = g + 1.0;
    const double gm1 = g - 1.0;
    const double M_post = equation.Machnumber(post_shock);
    const double M2 = (M_post - M_S) * (M_post - M_S);

    const double rho_c = gp1 * M2 / (gm1 * M2 + 2.0);
    const double rho_post = post_shock.density;
    const double rho_pre = rho_post * rho_c;

    const double p_c = (2.0 * g * M2 - gm1) / gp1;
    const double p_post = post_shock.pressure;
    const double p_pre = p_post * p_c;

    const double t = rho_post / rho_pre;
    const double a_post = post_shock.speed_of_sound;
    const double u_S = M_S * a_post;
    const double u_post = equation.Velocity(post_shock).matrix().norm();
    const double u_pre = u_S * (1.0 - t) + u_post * t;

    return equation.CompleteFromPrim(rho_pre, u_pre * normal, p_pre);
  }

  ShockMachnumber(fub::PerfectGas<2> equation, Geometry geometry,
                  const fub::PerfectGas<2>::Complete& post_shock, double M_S,
                  const fub::Array<double, 2, 1>& normal)
      : fub::amrex::cutcell::RiemannProblem<fub::PerfectGas<2>, Geometry>(
            equation, std::move(geometry),
            ComputePreShockState(equation, post_shock, M_S, normal),
            post_shock) {}
};

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
    HLLEM hllem{eq};
    FluxMethod flux_method{eq};
    fub::KbnCutCellMethod cutcell_method(flux_method, hllem);
    fub::amrex::cutcell::FluxMethod adapter(std::move(cutcell_method));
    return adapter;
  }
};

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};
  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::PerfectGas<2> equation;

  fub::ProgramOptions equation_options = fub::GetOptions(options, "Equation");
  equation.gamma = fub::GetOptionOr(equation_options, "gamma", equation.gamma);
  equation.Rspec = fub::GetOptionOr(equation_options, "Rpsec", equation.Rspec);
  equation.gamma_minus_1_inv = 1.0 / (equation.gamma - 1.0);
  equation.gamma_array_ = fub::Array1d::Constant(equation.gamma);
  equation.gamma_minus_1_inv_array_ = fub::Array1d::Constant(equation.gamma_minus_1_inv);

  BOOST_LOG(log) << "Equation:";
  BOOST_LOG(log) << fmt::format(" - gamma = {}", equation.gamma);
  BOOST_LOG(log) << fmt::format(" - Rspec = {}", equation.Rspec);

  fub::amrex::CartesianGridGeometry geometry =
      fub::GetOptions(options, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  geometry.Print(log);

  const int scratch_gcw = fub::GetOptionOr(options, "scratch_gcw", 4);
  const int flux_gcw = fub::GetOptionOr(options, "flux_gcw", 2);
  fub::amrex::cutcell::PatchHierarchyOptions hier_opts =
      fub::GetOptions(options, "PatchHierarchy");
  hier_opts.ngrow_eb_level_set =
      std::max(scratch_gcw + 1, hier_opts.ngrow_eb_level_set);
  BOOST_LOG(log) << "PatchHierarchy:";
  hier_opts.Print(log);

  BOOST_LOG(log) << "Read Wall Data...";
  std::vector<std::string> wall_filenames{};
  wall_filenames = fub::GetOptionOr(options, "wall_filenames", wall_filenames);

  std::vector<fub::PolymorphicGeometry<2>> geometries;
  std::transform(wall_filenames.begin(), wall_filenames.end(),
                 std::back_inserter(geometries),
                 [](const std::string& filename) {
                   std::ifstream ifs(filename);
                   return fub::PolymorphicGeometry<2>(ReadPolygonData(ifs));
                 });
  BOOST_LOG(log) << "Compute EB level set data...";
  auto embedded_boundary =
      fub::amrex::Geometry(fub::PolymorphicUnion(geometries));
  auto shop = amrex::EB2::makeShop(embedded_boundary);
  hier_opts.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, geometry, hier_opts);

  using namespace fub::amrex::cutcell;
  using State = fub::Complete<fub::PerfectGas<2>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.22;
  cons.momentum << 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> post_shock_state;
  fub::CompleteFromCons(equation, post_shock_state, cons);

  double shock_mach_number = fub::GetOptionOr(options, "Mach_number", 1.1);
  shock_mach_number =
      fub::GetOptionOr(options, "shock_mach_number", shock_mach_number);
  const fub::Array<double, 2, 1> normal{1.0, 0.0};

  ShockMachnumber<fub::Halfspace> initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.01), post_shock_state,
      shock_mach_number, normal);

  const State& pre_shock_state = initial_data.left_;
  BOOST_LOG(log) << "Post-Shock-State:\n"
                 << "\tdensity: " << post_shock_state.density << " [kg/m^3]\n"
                 << "\tvelocity: "
                 << equation.Velocity(post_shock_state).transpose()
                 << " [m/s]\n"
                 << "\tpressure: " << post_shock_state.pressure << " [Pa]";

  BOOST_LOG(log) << "Calculated Pre-Shock-State:\n"
                 << "\tdensity: " << pre_shock_state.density << " [kg/m^3]\n"
                 << "\tvelocity: "
                 << equation.Velocity(pre_shock_state).transpose() << " [m/s]\n"
                 << "\tpressure: " << pre_shock_state.pressure << " [Pa]";

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, hier_opts), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  using namespace std::literals;
  auto flux_method_factory = GetFluxMethodFactory(
      std::pair{"NoReconstruct"s, MakeFlux<HLLEM>()},
      std::pair{"Conservative"s, MakeFlux<ConservativeReconstruction>()},
      std::pair{"Primitive"s, MakeFlux<PrimitiveReconstruction>()},
      std::pair{"Characteristics"s, MakeFlux<CharacteristicsReconstruction>()});

  std::string reconstruction =
      fub::GetOptionOr(options, "reconstruction", "Characteristics"s);
  BOOST_LOG(log) << "Reconstruction: " << reconstruction;
  auto flux_method = flux_method_factory.at(reconstruction)(equation);

  HyperbolicMethod method{flux_method, TimeIntegrator{},
                          Reconstruction{equation}};

  BOOST_LOG(log) << fmt::format("scratch_gcw = {}", scratch_gcw);
  BOOST_LOG(log) << fmt::format("flux_gcw = {}", flux_gcw);
  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<2>,
      IntegratorContext(gridding, method, scratch_gcw, flux_gcw));

  // fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
  fub::NoSubcycleSolver solver(std::move(level_integrator));

  using CounterOutput = fub::CounterOutput<GriddingAlgorithm,
                                           std::chrono::milliseconds>;
  using Plotfile = PlotfileOutput<fub::PerfectGas<2>>;
  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<WriteHdf5>("HDF5");
  factory.RegisterOutput<Plotfile>("Plotfiles", equation);
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
