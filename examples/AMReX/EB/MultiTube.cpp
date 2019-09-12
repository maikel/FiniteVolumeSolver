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
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <boost/filesystem.hpp>

#include <cmath>
#include <iostream>

#include <xmmintrin.h>

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 3;

static constexpr double r_tube = 0.015;
static constexpr double r_inner = 0.5 * 0.130;
static constexpr double r_outer = 0.5 * 0.385;
static constexpr double r_tube_center = 0.5 * r_inner + 0.5 * r_outer;
static constexpr double alpha = 2. * M_PI / 5.;

auto Center(double x, double phi) -> ::amrex::RealArray {
  using std::cos;
  using std::sin;
  return {x, r_tube_center * sin(phi), r_tube_center * cos(phi)};
}

auto DomainAroundCenter(const ::amrex::RealArray& x, double rx)
    -> ::amrex::RealBox {
  return ::amrex::RealBox{{x[0] - rx, x[1] - r_tube, x[2] - r_tube},
                          {x[0] + rx, x[1] + r_tube, x[2] + r_tube}};
}

struct TubeSolverOptions {
  int n_cells{200};
  int max_refinement_level{1};
  std::array<double, 2> x_domain{-1.5, -0.03};
  double phi{0.0};
};

auto MakeTubeSolver(fub::Burke2012& mechanism, const TubeSolverOptions& opts,
                    const boost::program_options::variables_map& vm) {
  const std::array<int, AMREX_SPACEDIM> n_cells{opts.n_cells, 1, 1};
  const double x_lo = opts.x_domain[0];
  const double x_up = opts.x_domain[1];
  const double x_len = x_up - x_lo;
  const double r_len = 0.5 * x_len;
  const double x_mid = 0.5 * x_lo + 0.5 * x_up;
  amrex::RealBox xbox = DomainAroundCenter(Center(x_mid, opts.phi), r_len);
  const std::array<int, AMREX_SPACEDIM> periodicity{0, 0, 0};

  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};

  using namespace fub::amrex;

  CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = xbox;

  DataDescription desc = MakeDataDescription(equation);

  PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = opts.max_refinement_level;
  hier_opts.refine_ratio = amrex::IntVect{2, 1, 1};

  amrex::Geometry geom(
      amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1}}, &xbox,
      -1, periodicity.data());
  geom.refine(hier_opts.refine_ratio);

  using Complete = fub::IdealGasMix<1>::Complete;
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 1e-3),
                            std::make_pair(&Complete::pressure, 1e-2),
                            std::make_pair(&Complete::temperature, 1e-1)};

  ::amrex::Box refine_box{{opts.n_cells - 5, 0, 0}, {opts.n_cells - 1, 0, 0}};
  ConstantBox constant_box{refine_box};

  equation.GetReactor().SetMoleFractions("N2:79,O2:21");
  equation.GetReactor().SetTemperature(300.0);
  equation.GetReactor().SetPressure(101325.0);
  fub::Complete<fub::IdealGasMix<Tube_Rank>> state(equation);
  equation.CompleteFromReactor(state);
  ConstantData initial_data{equation, state};

  PressureValveBoundary valve{equation, vm};
  BoundarySet boundaries{{valve}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(desc, geometry, hier_opts), initial_data,
      TagAllOf(gradient, constant_box), boundaries);
  gridding->InitializeHierarchy(0.0);

  fub::ideal_gas::MusclHancockPrimMethod<Tube_Rank> flux_method{equation};
  HyperbolicMethod method{FluxMethod(fub::execution::openmp, flux_method),
                          ForwardIntegrator(fub::execution::openmp),
                          Reconstruction(fub::execution::openmp, equation)};

  return fub::amrex::IntegratorContext(gridding, method);
}

::amrex::Box BoxWhichContains(const ::amrex::RealBox& xbox,
                              const ::amrex::Geometry& geom) {
  ::amrex::Box domain = geom.Domain();
  ::amrex::IntVect lo = domain.smallEnd();
  ::amrex::IntVect up = domain.bigEnd();
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    for (int i = domain.smallEnd(d); i < domain.bigEnd(d); ++i) {
      const double x = geom.CellCenter(i, d);
      if (x < xbox.lo(d)) {
        lo[d] = std::max(lo[d], i);
      }
      if (x > xbox.hi(d)) {
        up[d] = std::min(up[d], i);
      }
    }
  }
  return ::amrex::Box{lo, up};
}

auto MakePlenumSolver(int num_cells, int n_level, fub::Burke2012& mechanism) {
  const std::array<int, Plenum_Rank> n_cells{num_cells, num_cells, num_cells};
  const std::array<double, Plenum_Rank> xlower{-0.03, -0.5 * 0.56, -0.5 * 0.56};
  const std::array<double, Plenum_Rank> xupper{+0.53, +0.5 * 0.56, +0.5 * 0.56};
  const std::array<int, Plenum_Rank> periodicity{0, 0, 0};

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(amrex::Box{{}, ::amrex::IntVect(num_cells - 1)},
                              &xbox, -1, periodicity.data());

  auto embedded_boundary = amrex::EB2::makeUnion(
      amrex::EB2::makeIntersection(
          amrex::EB2::CylinderIF(r_outer, 0.5, 0, {0.25, 0.0, 0.0}, true),
          amrex::EB2::CylinderIF(r_tube_center, 1.0, 0,
                                 {1.0 - 1.0e-6, 0.0, 0.0}, true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 0.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 1.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 2.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 3.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 4.0 * alpha),
                                 true)),
      amrex::EB2::CylinderIF(r_inner, 1.0, 0, {0.25, 0.0, 0.0}, false));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};

  // Make Gridding Algorithm

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  equation.GetReactor().SetMoleFractions("N2:79,O2:21");
  equation.GetReactor().SetTemperature(300.0);
  equation.GetReactor().SetPressure(101325.0);
  fub::Complete<fub::IdealGasMix<Plenum_Rank>> right(equation);
  equation.CompleteFromReactor(right);

  using namespace fub::amrex::cutcell;

  fub::amrex::cutcell::RiemannProblem initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04), right, right);

  PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces = MakeIndexSpaces(shop, coarse_geom, n_level);

  using State = fub::Complete<fub::IdealGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

  ::amrex::RealBox inlet{{-0.1, -0.015, -0.015}, {0.05, +0.015, +0.015}};
  const ::amrex::Box refine_box = BoxWhichContains(inlet, coarse_geom);
  ConstantBox constant_box{refine_box};

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1},
                                  TransmissiveBoundary{fub::Direction::Z, 0},
                                  TransmissiveBoundary{fub::Direction::Z, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, options), initial_data,
      TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
      boundary_condition);
  gridding->InitializeHierarchy(0.0);

  // Make Solver

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::ideal_gas::MusclHancockPrimMethod<Plenum_Rank> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{fub::execution::openmp, cutcell_method},
                          fub::amrex::cutcell::TimeIntegrator{},
                          Reconstruction{fub::execution::openmp, equation}};

  return fub::amrex::cutcell::IntegratorContext(gridding, method);
}

struct ProgramOptions {
  ProgramOptions() = default;

  explicit ProgramOptions(const boost::program_options::variables_map& vm) {
    plenum_n_cells = vm["plenum_n_cells"].as<int>();
    max_refinement_level = vm["max_number_of_levels"].as<int>();
    plenum_checkpoint = vm["plenum_checkpoint"].as<std::string>();
    tube_checkpoint = vm["tube_checkpoint"].as<std::string>();
  }

  static boost::program_options::options_description GetCommandLineOptions() {
    namespace po = boost::program_options;
    po::options_description desc{};
    // clang-format off
    desc.add_options()
        ("plenum_n_cells", po::value<int>()->default_value(128), "Base number of cells in the plenum for the coarsest level")
        ("max_number_of_levels", po::value<int>()->default_value(1), "Maximal number of refinement levels across all domains.")
        ("plenum_checkpoint", po::value<std::string>()->default_value(""), "The path to the checkpoint files for the plenum.")
        ("tube_checkpoint", po::value<std::string>()->default_value(""), "The path to the checkpoint files for the tubes.");
    // clang-format on
    return desc;
  }

  int plenum_n_cells{128};
  int max_refinement_level{1};
  std::string plenum_checkpoint{};
  std::string tube_checkpoint{};
};

void MyMain(const ProgramOptions& po,
            const boost::program_options::variables_map& vm);

std::optional<boost::program_options::variables_map>
ParseCommandLine(int argc, char** argv) {
  namespace po = boost::program_options;
  po::options_description desc = fub::GetCommandLineRunOptions();
  desc.add(ProgramOptions::GetCommandLineOptions());
  desc.add(fub::amrex::PressureValveOptions::GetCommandLineOptions());
  desc.add(fub::amrex::IgniteDetonationOptions::GetCommandLineOptions());
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (std::exception& e) {
    amrex::Print()
        << "[Error] An Error occured while reading program options:\n";
    amrex::Print() << e.what();
    return {};
  }

  if (vm.count("help")) {
    amrex::Print() << desc << "\n";
    return {};
  }

  std::ostringstream sout{};
  fub::PrintRunOptions(sout, fub::GetRunOptions(vm));
  amrex::Print() << sout.str();

  ProgramOptions o(vm);
  constexpr double tube_len_over_plenum_len = 1.47 / 0.56;
  int tube_n_cells = static_cast<int>(tube_len_over_plenum_len *
                                      static_cast<double>(o.plenum_n_cells));
  tube_n_cells = tube_n_cells - tube_n_cells % 8;
  ::amrex::Print() << "[Info] plenum_n_cells = " << o.plenum_n_cells << '\n';
  ::amrex::Print() << "[Info] tube_n_cells = " << tube_n_cells << '\n';
  ::amrex::Print() << "[Info] max_refinement_level = " << o.max_refinement_level
                   << '\n';

  if (!o.plenum_checkpoint.empty()) {
    ::amrex::Print() << "[Info] plenum_checkpoint = " << o.plenum_checkpoint
                     << '\n';
    ::amrex::Print() << "[Info] tube_checkpoint = " << o.tube_checkpoint
                     << '\n';
  }
  return vm;
}

void MyMain(const boost::program_options::variables_map& vm) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::Burke2012 mechanism{};

  ProgramOptions po(vm);

  const int n_level = po.max_refinement_level;
  auto plenum = MakePlenumSolver(po.plenum_n_cells, n_level, mechanism);

  constexpr double tube_len_over_plenum_len = 1.47 / 0.56;
  int tube_n_cells = static_cast<int>(tube_len_over_plenum_len *
                                      static_cast<double>(po.plenum_n_cells));
  tube_n_cells = tube_n_cells - tube_n_cells % 8;
  std::vector<fub::amrex::IntegratorContext> tubes;

  //  ::amrex::RealBox inlet{{-0.1, -0.015, -0.015}, {0.05, +0.015, +0.015}};

  std::vector<fub::amrex::BlockConnection> connectivity{};

  auto MakeConnection = [&](int k) {
    TubeSolverOptions opts{};
    opts.max_refinement_level = po.max_refinement_level;
    opts.n_cells = tube_n_cells;
    opts.phi = k * alpha;
    tubes.push_back(MakeTubeSolver(mechanism, opts, vm));
    fub::amrex::BlockConnection connection;
    connection.direction = fub::Direction::X;
    connection.side = 0;
    connection.plenum.id = 0;
    connection.tube.id = k;
    connection.tube.mirror_box = tubes[k]
                                     .GetGriddingAlgorithm()
                                     ->GetPatchHierarchy()
                                     .GetGeometry(0)
                                     .Domain();
    connection.plenum.mirror_box =
        BoxWhichContains(DomainAroundCenter(Center(-0.03, k * alpha), 0.03),
                         plenum.GetGeometry(0));
    return connection;
  };

  connectivity.push_back(MakeConnection(0));
  connectivity.push_back(MakeConnection(1));
  connectivity.push_back(MakeConnection(2));
  connectivity.push_back(MakeConnection(3));
  connectivity.push_back(MakeConnection(4));

  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};
  fub::IdealGasMix<Plenum_Rank> equation{mechanism};

  fub::amrex::MultiBlockIntegratorContext context(
      fub::FlameMasterReactor(mechanism), std::move(tubes), {std::move(plenum)},
      std::move(connectivity));

  fub::DimensionalSplitLevelIntegrator system_solver(fub::int_c<Plenum_Rank>,
                                                     context);

  fub::amrex::MultiBlockIgniteDetonation ignition{
      tube_equation, context.GetGriddingAlgorithm(),
      fub::amrex::IgniteDetonationOptions(vm)};

  fub::DimensionalSplitSystemSourceSolver ign_solver(system_solver, ignition);

  fub::amrex::MultiBlockKineticSouceTerm source_term{
      fub::IdealGasMix<Tube_Rank>{mechanism}, context.GetGriddingAlgorithm()};

  fub::DimensionalSplitSystemSourceSolver solver{ign_solver, source_term};

  std::string base_name = "MultiTube";

  auto output =
      [&](const std::shared_ptr<fub::amrex::MultiBlockGriddingAlgorithm>&
              gridding,
          std::ptrdiff_t cycle, fub::Duration, int = 0) {
        auto tubes = gridding->GetTubes();
        int k = 0;
        for (auto& tube : tubes) {
          std::string name =
              fmt::format("{}/Tube_{}/plt{:05}", base_name, k, cycle);
          ::amrex::Print() << "Start output to '" << name << "'.\n";
          fub::amrex::WritePlotFile(name, tube->GetPatchHierarchy(),
                                    tube_equation);
          ::amrex::Print() << "Finished output to '" << name << "'.\n";
          k = k + 1;
        }
        std::string name = fmt::format("{}/Plenum/plt{:05}", base_name, cycle);
        ::amrex::Print() << "Start output to '" << name << "'.\n";
        fub::amrex::cutcell::WritePlotFile(
            name, gridding->GetPlena()[0]->GetPatchHierarchy(), equation);
        ::amrex::Print() << "Finished output to '" << name << "'.\n";
      };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint());
  fub::RunOptions run_options = fub::GetRunOptions(vm);
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  {
    fub::amrex::ScopeGuard _{};
    auto vm = ParseCommandLine(argc, argv);
    if (vm) {
      MyMain(*vm);
    }
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}
