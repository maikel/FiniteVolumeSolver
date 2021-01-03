// Copyright (c) 2020 Maikel Nadolski
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

#include "fub/AMReX/boundary_condition/GenericPressureValveBoundary.hpp"
#include "fub/equations/PerfectGasMix.hpp"
#include "fub/equations/perfect_gas_mix/IgnitionDelayKinetics.hpp"
#include "fub/flux_method/MusclHancockMethod2.hpp"

#include "fub/AMReX/multi_block/MultiBlockBoundary2.hpp"
#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm2.hpp"
#include "fub/AMReX/multi_block/MultiBlockIntegratorContext2.hpp"

#include "fub/AMReX/cutcell/boundary_condition/IsentropicPressureExpansion.hpp"
#include "fub/AMReX/cutcell/boundary_condition/MachnumberBoundary.hpp"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <boost/filesystem.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Translation.H>

#include <cmath>
#include <iostream>

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 2;

static_assert(AMREX_SPACEDIM == 2);

static constexpr double r_tube = 0.015;

struct InitialDataInTube {
  using Equation = fub::PerfectGasMix<1>;
  using Complete = fub::Complete<Equation>;
  using KineticState = fub::KineticState<Equation>;
  using Conservative = fub::Conservative<Equation>;

  Equation equation_;
  double x_0_;

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration /*time*/) const {
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& data = patch_level.data;
    fub::amrex::ForEachFab(
        fub::execution::openmp, data, [&](amrex::MFIter& mfi) {
          KineticState state(equation_);
          Complete complete(equation_);
          fub::Array<double, 1, 1> velocity{0.0};
          fub::View<Complete> states = fub::amrex::MakeView<Complete>(
              data[mfi], equation_, mfi.tilebox());
          fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
            const double x = geom.CellCenter(int(i), 0);
            const double rel_x = x - x_0_;
            const double pressure = 1.0;
            state.temperature = (rel_x < 0.5) ? 1.0 : 2.5;
            state.density = pressure / state.temperature / equation_.Rspec;
            state.mole_fractions[0] = 0.0;
            state.mole_fractions[1] = (rel_x < 0.4) ? 1.0 : 0.0;
            state.mole_fractions[2] = !(rel_x < 0.4) ? 1.0 : 0.0;
            fub::euler::CompleteFromKineticState(equation_, complete, state,
                                                 velocity);
            fub::Store(states, complete, {i});
          });
        });
  }
};

auto MakeTubeSolver(const fub::ProgramOptions& options,
                    const std::shared_ptr<fub::CounterRegistry>& counters) {
  using namespace fub::amrex;

  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << "==================== Tube =========================";

  CartesianGridGeometry grid_geometry(fub::GetOptions(options, "GridGeometry"));
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  PatchHierarchyOptions hierarchy_options(
      fub::GetOptions(options, "PatchHierarchy"));
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  using Eq = fub::PerfectGasMix<Tube_Rank>;
  using Complete = Eq::Complete;
  fub::PerfectGasMix<Tube_Rank> equation{};
  equation.n_species = 2;
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 1e-2),
                            std::make_pair(&Complete::pressure, 1e-2)};

  DataDescription desc = MakeDataDescription(equation);

  ::amrex::Box refine_box{{grid_geometry.cell_dimensions[0] - 5, 0},
                          {grid_geometry.cell_dimensions[0] - 1, 0}};
  ConstantBox constant_box{refine_box};

  InitialDataInTube initial_data{equation, grid_geometry.coordinates.lo()[0]};

  fub::perfect_gas_mix::IgnitionDelayKinetics<1> source_term{equation};

  static constexpr double buffer = 0.5;
  static constexpr double pbufwidth = 1e-10;
  const double lambda =
      -std::log(source_term.options.Yign / source_term.options.Yinit /
                source_term.options.tau);
  auto fill_f = [lambda, Yign = source_term.options.Yign,
                 tau = source_term.options.tau](double t) {
    return Yign * std::exp(lambda * (tau - t));
  };
  auto fill_f_val = [fill_f](double t) {
    const double ff = fill_f(t);
    if (ff == 1.0) {
      return std::numeric_limits<double>::max();
    }
    return ff / (1.0 - ff);
  };
  static constexpr double t_ignite = 1.1753;
  static constexpr double t_ignite_diff = t_ignite - 1.0;
  auto inflow_function = [fill_f_val](
                             const fub::PerfectGasMix<1>&,
                             fub::KineticState<fub::PerfectGasMix<1>>& kin,
                             fub::Duration tp, const amrex::MultiFab&,
                             const fub::amrex::GriddingAlgorithm&, int) {
    const double dt = tp.count();
    kin.temperature = 1.0;
    kin.density = 1.0;
    if (dt > buffer) {
      kin.mole_fractions[0] = fill_f_val(dt - t_ignite_diff);
      kin.mole_fractions[1] = 1.0;
    } else {
      // FR
      kin.mole_fractions[0] = 0.0;
      kin.mole_fractions[1] = 0.0;
    }
    kin.mole_fractions[2] = std::max(
        0.0, 10.0 * std::min(1.0, 1.0 - (dt - buffer) / pbufwidth / buffer));
    const double sum = kin.mole_fractions.sum();
    kin.mole_fractions /= sum;
  };

  fub::amrex::GenericPressureValveBoundary valve(equation, inflow_function, {});

  BoundarySet boundaries{{valve}};

  // If a checkpoint path is specified we will fill the patch hierarchy with
  // data from the checkpoint file, otherwise we will initialize the data by
  // the initial data function.
  std::shared_ptr<GriddingAlgorithm> gridding = [&] {
    std::string checkpoint =
        fub::GetOptionOr(options, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      BOOST_LOG(log) << "Initialize grid by initial condition...";
      std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(desc, grid_geometry, hierarchy_options), initial_data,
          TagAllOf(gradient, constant_box), boundaries);
      gridding->GetPatchHierarchy().SetCounterRegistry(counters);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
      checkpoint = fmt::format("{}/Tube", checkpoint);
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
      PatchHierarchy h = ReadCheckpointFile(checkpoint, desc, grid_geometry,
                                            hierarchy_options);
      std::shared_ptr<GriddingAlgorithm> gridding =
          std::make_shared<GriddingAlgorithm>(std::move(h), initial_data,
                                              TagAllOf(gradient, constant_box),
                                              boundaries);
      return gridding;
    }
  }();

  using HLLEM = fub::perfect_gas::HllemMethod<Eq>;
  using CharacteristicsReconstruction = fub::FluxMethod<fub::MusclHancock2<
      Eq,
      fub::CharacteristicsGradient<
          Eq, fub::CentralDifferenceGradient<fub::MinModLimiter>>,
      fub::CharacteristicsReconstruction<Eq>, HLLEM>>;

  CharacteristicsReconstruction flux_method{equation};
  HyperbolicMethod method{FluxMethodAdapter(flux_method),
                          EulerForwardTimeIntegrator(),
                          Reconstruction(equation)};

  const int scratch_gcw = 8;
  const int flux_gcw = 6;

  IntegratorContext context(gridding, method, scratch_gcw, flux_gcw);

  BOOST_LOG(log) << "==================== End Tube =========================";

  return context;
}

auto MakePlenumSolver(const std::map<std::string, pybind11::object>& options) {
  fub::SeverityLogger log = fub::GetInfoLogger();
  BOOST_LOG(log) << "==================== Plenum =========================";
  using namespace fub::amrex::cutcell;
  auto MakePolygon = [](auto... points) {
    fub::Polygon::Vector xs{};
    fub::Polygon::Vector ys{};

    (xs.push_back(std::get<0>(points)), ...);
    (ys.push_back(std::get<1>(points)), ...);

    return fub::Polygon(std::move(xs), std::move(ys));
  };

  double r_inlet_start = r_tube;
  double r_inlet_end = 0.5 * r_tube;
  {
    const fub::ProgramOptions inlet_options =
        fub::GetOptions(options, "InletGeometry");
    r_inlet_start = fub::GetOptionOr(inlet_options, "r_start", r_inlet_start);
    r_inlet_end = fub::GetOptionOr(inlet_options, "r_end", r_inlet_end);
  }

  auto ConvergentInlet = [&](double height,
                             const std::array<double, Plenum_Rank>& center) {
    const double xlo = center[0] - height;
    const double xhi = center[0];
    const double r = r_inlet_start;
    const double r2 = r_inlet_end;
    const double xdiv = xhi - 4.0 * r;
    auto polygon =
        MakePolygon(std::pair{xlo, r}, std::pair{xdiv, r}, std::pair{xhi, r2},
                    std::pair{xhi, -r2}, std::pair{xdiv, -r},
                    std::pair{xlo, -r}, std::pair{xlo, r});
    auto tube_in_zero = fub::amrex::Geometry(fub::Invert(polygon));
    amrex::RealArray real_center{center[0], center[1]};
    auto tube_in_center = amrex::EB2::TranslationIF(tube_in_zero, real_center);
    return tube_in_center;
  };

  const fub::ProgramOptions blade_options =
      fub::GetOptions(options, "BladeGeometry");
  const double x0 = fub::GetOptionOr(blade_options, "x0", 4.0 * r_tube);
  const double dx = fub::GetOptionOr(blade_options, "dx", 4.0 * r_tube);
  const double y0 = fub::GetOptionOr(blade_options, "y0", 0.0);
  const double dy = fub::GetOptionOr(blade_options, "dy", r_tube);
  const double dBox = fub::GetOptionOr(blade_options, "dBox", r_tube);
  auto Box = [&](double y) {
    amrex::RealArray lo{x0, y0 + y};
    amrex::RealArray hi{x0 + dx, y0 + y + dy};
    return amrex::EB2::BoxIF(lo, hi, false);
  };

  auto embedded_boundary = amrex::EB2::makeUnion(
      amrex::EB2::makeIntersection(
          amrex::EB2::PlaneIF({0.0, 0.0}, {1.0, 0.0}, false),
          ConvergentInlet(2.0, {0.0, 0.0})),
      Box(0.0), Box(dy + dBox), Box(2 * (dy + dBox)), Box(3 * (dy + dBox)),
      Box(4 * (dy + dBox)), Box(5 * (dy + dBox)), Box(6 * (dy + dBox)),
      Box(-1.0 * (dy + dBox)), Box(-2.0 * (dy + dBox)), Box(-3.0 * (dy + dBox)),
      Box(-4.0 * (dy + dBox)), Box(-5.0 * (dy + dBox)),
      Box(-6.0 * (dy + dBox)));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::amrex::CartesianGridGeometry grid_geometry =
      fub::GetOptions(options, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  PatchHierarchyOptions hierarchy_options =
      fub::GetOptions(options, "PatchHierarchy");
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  BOOST_LOG(log) << "Compute EB level set data...";
  std::shared_ptr<fub::CounterRegistry> registry =
      std::make_shared<fub::CounterRegistry>();
  {
    fub::Timer timer = registry->get_timer("MakeIndexSpaces");
    hierarchy_options.index_spaces =
        MakeIndexSpaces(shop, grid_geometry, hierarchy_options);
  }

  fub::PerfectGasMix<Plenum_Rank> equation{};
  equation.n_species = 2;

  fub::Complete<fub::PerfectGasMix<Plenum_Rank>> left(equation);
  fub::Complete<fub::PerfectGasMix<Plenum_Rank>> right(equation);
  {
    using namespace std::literals;
    const fub::ProgramOptions initial_options =
        fub::GetOptions(options, "InitialCondition");
    const fub::ProgramOptions left_options =
        fub::GetOptions(initial_options, "left");
    std::array<double, 3> moles{0.0, 0.0, 1.0};
    moles = fub::GetOptionOr(left_options, "moles", moles);
    double temperature = fub::GetOptionOr(left_options, "temperature", 1.0);
    double density = fub::GetOptionOr(left_options, "density", 1.0);
    std::array<double, Plenum_Rank> velocity{};
    velocity = fub::GetOptionOr(left_options, "velocity", velocity);
    fub::KineticState<fub::PerfectGasMix<Plenum_Rank>> kin(equation);
    kin.temperature = temperature;
    kin.density = density;
    kin.mole_fractions = fub::AsEigenVector(moles);
    fub::euler::CompleteFromKineticState(equation, left, kin,
                                         fub::AsEigenVector(velocity));

    const fub::ProgramOptions right_options =
        fub::GetOptions(initial_options, "right");
    moles = std::array<double, 3>{0.0, 0.0, 1.0};
    moles = fub::GetOptionOr(right_options, "moles", moles);
    temperature = fub::GetOptionOr(right_options, "temperature", 1.0);
    density = fub::GetOptionOr(right_options, "pressure", 1.0);
    velocity = std::array<double, Plenum_Rank>{};
    velocity = fub::GetOptionOr(right_options, "velocity", velocity);
    kin.temperature = temperature;
    kin.density = density;
    kin.mole_fractions = fub::AsEigenVector(moles);
    fub::euler::CompleteFromKineticState(equation, right, kin,
                                         fub::AsEigenVector(velocity));
  }
  RiemannProblem initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.0),
                              left, right);

  //  using Complete = fub::Complete<fub::PerfectGas<Plenum_Rank>>;
  using Complete = fub::Complete<fub::PerfectGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&Complete::pressure, 0.05},
                             std::pair{&Complete::density, 0.05}};

  const int scratch_gcw = 8;
  const int flux_gcw = 6;

  fub::amrex::cutcell::MachnumberBoundaryOptions pb{};
  pb.dir = fub::Direction::X;
  pb.side = 1;
  const double y_len = grid_geometry.coordinates.length(1);
  const double rel_y_width = 0.5 * r_inlet_start / y_len;
  ::amrex::IntVect pblo{grid_geometry.cell_dimensions[0], 0};
  ::amrex::IntVect pbhi{grid_geometry.cell_dimensions[0] + 7, int(rel_y_width * grid_geometry.cell_dimensions[1])};
  pb.boundary_section = ::amrex::Box(pblo, pbhi);

  fub::amrex::cutcell::MachnumberBoundary<fub::PerfectGasMix<2>> pbdry(equation, pb);
  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::X, 1},
                                  ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::Y, 0},
                                  ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::Y, 1},
                                  pbdry}};

  ::amrex::RealBox xbox = grid_geometry.coordinates;
  ::amrex::Geometry coarse_geom = fub::amrex::GetCoarseGeometry(grid_geometry);

  ::amrex::RealBox inlet{{xbox.lo(0), -0.15}, {0.18, +0.15}};
  ::amrex::Box refine_box = fub::amrex::BoxWhichContains(inlet, coarse_geom);
  ConstantBox constant_refinebox{refine_box};

  std::shared_ptr gridding = [&] {
    std::string checkpoint =
        fub::GetOptionOr(options, "checkpoint", std::string{});
    if (checkpoint.empty()) {
      BOOST_LOG(log) << "Initialize grid by initial condition...";
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          PatchHierarchy(equation, grid_geometry, hierarchy_options),
          initial_data, TagAllOf(TagCutCells(), constant_refinebox),
          boundary_condition);
      grid->InitializeHierarchy(0.0);
      return grid;
    } else {
      checkpoint += "/Plenum";
      BOOST_LOG(log) << "Initialize grid from given checkpoint '" << checkpoint
                     << "'.";
      PatchHierarchy h = ReadCheckpointFile(
          checkpoint, fub::amrex::MakeDataDescription(equation), grid_geometry,
          hierarchy_options);
      std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
          std::move(h), initial_data,
          TagAllOf(TagCutCells(), constant_refinebox), boundary_condition);
      return grid;
    }
  }();

  using HLLEM = fub::perfect_gas::HllemMethod<fub::PerfectGasMix<2>>;

  using CharacteristicsReconstruction = fub::FluxMethod<fub::MusclHancock2<
      fub::PerfectGasMix<Plenum_Rank>,
      fub::CharacteristicsGradient<
          fub::PerfectGasMix<Plenum_Rank>,
          fub::CentralDifferenceGradient<fub::MinModLimiter>>,
      fub::CharacteristicsReconstruction<fub::PerfectGasMix<Plenum_Rank>>,
      HLLEM>>;

  HLLEM hllem_method{equation};
  CharacteristicsReconstruction flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hllem_method);

  HyperbolicMethod method{FluxMethod{cutcell_method}, TimeIntegrator{},
                          Reconstruction{equation}};

  IntegratorContext context(gridding, method, scratch_gcw, flux_gcw);

  BOOST_LOG(log) << "==================== End Plenum =========================";

  return context;
}

void MyMain(const std::map<std::string, pybind11::object>& vm);

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

void WriteCheckpoint(const std::string& path,
                     const fub::amrex::MultiBlockGriddingAlgorithm2& grid) {
  auto tubes = grid.GetTubes();
  std::string name = fmt::format("{}/Tube", path);
  fub::amrex::WriteCheckpointFile(name, tubes[0]->GetPatchHierarchy());
  name = fmt::format("{}/Plenum", path);
  fub::amrex::cutcell::WriteCheckpointFile(
      name, grid.GetPlena()[0]->GetPatchHierarchy());
}

void MyMain(const std::map<std::string, pybind11::object>& vm) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};

  std::vector<fub::amrex::cutcell::IntegratorContext> plenum{};
  std::vector<fub::amrex::IntegratorContext> tubes{};
  std::vector<fub::amrex::BlockConnection> connectivity{};

  plenum.push_back(MakePlenumSolver(fub::GetOptions(vm, "Plenum")));
  auto counter_database = plenum[0].GetCounterRegistry();

  auto&& tube = MakeTubeSolver(fub::GetOptions(vm, "Tube"), counter_database);
  fub::amrex::BlockConnection connection;
  connection.direction = fub::Direction::X;
  connection.side = 0;
  connection.ghost_cell_width = 8;
  connection.plenum.id = 0;
  connection.tube.id = 0;
  connection.tube.mirror_box =
      tube.GetGriddingAlgorithm()->GetPatchHierarchy().GetGeometry(0).Domain();
  connection.plenum.mirror_box = plenum[0]
                                     .GetGriddingAlgorithm()
                                     ->GetPatchHierarchy()
                                     .GetGeometry(0)
                                     .Domain();
  tubes.push_back(std::move(tube));
  connectivity.push_back(connection);

  fub::PerfectGasMix<Plenum_Rank> plenum_equation{};
  plenum_equation.n_species = 2;
  fub::PerfectGasMix<Tube_Rank> tube_equation{};
  tube_equation.n_species = 2;

  fub::amrex::MultiBlockIntegratorContext2 context(
      tube_equation, plenum_equation, std::move(tubes), std::move(plenum),
      std::move(connectivity));

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<Plenum_Rank>, std::move(context), fub::GodunovSplitting{});

  fub::amrex::MultiBlockSourceTerm<
      fub::perfect_gas_mix::IgnitionDelayKinetics<1>>
      source_term(
          {fub::perfect_gas_mix::IgnitionDelayKinetics<1>(tube_equation)});

  fub::SplitSystemSourceLevelIntegrator level_integrator(
      std::move(system_solver), std::move(source_term),
      fub::GodunovSplitting{});

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  using namespace fub::amrex;
  struct MakeCheckpoint
      : public fub::OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm2> {
    std::string directory_ = "ConvergentNozzle";
    MakeCheckpoint(const fub::ProgramOptions& options)
        : OutputAtFrequencyOrInterval(options) {
      directory_ = fub::GetOptionOr(options, "directory", directory_);
    }
    void operator()(const MultiBlockGriddingAlgorithm2& grid) override {
      std::string name = fmt::format("{}/{:09}", directory_, grid.GetCycles());
      fub::SeverityLogger log = fub::GetInfoLogger();
      BOOST_LOG(log) << fmt::format("Write Checkpoint to '{}'!", name);
      WriteCheckpoint(name, grid);
    }
  };

  fub::OutputFactory<MultiBlockGriddingAlgorithm2> factory{};
  factory.RegisterOutput<MakeCheckpoint>("Checkpoint");
  using CounterOutput =
      fub::CounterOutput<fub::amrex::MultiBlockGriddingAlgorithm2,
                         std::chrono::milliseconds>;
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<
      MultiBlockPlotfileOutput2<fub::PerfectGasMix<1>, fub::PerfectGasMix<2>>>(
      "Plotfiles", tube_equation, plenum_equation);
  fub::MultipleOutputs<MultiBlockGriddingAlgorithm2> outputs(
      std::move(factory), fub::GetOptions(vm, "Output"));

  outputs(*solver.GetGriddingAlgorithm());
  fub::SeverityLogger log = fub::GetInfoLogger();
  fub::RunOptions run_options = fub::GetOptions(vm, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);
  fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
}
