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
#include <boost/program_options.hpp>

#include <cmath>
#include <iostream>

#include <xmmintrin.h>

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 3;

struct TemperatureRamp {
  using Complete = fub::IdealGasMix<Tube_Rank>::Complete;

  fub::IdealGasMix<Tube_Rank> equation_;
  double ignition_pos_{-1.0};
  double air_position_{0.0};
  double equiv_raito_{1.0};

  void InitializeData(::amrex::MultiFab& data, const ::amrex::Geometry& geom) {
    fub::FlameMasterReactor& reactor = equation_.GetReactor();
    const double high_temp = 1450.0;
    const double low_temp = 300.0;
    Complete complete(equation_);
    Complete burnt{equation_};
    Complete cold(equation_);
    Complete air(equation_);

    double h2_moles = 42.0 * equiv_raito_;
    double o2_moles = std::max(21.0 - 0.5 * h2_moles, 0.0);

    std::string fuel_moles = fmt::format("N2:79,O2:21,H2:{}", h2_moles);
    reactor.SetMoleFractions(fuel_moles);
    reactor.SetTemperature(low_temp);
    reactor.SetPressure(101325.0);
    equation_.CompleteFromReactor(cold);

    reactor.SetMoleFractions("N2:79,O2:21");
    reactor.SetTemperature(low_temp);
    reactor.SetPressure(101325.0);
    equation_.CompleteFromReactor(air);

    reactor.SetMoleFractions(
        fmt::format("N2:79,O2:{},H2O:{}", o2_moles, h2_moles));
    reactor.SetTemperature(high_temp);
    reactor.SetPressure(101325.0);
    equation_.CompleteFromReactor(burnt);

    fub::amrex::ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
        double x[AMREX_SPACEDIM] = {};
        geom.CellCenter(::amrex::IntVect{int(i)}, x);
        if (x[0] < ignition_pos_ - 0.05) {
          fub::Store(states, cold, {i});
        } else if (x[0] < ignition_pos_) {
          const double x0 = ignition_pos_ - x[0];
          const double d = std::clamp(x0 / 0.05, 0.0, 1.0);
          reactor.SetMoleFractions(fuel_moles);
          reactor.SetTemperature(d * low_temp + (1.0 - d) * high_temp);
          reactor.SetPressure(101325.0);
          equation_.CompleteFromReactor(complete);
          fub::Store(states, complete, {i});
        } else if (x[0] < ignition_pos_ + 0.05) {
          const double x0 = x[0] - ignition_pos_;
          const double d = std::clamp(x0 / 0.05, 0.0, 1.0);
          reactor.SetMoleFractions(fuel_moles);
          reactor.SetTemperature(d * low_temp + (1.0 - d) * high_temp);
          reactor.SetPressure(101325.0);
          equation_.CompleteFromReactor(complete);
          fub::Store(states, complete, {i});
        } else if (x[0] < air_position_) {
          fub::Store(states, cold, {i});
        } else {
          fub::Store(states, air, {i});
        }
      });
    });
  }
};

struct ProgramOptions {
  int max_cycles{-1};
  double plenum_domain_length{1.0};
  int plenum_n_cells{128};
  int max_refinement_level{1};
  double tube_domain_length{1.47};
  double tube_ignition_position{0.6 - 1.47};
  double tube_air_position{0.0};
  double tube_equiv_ratio{1.0};
  std::string plenum_checkpoint{};
  std::string tube_checkpoint{};
};

auto MakeTubeSolver(const ProgramOptions& po, fub::Burke2012& mechanism) {
  const int num_cells = [&] {
    const double len_ratio = po.tube_domain_length / po.plenum_domain_length;
    int cells =
        static_cast<int>(len_ratio * static_cast<double>(po.plenum_n_cells));
    cells -= cells % 8;
    return cells;
  }();

  const double length = po.tube_domain_length;
  const double tube_radius = 0.015;
  const int n_level = po.max_refinement_level;

  const std::array<int, AMREX_SPACEDIM> n_cells{num_cells, 1, 1};
  const std::array<double, AMREX_SPACEDIM> xlower{-length - 0.03, -tube_radius,
                                                  -tube_radius};
  const std::array<double, AMREX_SPACEDIM> xupper{-0.03, tube_radius,
                                                  tube_radius};
  const std::array<int, AMREX_SPACEDIM> periodicity{0, 0, 0};

  amrex::RealBox xbox(xlower, xupper);
  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};

  using namespace fub::amrex;

  CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = n_level;
  hier_opts.refine_ratio = amrex::IntVect{2, 1, 1};

  amrex::Geometry geom(
      amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1}}, &xbox,
      -1, periodicity.data());
  geom.refine(hier_opts.refine_ratio);

  using Complete = fub::IdealGasMix<1>::Complete;
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 5e-3),
                            std::make_pair(&Complete::pressure, 5e-3)};

  ::amrex::Box refine_box{{num_cells - 5, 0, 0}, {num_cells - 1, 0, 0}};
  ConstantBox constant_box{refine_box};

  TemperatureRamp initial_data{equation};
  initial_data.ignition_pos_ = po.tube_ignition_position;
  initial_data.air_position_ = po.tube_air_position;
  initial_data.equiv_raito_ = po.tube_equiv_ratio;

  BoundarySet boundaries{
      {ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::X, 0},
       TransmissiveBoundary{fub::Direction::X, 1}}};

  PatchHierarchy hierarchy = [&] {
    if (po.tube_checkpoint.empty()) {
      return PatchHierarchy(equation, geometry, hier_opts);
    }
    return ReadCheckpointFile(po.tube_checkpoint,
                              fub::amrex::MakeDataDescription(equation),
                              geometry, hier_opts);
  }();
  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      hierarchy, initial_data, TagAllOf(gradient, constant_box, TagBuffer(2)),
      boundaries);
  gridding->InitializeHierarchy(0.0);
  if (po.tube_checkpoint.empty()) {
    gridding->InitializeHierarchy(0.0);
  }

  //  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Tube_Rank>> signals{};
  //  fub::HllMethod hll_method(equation, signals);
  // fub::MusclHancockMethod flux_method{equation, hll_method};
  fub::ideal_gas::MusclHancockPrimMethod<Tube_Rank> flux_method(equation);

  HyperbolicMethod method{FluxMethod(fub::execution::seq, flux_method),
                          ForwardIntegrator(fub::execution::seq),
                          Reconstruction(fub::execution::seq, equation)};

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

auto MakePlenumSolver(const ProgramOptions& po, fub::Burke2012& mechanism) {
  const int num_cells = po.plenum_n_cells;
  const std::array<int, Plenum_Rank> n_cells{num_cells, num_cells, num_cells};
  const double length = po.plenum_domain_length;
  const double half = 0.5 * po.plenum_domain_length;
  const int n_level = po.max_refinement_level;

  const std::array<double, Plenum_Rank> xlower{-0.03, -half, -half};
  const std::array<double, Plenum_Rank> xupper{length - 0.03, +0.15, +0.15};
  const std::array<int, Plenum_Rank> periodicity{0, 0, 0};

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(
      amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1}}, &xbox,
      -1, periodicity.data());

  auto embedded_boundary = amrex::EB2::makeIntersection(
      amrex::EB2::PlaneIF({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, false),
      amrex::EB2::CylinderIF(0.015, -1.0, 0, {1e6, 0.0, 0.0}, true));
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
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.01},
                             std::pair{&State::density, 0.01}};

  ::amrex::RealBox inlet{{-0.1, -0.015, -0.015}, {0.01, +0.015, +0.015}};
  const ::amrex::Box refine_box = BoxWhichContains(inlet, coarse_geom);
  ConstantBox constant_box{refine_box};

  //  const double p0 = 101325.0;
  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1},
                                  TransmissiveBoundary{fub::Direction::Z, 0},
                                  TransmissiveBoundary{fub::Direction::Z, 1}}};

  PatchHierarchy hierarchy = [&] {
    if (po.plenum_checkpoint.empty()) {
      return PatchHierarchy(equation, geometry, options);
    }
    return ReadCheckpointFile(po.plenum_checkpoint,
                              fub::amrex::MakeDataDescription(equation),
                              geometry, options);
  }();
  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      std::move(hierarchy), initial_data,
      TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
      boundary_condition);
  if (po.plenum_checkpoint.empty()) {
    gridding->InitializeHierarchy(0.0);
  }

  // Make Solver

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  //  fub::MusclHancockMethod flux_method(equation, hll_method);
  fub::ideal_gas::MusclHancockPrimMethod<Plenum_Rank> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{
      FluxMethod{fub::execution::openmp_simd, cutcell_method},
      fub::amrex::cutcell::TimeIntegrator{},
      Reconstruction{fub::execution::openmp_simd, equation}};

  return fub::amrex::cutcell::IntegratorContext(gridding, method);
}

void MyMain(const ProgramOptions& po,
            const boost::program_options::variables_map& vm);

std::optional<std::pair<ProgramOptions, boost::program_options::variables_map>>
ParseCommandLine(int argc, char** argv) {
  namespace po = boost::program_options;
  ProgramOptions opts;
  po::options_description desc = fub::GetCommandLineRunOptions();
  desc.add_options()(
      "plenum_checkpoint",
      po::value<std::string>(&opts.plenum_checkpoint)->default_value(""),
      "Restart the simulation from a given plenum checkpoint.")(
      "tube_checkpoint",
      po::value<std::string>(&opts.tube_checkpoint)->default_value(""),
      "Restart the simulation from a given tube checkpoint.")(
      "n_cells",
      po::value<int>(&opts.plenum_n_cells)->default_value(opts.plenum_n_cells),
      "Set number of cells in each direction for the plenum")(
      "max_refinement_level",
      po::value<int>(&opts.max_refinement_level)
          ->default_value(opts.max_refinement_level),
      "Set the maximal number of refinement levels (>= 1)")(
      "plenum_len",
      po::value<double>(&opts.plenum_domain_length)
          ->default_value(opts.plenum_domain_length),
      "Set the base length for the plenum")(
      "tube_len",
      po::value<double>(&opts.tube_domain_length)
          ->default_value(opts.tube_domain_length),
      "Set the base length for the tube")(
      "tube_ignition",
      po::value<double>(&opts.tube_ignition_position)
          ->default_value(opts.tube_ignition_position),
      "Set the position for the ignition of the detonation inside the tube")(
      "tube_air_position",
      po::value<double>(&opts.tube_air_position)
          ->default_value(opts.tube_air_position),
      "Sets the starting position for an air buffer in the tube.")(
      "tube_equiv_ratio",
      po::value<double>(&opts.tube_equiv_ratio)
          ->default_value(opts.tube_equiv_ratio),
      "Sets the equivalence ratio of the fuel in the tube");
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
  amrex::Print() << fmt::format("[Info] Simulation will run with the following "
                                "additional options:\n[Info]\n");
  amrex::Print() << fmt::format("[Info] max_refinement_level = {}\n",
                                opts.max_refinement_level);
  std::array<double, 3> plo{-0.03, -0.5 * opts.plenum_domain_length,
                            -0.5 * opts.plenum_domain_length};
  std::array<double, 3> pup{opts.plenum_domain_length - 0.03,
                            +0.5 * opts.plenum_domain_length,
                            +0.5 * opts.plenum_domain_length};
  amrex::Print() << fmt::format(
      "[Info] plenum_domain = [{}, {}] x [{}, {}] x [{}, {}]\n", plo[0], pup[0],
      plo[1], pup[1], plo[2], pup[2]);
  amrex::Print() << fmt::format(
      "[Info] plenum_n_cells = {}\n[Info] plenum_dx = {}\n",
      opts.plenum_n_cells, opts.plenum_domain_length / opts.plenum_n_cells);
  std::array<double, 3> tlo{-opts.tube_domain_length - 0.03, -0.015, -0.015};
  std::array<double, 3> tup{-0.03, +0.015, +0.015};
  amrex::Print() << fmt::format(
      "[Info] tube_domain = [{}, {}] x [{}, {}] x [{}, {}]\n", tlo[0], tup[0],
      tlo[1], tup[1], tlo[2], tup[2]);
  const double ratio = opts.tube_domain_length / opts.plenum_domain_length;
  const int tube_n_cells =
      static_cast<int>(ratio * static_cast<double>(opts.plenum_n_cells));
  amrex::Print() << fmt::format(
      "[Info] tube_n_cells = {}\n[Info] tube_dx = {}\n", tube_n_cells,
      opts.tube_domain_length / tube_n_cells);
  amrex::Print() << fmt::format("[Info] tube_ignition_pos = {}\n",
                                opts.tube_ignition_position);
  amrex::Print() << fmt::format("[Info] tube_equiv_ratio = {}\n",
                                opts.tube_equiv_ratio);
  amrex::Print() << fmt::format("[Info] tube_air_position = {}\n",
                                opts.tube_air_position);

  if (!opts.plenum_checkpoint.empty() || !opts.tube_checkpoint.empty()) {
    if (opts.plenum_checkpoint.empty() || opts.tube_checkpoint.empty()) {
      amrex::Print() << "[Error] Only one Checkpoint file specified but you "
                        "need a checkpoint for each domain.\n";
      return {};
    }
    amrex::Print() << "[Info]\n[Info] Restart from a checkpoint!\n";
    amrex::Print() << "[Info] Plenum Checkpoint: " << opts.plenum_checkpoint
                   << '\n';
    amrex::Print() << "[Info] Tube Checkpoint: " << opts.tube_checkpoint
                   << '\n';
  }

  return std::make_pair(opts, vm);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  {
    fub::amrex::ScopeGuard _{};
    auto po = ParseCommandLine(argc, argv);
    if (po) {
      MyMain(po->first, po->second);
    }
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}

void MyMain(const ProgramOptions& po,
            const boost::program_options::variables_map& vm) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::Burke2012 mechanism{};

  auto plenum = MakePlenumSolver(po, mechanism);
  auto tube = MakeTubeSolver(po, mechanism);

  ::amrex::RealBox inlet{{-0.1, -0.015, -0.015}, {0.05, +0.015, +0.015}};

  fub::amrex::BlockConnection connection;
  connection.direction = fub::Direction::X;
  connection.side = 0;
  connection.plenum.id = 0;
  connection.plenum.mirror_box = BoxWhichContains(inlet, plenum.GetGeometry(0));
  connection.tube.id = 0;
  connection.tube.mirror_box =
      tube.GetGriddingAlgorithm()->GetPatchHierarchy().GetGeometry(0).Domain();

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};
  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};

  fub::amrex::MultiBlockIntegratorContext context(
      fub::FlameMasterReactor(mechanism), {std::move(tube)},
      {std::move(plenum)}, {connection});

  fub::DimensionalSplitLevelIntegrator system_solver(fub::int_c<Plenum_Rank>,
                                                     std::move(context));

  fub::amrex::MultiBlockKineticSouceTerm source_term{
      tube_equation, system_solver.GetGriddingAlgorithm()};

  fub::DimensionalSplitSystemSourceSolver solver{system_solver, source_term};

  std::string base_name = "LongLinearShock_3d";

  std::vector<double> probes_buffer(6 * 3);
  fub::basic_mdspan<double, fub::extents<3, fub::dynamic_extent>> probes(
      probes_buffer.data(), 6);
  probes(0, 0) = 1.0E-6 * 0.03;
  probes(0, 1) = 4.0 * 0.03;
  probes(0, 2) = 9.0 * 0.03;
  probes(0, 3) = 14.0 * 0.03;
  probes(0, 4) = 19.0 * 0.03;
  probes(0, 5) = 24.0 * 0.03;

  std::vector<double> tube_probes_buffer(5 * 3);
  fub::basic_mdspan<double, fub::extents<3, fub::dynamic_extent>> tube_probes(
      tube_probes_buffer.data(), 5);
  tube_probes(0, 0) = -4.0 * 0.03;
  tube_probes(0, 1) = -9.0 * 0.03;
  tube_probes(0, 2) = -14.0 * 0.03;
  tube_probes(0, 3) = -19.0 * 0.03;
  tube_probes(0, 4) = -24.0 * 0.03;

  ::amrex::Box slice_box = [&](double z0) {
    const auto& plenum =
        context.GetGriddingAlgorithm()->GetPlena()[0]->GetPatchHierarchy();
    const int finest_level = plenum.GetNumberOfLevels() - 1;
    const ::amrex::Geometry& geom = plenum.GetGeometry(finest_level);
    const ::amrex::RealBox& probDomain = geom.ProbDomain();
    const double xlo[] = {probDomain.lo(0), z0, probDomain.lo(2)};
    const double* xhi = probDomain.hi();
    const ::amrex::RealBox slice_x(xlo, xhi);
    ::amrex::Box slice_box = BoxWhichContains(slice_x, geom);
    slice_box.setBig(2, slice_box.smallEnd(2));
    return slice_box;
  }(0.0);

  // Write Checkpoints 0min + every 5min
  //  const fub::Duration checkpoint_offest = std::chrono::minutes(5);
  //  fub::Duration next_checkpoint = std::chrono::minutes(0);
  int rank = -1;
  MPI_Comm_rank(context.GetMpiCommunicator(), &rank);
  auto output =
      [&](std::shared_ptr<fub::amrex::MultiBlockGriddingAlgorithm> gridding,
          auto cycle, auto timepoint, int output_num) {
        if (output_num == 0) {
          ::amrex::Print() << "Checkpointing.\n";
          fub::amrex::WriteCheckpointFile(
              fmt::format("{}/Checkpoint/Tube_{:05}", base_name, cycle),
              gridding->GetTubes()[0]->GetPatchHierarchy());
          fub::amrex::cutcell::WriteCheckpointFile(
              fmt::format("{}/Checkpoint/Plenum_{:05}", base_name, cycle),
              gridding->GetPlena()[0]->GetPatchHierarchy());
          fub::amrex::WritePlotFile(
              fmt::format("{}/Plotfile/Tube_plt{}", base_name, cycle),
              gridding->GetTubes()[0]->GetPatchHierarchy(), tube_equation);
          fub::amrex::cutcell::WritePlotFile(
              fmt::format("{}/Plotfile/Plenum_plt{}", base_name, cycle),
              gridding->GetPlena()[0]->GetPatchHierarchy(), equation);
          fub::amrex::cutcell::Write2Dfrom3D(
              fmt::format("{}/Plenum_{:05}.dat", base_name, cycle),
              gridding->GetPlena()[0]->GetPatchHierarchy(), slice_box, equation,
              timepoint, cycle, context.GetMpiCommunicator());
          fub::amrex::WriteTubeData(
              fmt::format("{}/Tube_{:05}.dat", base_name, cycle),
              gridding->GetTubes()[0]->GetPatchHierarchy(), tube_equation,
              timepoint, cycle, context.GetMpiCommunicator());
        }
        if (output_num >= 0) {
          ::amrex::Print() << "Start Output for Probes.\n";
          {
            std::vector<double> buffer =
                GatherStates(gridding->GetTubes()[0]->GetPatchHierarchy(),
                             tube_probes, context.GetMpiCommunicator());
            fub::mdspan<const double, 2> states(
                buffer.data(), tube_probes.extent(1),
                buffer.size() / tube_probes.extent(1));
            ::amrex::Print()
                << fmt::format("{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}\n",
                               "Position", "Time", "Density", "VelocityX",
                               "SpeedOfSound", "Temperature", "Pressure");
            for (int i = tube_probes.extent(1) - 1; i >= 0; --i) {
              const double rho = states(i, 0);
              const double u = states(i, 1) / rho;
              const double T = states(i, 16);
              const double p = states(i, 14);
              const double a = states(i, 15);
              const double t = timepoint.count();
              const double x = tube_probes(0, i);
              ::amrex::Print()
                  << fmt::format("{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< "
                                 "24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}\n",
                                 x, t, rho, u, a, T, p);
            }
          }

          {
            std::vector<double> buffer =
                GatherStates(gridding->GetPlena()[0]->GetPatchHierarchy(),
                             probes, context.GetMpiCommunicator());
            fub::mdspan<const double, 2> states(buffer.data(), probes.extent(1),
                                                buffer.size() /
                                                    probes.extent(1));
            ::amrex::Print() << fmt::format(
                "{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}"
                "{:24s}{:24s}\n",
                "Position", "Time", "Density", "VelocityX", "VelocityY",
                "VelocityZ", "SpeedOfSound", "Temperature", "Pressure");
            for (int i = 0; i < probes.extent(1); ++i) {
              const double rho = states(i, 0);
              const double u = states(i, 1) / rho;
              const double v = states(i, 2) / rho;
              const double w = states(i, 3) / rho;
              const double a = states(i, 17);
              const double T = states(i, 18);
              const double p = states(i, 16);
              const double t = timepoint.count();
              const double x = probes(0, i);
              ::amrex::Print()
                  << fmt::format("{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< "
                                 "24.15g}{:< 24.15g}{:< 24.15g}{:< "
                                 "24.15g}{:< 24.15g}{:< 24.15g}\n",
                                 x, t, rho, u, v, w, a, T, p);
            }
          }
          ::amrex::Print() << "End Output for Probes.\n";
        }
      };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint(), 0);
  fub::RunOptions run_options = fub::GetRunOptions(vm);
  run_options.output_frequency.push_back(1);
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
