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

#include <boost/program_options.hpp>

#include <fmt/format.h>

#include <iostream>

#include <xmmintrin.h>

void WriteMatlabData(std::ostream& out, const amrex::FArrayBox& fab,
                     const fub::IdealGasMix<1>& eq,
                     const amrex::Geometry& geom) {
  using namespace fub;
  auto view = fub::amrex::MakeView<const Complete<IdealGasMix<1>>>(fab, eq);
  out << fmt::format("X Density VelocityX Temperature Pressure ");
  const int nspecies = eq.GetReactor().GetNSpecies();
  span<const std::string> names = eq.GetReactor().GetSpeciesNames();
  for (int s = 0; s < nspecies - 1; ++s) {
    out << names[s] << ' ';
  }
  out << names[nspecies - 1] << '\n';
  ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i) {
    double x[3] = {0.0, 0.0, 0.0};
    ::amrex::IntVect iv{int(i), 0, 0};
    geom.CellCenter(iv, x);
    const double density = view.density(i);
    const double velocity_x =
        density > 0.0 ? view.momentum(i, 0) / density : 0.0;
    const double temperature = view.temperature(i);
    const double pressure = view.pressure(i);
    out << fmt::format(
        "{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}", x[0],
        density, velocity_x, temperature, pressure);
    for (int s = 0; s < nspecies; ++s) {
      out << fmt::format("{:< 24.15e}", view.species(i, s));
    }
    out << '\n';
  });
}

struct TemperatureRamp {
  using Complete = fub::IdealGasMix<1>::Complete;

  fub::IdealGasMix<1> equation_;
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
        if (x[0] < ignition_pos_) {
          fub::Store(states, burnt, {i});
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
  double final_time{0.20};
  double cfl{0.8};
  int n_cells{200};
  int max_refinement_level{1};
  double domain_length{1.5};
  double ignition_position{0.6};
  double air_position{0.5 + 0.35};
  double equiv_ratio{1.0};
  double output_interval{1.0E-5};
};

namespace fub {

namespace amrex {
void WriteTubeData(std::ostream& out, const PatchHierarchy& hierarchy,
                   const IdealGasMix<1>& eq, fub::Duration time_point,
                   std::ptrdiff_t cycle_number, MPI_Comm comm) {
  const std::size_t n_level =
      static_cast<std::size_t>(hierarchy.GetNumberOfLevels());
  std::vector<::amrex::FArrayBox> fabs{};
  fabs.reserve(n_level);
  for (std::size_t level = 0; level < n_level; ++level) {
    const int ilvl = static_cast<int>(level);
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(ilvl);
    ::amrex::Box domain = level_geom.Domain();
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(ilvl).data;
    ::amrex::FArrayBox local_fab(domain, level_data.nComp());
    local_fab.setVal(0.0);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& patch_data = level_data[mfi];
      local_fab.copy(patch_data);
    });
    int rank = -1;
    ::MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      ::amrex::FArrayBox& fab = fabs.emplace_back(domain, level_data.nComp());
      fab.setVal(0.0);
      ::MPI_Reduce(local_fab.dataPtr(), fab.dataPtr(), local_fab.size(),
                   MPI_DOUBLE, MPI_SUM, 0, comm);
      if (level > 0) {
        for (int comp = 1; comp < level_data.nComp(); ++comp) {
          for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
            for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
              ::amrex::IntVect fine_i{i, j, domain.smallEnd(2)};
              ::amrex::IntVect coarse_i = fine_i;
              coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
              if (fab(fine_i, 0) == 0.0) {
                fab(fine_i, comp) = fabs[level - 1](coarse_i, comp);
              }
            }
          }
        }
        for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
          for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
            ::amrex::IntVect fine_i{i, j, domain.smallEnd(2)};
            ::amrex::IntVect coarse_i = fine_i;
            coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
            if (fab(fine_i, 0) == 0.0) {
              fab(fine_i, 0) = fabs[level - 1](coarse_i, 0);
            }
          }
        }
      }
      if (level == n_level - 1) {
        out << fmt::format("nx = {}\n", domain.length(0));
        out << fmt::format("t = {}\n", time_point.count());
        out << fmt::format("cycle = {}\n", cycle_number);
        WriteMatlabData(out, fab, eq, level_geom);
        out.flush();
      }
    } else {
      ::MPI_Reduce(local_fab.dataPtr(), nullptr, local_fab.size(), MPI_DOUBLE,
                   MPI_SUM, 0, comm);
    }
  }
}
} // namespace amrex

} // namespace fub

std::optional<ProgramOptions> ParseCommandLine(int argc, char** argv) {
  namespace po = boost::program_options;
  ProgramOptions opts;
  po::options_description desc("Program Options");
  desc.add_options()("help", "Print help messages")(
      "cfl", po::value<double>(&opts.cfl)->default_value(opts.cfl),
      "Set the CFL condition")("max_refinement_level",
                               po::value<int>(&opts.max_refinement_level)
                                   ->default_value(opts.max_refinement_level),
                               "Set the maximal refinement level")(
      "n_cells",
      po::value<int>(&opts.n_cells)->default_value(opts.n_cells),
      "Set number of cells in each direction for the plenum")(
      "final_time",
      po::value<double>(&opts.final_time)->default_value(opts.final_time),
      "Set the final simulation time")(
      "domain_length",
      po::value<double>(&opts.domain_length)
          ->default_value(opts.domain_length),
      "Set the base length for the tube")(
      "ignition_pos",
      po::value<double>(&opts.ignition_position)
          ->default_value(opts.ignition_position),
      "Set the position for the ignition of the detonation inside the tube")(
      "air_buffer_start",
      po::value<double>(&opts.air_position)
          ->default_value(opts.air_position),
      "Sets the starting position for an air buffer in the tube.")(
      "equiv_ratio",
      po::value<double>(&opts.equiv_ratio)
          ->default_value(opts.equiv_ratio),
      "Sets the equivalence ratio of the fuel in the tube")(
      "output_interval",
      po::value<double>(&opts.output_interval)
          ->default_value(opts.output_interval),
      "Sets the output interval");
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (std::exception& e) {
    amrex::Print()
        << "[Error] An Error occured while reading program options:\n";
    amrex::Print() << e.what() << '\n';
    return {};
  }

  if (vm.count("help")) {
    amrex::Print() << desc << "\n";
    return {};
  }

  opts.max_refinement_level = std::max(1, opts.max_refinement_level);

  amrex::Print() << fmt::format(
      "[Info] Simulation will run with the following options:\n[Info]\n");
  amrex::Print() << fmt::format("[Info] final_time = {}s\n", opts.final_time);
  amrex::Print() << fmt::format("[Info] output_interval = {}s\n",
                                opts.output_interval);
  amrex::Print() << fmt::format("[Info] cfl = {}\n", opts.cfl);
  amrex::Print() << fmt::format("[Info] max_refinement_level = {}\n",
                                opts.max_refinement_level);
  std::array<double, 3> xlo{0.0, 0.0, 0.0};
  std::array<double, 3> xup{opts.domain_length, +0.03, +0.03};
  amrex::Print() << fmt::format(
      "[Info] domain = [{}, {}] x [{}, {}] x [{}, {}]\n", xlo[0], xup[0],
      xlo[1], xup[1], xlo[2], xup[2]);
  amrex::Print() << fmt::format(
      "[Info] n_cells = {}\n[Info] dx = {}\n",
      opts.n_cells, opts.domain_length / opts.n_cells);
  amrex::Print() << fmt::format("[Info] ignition_pos = {}\n",
                                opts.ignition_position);
  amrex::Print() << fmt::format("[Info] equiv_ratio = {}\n",
                                opts.equiv_ratio);
  amrex::Print() << fmt::format("[Info] air_position = {}\n",
                                opts.air_position);

  return opts;
}

void MyMain(const ProgramOptions& opts) {
  // Store a reference timepoint to measure the wall time duration
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  // Enable floating point exceptions.
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO |
                         _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW |
                         _MM_MASK_INVALID);

  constexpr int Dim = AMREX_SPACEDIM;

  // Setup the domain parameters
  const std::array<int, Dim> n_cells{AMREX_D_DECL(opts.n_cells, 1, 1)};
  const int nlevels = opts.max_refinement_level;
  const std::array<double, Dim> xlower{AMREX_D_DECL(0.0, 0.0, 0.0)};
  const std::array<double, Dim> xupper{AMREX_D_DECL(opts.domain_length, +0.03, +0.03)};

  // Define the equation which will be solved
  fub::Burke2012 mechanism{};
  fub::IdealGasMix<1> equation{fub::FlameMasterReactor(mechanism)};

  // Define the GriddingAlgorithm for this simulation and initialize data.
  // {{{
  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  TemperatureRamp initial_data{equation};
  initial_data.equiv_raito_ = opts.equiv_ratio;
  initial_data.ignition_pos_ = opts.ignition_position;
  initial_data.air_position_ = opts.air_position;

  using Complete = fub::IdealGasMix<1>::Complete;
  fub::amrex::GradientDetector gradient{
      equation, std::make_pair(&Complete::pressure, 1e-3),
      std::make_pair(&Complete::density, 1e-3),
      std::make_pair(&Complete::temperature, 1e-1)};

  fub::amrex::BoundarySet boundary;
  using fub::amrex::IsentropicBoundary;
  using fub::amrex::ReflectiveBoundary;
  using fub::amrex::TransmissiveBoundary;
  boundary.conditions.push_back(
      ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::X, 0});
  //  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 1});
  boundary.conditions.push_back(
      IsentropicBoundary{equation, 101325.0, fub::Direction::X, 1});

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = nlevels;
  hier_opts.refine_ratio = ::amrex::IntVect{AMREX_D_DECL(2, 1, 1)};

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts), initial_data,
      gradient, boundary);
  gridding->InitializeHierarchy(0.0);
  // }}}

  // Setup the numerical Method used to solve this problem.
  // {{{
  fub::EinfeldtSignalVelocities<fub::IdealGasMix<1>> signals{};
  fub::HllMethod hll_method(equation, signals);

  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethod(fub::execution::simd, hll_method),
      fub::amrex::ForwardIntegrator(fub::execution::simd),
      fub::amrex::Reconstruction(fub::execution::simd, equation)};

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<1>, fub::amrex::IntegratorContext(gridding, method),
      fub::GodunovSplitting());

  fub::ideal_gas::KineticSourceTerm<1> source_term(equation, gridding);

  fub::DimensionalSplitSystemSourceSolver solver(system_solver, source_term,
                                                 fub::StrangSplitting());
  // }}}

  // Run the simulation with given feedback functions

  std::string base_name = "IdealGasMix/";

  auto output =
      [&](const std::shared_ptr<fub::amrex::GriddingAlgorithm>& gridding,
          std::ptrdiff_t cycle, fub::Duration timepoint, int) {
        std::string name = fmt::format("{}plt{:05}", base_name, cycle);
        std::ofstream out(name + ".dat");
        amrex::Print() << "Start output to '" << name << "'.\n";
        fub::amrex::WritePlotFile(name, gridding->GetPatchHierarchy(),
                                  equation);
        fub::amrex::WriteTubeData(out, gridding->GetPatchHierarchy(), equation,
                                  timepoint, cycle,
                                  solver.GetContext().GetMpiCommunicator());
        amrex::Print() << "Finished output to '" << name << "'.\n";
      };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), 0, 0.0s, 0);
  fub::RunOptions run_options{};
  run_options.cfl = opts.cfl;
  run_options.final_time = fub::Duration(opts.final_time);
  run_options.output_interval = {fub::Duration(opts.output_interval)};
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  {
    fub::amrex::ScopeGuard _{};
    std::optional<ProgramOptions> po = ParseCommandLine(argc, argv);
    if (po) {
      MyMain(*po);
    }
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}
