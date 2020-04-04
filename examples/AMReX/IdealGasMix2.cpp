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

struct TemperatureRamp {
  using Complete = fub::IdealGasMix<1>::Complete;

  TemperatureRamp(const fub::IdealGasMix<1>& equation,
                  const fub::ProgramOptions& opts)
      : equation_(equation) {
    ignition_pos_ = fub::GetOptionOr(opts, "ignition_pos", ignition_pos_);
    ramp_width_ = fub::GetOptionOr(opts, "ramp_width", ramp_width_);
    high_temperature_ =
        fub::GetOptionOr(opts, "high_temperature", high_temperature_);
    fill_fraction_ = fub::GetOptionOr(opts, "fill_fraction", fill_fraction_);
    using namespace std::literals;
    fub::ProgramOptions left = fub::GetOptions(opts, "left");
    left_moles_ = fub::GetOptionOr(left, "moles", left_moles_);
    double temperature = fub::GetOptionOr(left, "temperature", 300.0);
    double pressure = fub::GetOptionOr(left, "pressure", 101325.0);
    double velocity = fub::GetOptionOr(left, "velocity", 0.0);
    auto& reactor = equation_.GetReactor();
    reactor.SetMoleFractions(left_moles_);
    reactor.SetTemperature(temperature);
    reactor.SetPressure(pressure);
    equation_.CompleteFromReactor(left_, Eigen::Array<double, 1, 1>{velocity});
    equation_.CompleteFromCons(left_, left_);

    fub::ProgramOptions right = fub::GetOptions(opts, "right");
    right_moles_ = fub::GetOptionOr(right, "moles", right_moles_);
    temperature = fub::GetOptionOr(right, "temperature", 300.0);
    pressure = fub::GetOptionOr(right, "pressure", 101325.0);
    velocity = fub::GetOptionOr(right, "velocity", 0.0);
    reactor.SetMoleFractions(right_moles_);
    reactor.SetTemperature(temperature);
    reactor.SetPressure(pressure);
    equation_.CompleteFromReactor(right_, Eigen::Array<double, 1, 1>{velocity});
    equation_.CompleteFromCons(right_, right_);
  }

  template <typename Log> void Print(Log& log) {
    BOOST_LOG(log) << fmt::format(" - ignition_pos = {} [m]", ignition_pos_);
    BOOST_LOG(log) << fmt::format(" - ramp_width = {} [m]", ramp_width_);
    BOOST_LOG(log) << fmt::format(" - high_temperature = {} [K]",
                                  high_temperature_);
    BOOST_LOG(log) << fmt::format(" - left.moles = {} [-]", left_moles_);
    BOOST_LOG(log) << fmt::format(" - left.temperature = {} [K]",
                                  left_.temperature);
    BOOST_LOG(log) << fmt::format(" - left.pressure = {} [Pa]", left_.pressure);
    BOOST_LOG(log) << fmt::format(" - left.velocity = {} [m/s]",
                                  left_.momentum[0] / left_.density);
    BOOST_LOG(log) << fmt::format(" - right.moles = {} [-]", right_moles_);
    BOOST_LOG(log) << fmt::format(" - right.temperature = {} [K]",
                                  right_.temperature);
    BOOST_LOG(log) << fmt::format(" - right.pressure = {} [Pa]",
                                  right_.pressure);
    BOOST_LOG(log) << fmt::format(" - right.velocity = {} [m/s]",
                                  right_.momentum[0] / right_.density);
    BOOST_LOG(log) << fmt::format(" - fill_fraction = {} [m]", fill_fraction_);
  }

  fub::IdealGasMix<1> equation_;
  double ignition_pos_{-1.5};
  double air_position_{0.0};
  double ramp_width_{0.05};
  double high_temperature_{2000.0};
  double fill_fraction_{0.7};
  std::string left_moles_{"N2:79,O2:21,H2:42"};
  std::string right_moles_{"N2:79,O2:21"};
  Complete left_{equation_};
  Complete right_{equation_};
  Complete state_{equation_};

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration /* t */) {
    ::amrex::MultiFab& data = patch_level.data;
    const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    fub::amrex::ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
        double x = geom.CellCenter(i, 0);
        if (x < fill_fraction_) {
          fub::Store(states, left_, {i});
        } else {
          fub::Store(states, right_, {i});
        }
      });
    });

    fub::amrex::ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
        double x = geom.CellCenter(i, 0);
        if (ignition_pos_ - ramp_width_ < x && x < ignition_pos_) {
          fub::Load(state_, states, {i});
          equation_.SetReactorStateFromComplete(state_);
          const double lambda = std::clamp(
              (x - ignition_pos_ + ramp_width_) / ramp_width_, 0.0, 1.0);
          const double T_soll =
              (1.0 - lambda) * state_.temperature + lambda * high_temperature_;
          equation_.GetReactor().SetTemperature(T_soll);
          equation_.CompleteFromReactor(state_,
                                        state_.momentum / state_.density);
          equation_.CompleteFromCons(state_, state_);
          fub::Store(states, state_, {i});
        } else if (ignition_pos_ < x && x < ignition_pos_ + ramp_width_) {
          fub::Load(state_, states, {i});
          equation_.SetReactorStateFromComplete(state_);
          const double lambda =
              std::clamp((x - ignition_pos_) / ramp_width_, 0.0, 1.0);
          const double T_soll =
              (1.0 - lambda) * high_temperature_ + lambda * state_.temperature;
          equation_.GetReactor().SetTemperature(T_soll);
          equation_.CompleteFromReactor(state_,
                                        state_.momentum / state_.density);
          equation_.CompleteFromCons(state_, state_);
          fub::Store(states, state_, {i});
        }
      });
    });
  }
};

void MyMain(const fub::ProgramOptions& opts) {
  // Store a reference timepoint to measure the wall time duration
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard scope_guard{};

  using namespace fub::amrex;

  constexpr int TubeDim = 1;

  // Define the equation which will be solved
  fub::Burke2012 mechanism{};
  fub::IdealGasMix<TubeDim> equation{fub::FlameMasterReactor(mechanism)};

  fub::SeverityLogger info = fub::GetInfoLogger();

  // Define the GriddingAlgorithm for this simulation and initialize data.
  // {{{
  fub::amrex::CartesianGridGeometry grid_geometry =
      fub::GetOptions(opts, "CartesianGridGeometry");
  BOOST_LOG(info) << "CartesianGridGeometry: ";
  grid_geometry.Print(info);

  PatchHierarchyOptions hierarchy_options =
      fub::GetOptions(opts, "PatchHierarchy");
  BOOST_LOG(info) << "PatchHierarchy: ";
  hierarchy_options.Print(info);

  using Complete = fub::IdealGasMix<1>::Complete;
  fub::amrex::GradientDetector gradient{
      equation, std::make_pair(&Complete::pressure, 1e-3),
      std::make_pair(&Complete::density, 1e-3),
      std::make_pair(&Complete::temperature, 1e-1)};

  BoundarySet boundary;

  IsentropicPressureBoundaryOptions right_boundary =
      fub::GetOptions(opts, "IsentropicPressureBoundary");
  BOOST_LOG(info) << "IsentropicPressureBoundary: ";
  right_boundary.Print(info);

  boundary.conditions.push_back(
      ReflectiveBoundary{equation, fub::Direction::X, 0});
  boundary.conditions.push_back(
      IsentropicPressureBoundary{equation, right_boundary});

  TemperatureRamp initial_data(equation,
                               fub::GetOptions(opts, "TemperatureRamp"));
  BOOST_LOG(info) << "TemperatureRamp: ";
  initial_data.Print(info);

  std::shared_ptr gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, grid_geometry, hierarchy_options),
      initial_data, gradient, boundary);
  gridding->InitializeHierarchy(0.0);
  // }}}

  // Setup the numerical Method used to solve this problem.
  // {{{
  fub::ideal_gas::MusclHancockPrimMethod<1> flux_method(equation);

  fub::amrex::HyperbolicMethod method{fub::amrex::FluxMethodAdapter(flux_method),
                                      fub::amrex::EulerForwardTimeIntegrator(),
                                      fub::amrex::Reconstruction(equation)};

  const int scratch_gcw = 2 * flux_method.GetStencilWidth();
  const int flux_gcw = flux_method.GetStencilWidth();

  fub::DimensionalSplitLevelIntegrator system_solver(
      fub::int_c<1>,
      fub::amrex::IntegratorContext(gridding, method, scratch_gcw, flux_gcw),
      fub::GodunovSplitting());

  // fub::amrex::IgniteDetonationOptions io{};
  // io.ignite_interval = fub::Duration(0.01);
  // fub::amrex::IgniteDetonation ignite(equation,
  // hier_opts.max_number_of_levels,
  //                                     io);

  // fub::SplitSystemSourceLevelIntegrator ign_solver(
  //     std::move(system_solver), std::move(ignite), fub::GodunovSplitting{});

  fub::ideal_gas::KineticSourceTerm<1> source_term(equation);

  fub::SplitSystemSourceLevelIntegrator level_integrator(
      std::move(system_solver), std::move(source_term),
      fub::StrangSplittingLumped());

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));
  // }}}

  // Run the simulation with given feedback functions

  using namespace std::literals::chrono_literals;
  using Plotfile = PlotfileOutput<fub::IdealGasMix<1>>;
  using CounterOutput =
      fub::CounterOutput<GriddingAlgorithm, std::chrono::milliseconds>;
  fub::OutputFactory<GriddingAlgorithm> factory{};
  factory.RegisterOutput<Plotfile>("Plotfile", equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<WriteHdf5>("HDF5");
  fub::MultipleOutputs<GriddingAlgorithm> output(
      std::move(factory), fub::GetOptions(opts, "Output"));

  fub::RunOptions run_options = fub::GetOptions(opts, "RunOptions");
  BOOST_LOG(info) << "RunOptions: ";
  run_options.Print(info);

  output(*solver.GetGriddingAlgorithm());
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
