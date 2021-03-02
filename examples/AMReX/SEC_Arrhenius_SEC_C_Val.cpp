// Copyright (c) 2021 Maikel Nadolski
// Copyright (c) 2021 Christian Zenker
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
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/Solver.hpp"

#include "fub/equations/PerfectGasMix.hpp"
#include "fub/flux_method/MusclHancockMethod2.hpp"

#include "fub/AMReX/boundary_condition/GenericPressureValveBoundary.hpp"
#include "fub/AMReX/boundary_condition/IsentropicPressureExpansion.hpp"
#include "fub/AMReX/boundary_condition/TurbinePlenumBoundaryCondition.hpp"
#include "fub/equations/perfect_gas_mix/ArrheniusKinetics.hpp"
#include "fub/equations/perfect_gas_mix/PlenaControl.hpp"
#include "fub/equations/perfect_gas_mix/PlenaControlFeedback.hpp"
#include "fub/equations/perfect_gas_mix/PlenaControlOutput.hpp"

#include "fub/AMReX/AxialFluxMethodAdapter.hpp"
#include "fub/AMReX/AxialTimeIntegrator.hpp"
#include "fub/AMReX/DiffusionSourceTerm.hpp"

#include <cmath>
#include <fmt/format.h>
#include <iostream>

#include <fenv.h>

#include <boost/log/utility/manipulators/add_value.hpp>

#include <boost/filesystem.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace GT = fub::perfect_gas_mix::gt;

static constexpr int Tube_Rank = 1;

static_assert(AMREX_SPACEDIM == 2);

static constexpr double r_tube = 0.015;

static constexpr int n_species = 2;
static constexpr int n_passive_scalars = 0;

struct ComputeStableDt {
  fub::PerfectGasMix<1> eq;
  fub::Duration operator()(const fub::amrex::IntegratorContext& data, int level,
                           fub::Direction) const {
    const ::amrex::MultiFab& scratch = data.GetScratch(level);
    const ::amrex::Geometry& geom = data.GetGeometry(level);
    const double dx = geom.CellSize(0);
    double amax = 0.0;
    const double Minv = 1.0 / std::sqrt(eq.Msq);
    fub::IndexMapping<fub::PerfectGasMix<1>> c{eq};
    fub::amrex::ForEachFab(scratch, [&](const amrex::MFIter& mfi) {
      const amrex::Array4<const double> array = scratch[mfi].array();
      amrex::Dim3 lo = array.begin;
      amrex::Dim3 hi = array.end;
      for (int i = lo.x + scratch.nGrowVect()[0];
           i < hi.x - scratch.nGrowVect()[0]; ++i) {
        const double rho = array(i, 0, 0, c.density);
        const double rhou = array(i, 0, 0, c.momentum[0]);
        const double invrho = 1.0 / rho;
        const double u = std::abs(rhou * invrho) * Minv;
        const double p = array(i, 0, 0, c.pressure);
        // use exactly ruperts formula
        const double c = std::sqrt(eq.gamma * p * invrho);
        amax = std::max(amax, u + c * Minv);
      }
    });
    // fub::Duration dt(dx / amax);
    fub::Duration dt(5.0e-05 / 0.1125);
    return dt;
  }
};

struct InitialDataInTube {
  using Equation = fub::PerfectGasMix<1>;
  using Complete = fub::Complete<Equation>;
  using KineticState = fub::KineticState<Equation>;
  using Conservative = fub::Conservative<Equation>;

  Equation equation_;
  fub::perfect_gas_mix::ArrheniusKineticsOptions kinetics_;
  fub::perfect_gas_mix::gt::ControlOptions control_;
  double x_0_;
  double initially_filled_x_{0.4};

  double CompressorPressureRatio(double rpm) const noexcept {
    double pratio = control_.pratiomean +
                    control_.c_0 * control_.pratiovar *
                        std::atan(M_PI * (rpm / control_.rpmmean - 1.0)) / M_PI;
    return pratio;
  }

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration /*time*/) const {
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& data = patch_level.data;
    fub::amrex::ForEachFab(fub::execution::seq, data, [&](amrex::MFIter& mfi) {
      fub::Conservative<fub::PerfectGasMix<1>> state(equation_);
      Complete complete(equation_);
      // fub::Array<double, 1, 1> velocity{0.0};
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
        const double x = geom.CellCenter(int(i), 0);
        const double rel_x = x - x_0_;
        // const double pressure = 1.0;

        const double pV = CompressorPressureRatio(control_.rpmmin);
        const double TV =
            1.0 + (std::pow(pV, equation_.gamma_minus_one_over_gamma) - 1.0) /
                      control_.efficiency_compressor;
        const double p0 = 2.0;
        const double p = 0.95 * p0;
        const double T = TV + kinetics_.Q * equation_.gamma_minus_one;
        const double rho = p / T; //* equation_.ooRspec;

        state.density = rho;
        state.species[0] = (rel_x < initially_filled_x_) ? 1.0 : 0.0;
        state.energy = p * equation_.gamma_minus_one_inv;
        equation_.CompleteFromCons(complete, state);
        fub::Store(states, complete, {i});
      });
    });
  }
};

using FactoryFunction =
    std::function<fub::AnyFluxMethod<fub::amrex::IntegratorContext>(
        const fub::PerfectGasMix<1>&)>;

template <typename... Pairs> auto GetFluxMethodFactory(Pairs... ps) {
  std::map<std::string, FactoryFunction> factory;
  ((factory[ps.first] = ps.second), ...);
  return factory;
}

template <typename FluxMethod> struct MakeFlux {
  fub::AnyFluxMethod<fub::amrex::IntegratorContext>
  operator()(const fub::PerfectGasMix<1>& eq) const {
    FluxMethod flux_method{eq};
    fub::amrex::FluxMethodAdapter adapter(std::move(flux_method));
    return adapter;
  }
};

void WriteCheckpoint(
    const std::string& path, const fub::amrex::GriddingAlgorithm& grid,
    const std::shared_ptr<const fub::perfect_gas_mix::gt::ControlState>&
        control_state) {
  std::string name = fmt::format("{}/Tube", path);
  fub::amrex::WriteCheckpointFile(name, grid.GetPatchHierarchy());
  int rank = -1;
  ::MPI_Comm_rank(amrex::ParallelDescriptor::Communicator(), &rank);
  if (rank == 0) {
    name = fmt::format("{}/ControlState", path);
    std::ofstream valve_checkpoint(name);
    boost::archive::text_oarchive oa(valve_checkpoint);
    oa << *control_state;
  }
}

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard scope_guard{};

  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::ProgramOptions eq_options = fub::GetOptions(options, "Equation");
  double Rspec = fub::GetOptionOr(eq_options, "Rspec", 1.0);
  double gamma = fub::GetOptionOr(eq_options, "gamma", 1.28);
  fub::PerfectGasConstants constants{Rspec, gamma};

  BOOST_LOG(log) << "Equation:";
  BOOST_LOG(log) << fmt::format("  - Rspec = {}", constants.Rspec);
  BOOST_LOG(log) << fmt::format("  - gamma = {}", constants.gamma);

  fub::PerfectGasMix<Tube_Rank> tube_equation{constants, n_species - 1,
                                              n_passive_scalars};
  fub::ProgramOptions control_options_map =
      fub::GetOptions(options, "ControlOptions");
  GT::ControlOptions control_options(control_options_map);
  BOOST_LOG(log) << "ControlOptions:";
  control_options.Print(log);
  GT::Control control(tube_equation, control_options);
  std::shared_ptr<GT::ControlState> control_state = control.GetSharedState();
  std::string checkpoint =
      fub::GetOptionOr(control_options_map, "checkpoint", std::string{});
  if (!checkpoint.empty()) {
    std::string input =
        fub::ReadAndBroadcastFile(checkpoint + "/ControlState",
                                  ::amrex::ParallelDescriptor::Communicator());
    std::istringstream ifs(input);
    boost::archive::text_iarchive ia(ifs);
    ia >> *control_state;
  }

  fub::ProgramOptions tube_options = fub::GetOptions(options, "Tubes");

  // Make tube
  BOOST_LOG(log) << "==================== Tube =========================";
  using namespace fub::amrex;

  CartesianGridGeometry grid_geometry(
      fub::GetOptions(tube_options, "GridGeometry"));
  BOOST_LOG(log) << "GridGeometry:";
  grid_geometry.Print(log);

  fub::ProgramOptions hier_opts =
      fub::GetOptions(tube_options, "PatchHierarchy");
  PatchHierarchyOptions hierarchy_options(hier_opts);
  BOOST_LOG(log) << "PatchHierarchy:";
  hierarchy_options.Print(log);

  using Complete = fub::PerfectGasMix<1>::Complete;
  GradientDetector gradient{tube_equation,
                            std::make_pair(&Complete::density, 5e-3),
                            std::make_pair(&Complete::pressure, 5e-2)};

  DataDescription desc = MakeDataDescription(tube_equation);

  ::amrex::Box refine_box{{grid_geometry.cell_dimensions[0] - 5, 0},
                          {grid_geometry.cell_dimensions[0] - 1, 0}};
  ConstantBox constant_box{refine_box};

  fub::perfect_gas_mix::ArrheniusKinetics<1> arrhenius_source_term{
      tube_equation, fub::GetOptions(tube_options, "ArrheniusKinetics")};
  BOOST_LOG(log) << "ArrheniusKinetics:";
  arrhenius_source_term.options.Print(log);
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());

  const double initially_filled_x =
      fub::GetOptionOr(tube_options, "initially_filled_x", 0.4);
  BOOST_LOG(log) << "InitialData:";
  BOOST_LOG(log) << "  - initially_filled_x = " << initially_filled_x << " [m]";
  InitialDataInTube initial_data{
      tube_equation, arrhenius_source_term.options, control_options,
      grid_geometry.coordinates.lo()[0], initially_filled_x};

  FUB_ASSERT(control_state);
  fub::ProgramOptions SEC_options =
      fub::GetOptions(tube_options, "InflowOptionsSEC");
  double SEC_buffer = fub::GetOptionOr(SEC_options, "SEC_buffer", 0.06);
  double SEC_tti = fub::GetOptionOr(SEC_options, "SEC_tti", 1.2);
  double SEC_timin = fub::GetOptionOr(SEC_options, "SEC_timin", 0.1);
  BOOST_LOG(log) << "InflowOptionsSEC:";
  BOOST_LOG(log) << "  - SEC_buffer = " << SEC_buffer;
  BOOST_LOG(log) << "  - SEC_tti = " << SEC_tti;
  BOOST_LOG(log) << "  - SEC_timin = " << SEC_timin;

  // const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  auto switch_to_SEC_inflow_function =
      [eps, arrhenius_source_term, SEC_buffer, SEC_tti, SEC_timin,
       prim = fub::Primitive<fub::PerfectGasMix<1>>(tube_equation)](
          const fub::PerfectGasMix<1>& eq,
          fub::Complete<fub::PerfectGasMix<1>>& boundary_state,
          const fub::perfect_gas_mix::gt::PlenumState& compressor_state,
          double inner_pressure, fub::Duration t_diff, const amrex::MultiFab&,
          const fub::amrex::GriddingAlgorithm&, int) mutable {
        if (!compressor_state.SEC_Mode) {
          // Deflagration Mode
          const double X_inflow_left = 1.0;

          const double p_inflow_left = compressor_state.pressure;
          const double T_inflow_left = compressor_state.temperature;
          const double rho_inflow_left =
              p_inflow_left / T_inflow_left;// * eq.ooRspec;

          const double p = inner_pressure;
          const double ppv = p_inflow_left;
          const double rhopv = rho_inflow_left;
          const double Tpv = T_inflow_left;
          const double pin = p;
          const double Tin =
              Tpv * pow(pin / ppv, eq.gamma_minus_one_over_gamma);
          const double uin = std::sqrt(2.0 * eq.gamma_over_gamma_minus_one *
                                       std::max(0.0, Tpv - Tin));
          double rhoin = rhopv * pow(pin / ppv, eq.gamma_inv);

          FUB_ASSERT(rhoin > 0.0);
          FUB_ASSERT(pin > 0.0);
          prim.density = rhoin;
          prim.velocity[0] = uin;
          prim.pressure = pin;
          prim.species[0] = X_inflow_left;

          fub::CompleteFromPrim(eq, boundary_state, prim);
        } else {
          // SEC Mode
          const double fuel_retardatation =
              SEC_buffer; /* Reference: 0.06;   for icx = 256:  0.1*/
          const double tti = SEC_tti;     /* Reference: 0.75; best 1.2 */
          const double timin = SEC_timin; /* Reference: 0.1 */
          const double X_inflow_left = 1.0;

          const double p_inflow_left = compressor_state.pressure;
          const double T_inflow_left = compressor_state.temperature;
          const double rho_inflow_left =
              p_inflow_left / T_inflow_left;// * eq.ooRspec;

          const double p = inner_pressure;
          const double ppv = p_inflow_left;
          const double rhopv = rho_inflow_left;
          const double Tpv = T_inflow_left;
          const double pin = p;
          const double Tin =
              Tpv * pow(pin / ppv, eq.gamma_minus_one_over_gamma);
          const double uin = std::sqrt(2.0 * eq.gamma_over_gamma_minus_one *
                                       std::max(0.0, Tpv - Tin));
          double rhoin = rhopv * pow(pin / ppv, eq.gamma_inv);

          const double tign = std::max(timin, tti - t_diff.count());
          const double Xin = X_inflow_left;
          const double Tin1 = fub::perfect_gas_mix::TemperatureForIgnitionDelay(
              eq, arrhenius_source_term.options, tign, Xin, Tin, eps);

          /* adjust density to match desired temperature */
          rhoin *= Tin / Tin1;

          prim.density = rhoin;
          prim.velocity[0] = uin;
          prim.pressure = pin;
          auto heaviside = [](double x) { return (x > 0); };
          prim.species[0] = std::clamp(
              Xin * heaviside(t_diff.count() - fuel_retardatation), 0.0, 1.0);

          fub::CompleteFromPrim(eq, boundary_state, prim);
        }
      };

  fub::amrex::SwitchDeflagrationToSECValve valve(tube_equation, control_state,
                                                 switch_to_SEC_inflow_function);

  fub::amrex::TurbinePlenumBoundaryCondition plenum_bc{tube_equation,
                                                       control_state};
  BoundarySet boundaries{{valve, plenum_bc}};

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
      // gridding->GetPatchHierarchy().SetCounterRegistry(counters);
      gridding->InitializeHierarchy(0.0);
      return gridding;
    } else {
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

  auto [flux_method, time_integrator] =
      fub::amrex::GetFluxMethod(fub::GetOptions(options, "FluxMethod"),
                                gridding->GetPatchHierarchy(), tube_equation);
  HyperbolicMethod method{flux_method, time_integrator,
                          Reconstruction(tube_equation)};

  IntegratorContext context(gridding, method,
                            fub::GetOptions(tube_options, "IntegratorContext"));
  BOOST_LOG(log) << "IntegratorContext:";
  context.GetOptions().Print(log);
  context.SetComputeStableDt(ComputeStableDt{tube_equation});

  fub::amrex::DiffusionSourceTermOptions diff_opts =
      fub::GetOptions(options, "DiffusionSourceTerm");
  fub::amrex::DiffusionSourceTerm<fub::PerfectGasMix<1>> diffusion_source_term{
      tube_equation, diff_opts};
  // const fub::Duration good_guess_dt =
  //     diffusion_source_term.ComputeStableDt(context, 0);
  // TracePassiveScalarBoundary passive_scalar_boundary{tube_equation,
  //                                                    good_guess_dt};
  // context.GetGriddingAlgorithm()->GetBoundaryCondition() =
  //     BoundarySet{{valve, passive_scalar_boundary}};

  BOOST_LOG(log) << "==================== End Tube =========================";

  GT::ControlFeedbackOptions feedback_options =
      fub::GetOptions(options, "ControlFeedback");
  BOOST_LOG(log) << "ControlFeedback:";
  feedback_options.Print(log);
  GT::ControlFeedback<Tube_Rank> feedback(
      tube_equation, control, feedback_options); // TODO make feedback function
  context.SetPostAdvanceHierarchyFeedback(feedback);

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<Tube_Rank>, std::move(context), fub::GodunovSplitting{});

  // fub::SplitSystemSourceLevelIntegrator diffusive_integrator(
  //     std::move(level_integrator), std::move(diffusion_source_term),
  //     fub::StrangSplittingLumped());

  // fub::SplitSystemSourceLevelIntegrator reactive_integrator(
  //     std::move(diffusive_integrator), std::move(arrhenius_source_term),
  //     fub::StrangSplittingLumped());

  fub::SplitSystemSourceLevelIntegrator reactive_integrator(
      std::move(level_integrator), std::move(arrhenius_source_term),
      fub::StrangSplittingLumped());

  fub::SubcycleFineFirstSolver solver(std::move(reactive_integrator));

  // using namespace fub::amrex;
  using namespace std::literals::chrono_literals;

  using namespace fub::amrex;
  struct MakeCheckpoint
      : public fub::OutputAtFrequencyOrInterval<GriddingAlgorithm> {
    std::string directory_ = "SEC_Plenum_Arrhenius";
    std::shared_ptr<const fub::perfect_gas_mix::gt::ControlState>
        control_state_;
    MakeCheckpoint(
        const fub::ProgramOptions& options,
        std::shared_ptr<const fub::perfect_gas_mix::gt::ControlState> cs)
        : OutputAtFrequencyOrInterval(options), control_state_{cs} {
      directory_ = fub::GetOptionOr(options, "directory", directory_);
    }
    void operator()(const GriddingAlgorithm& grid) override {
      std::string name = fmt::format("{}/{:09}", directory_, grid.GetCycles());
      fub::SeverityLogger log = fub::GetInfoLogger();
      BOOST_LOG(log) << fmt::format("Write Checkpoint to '{}'!", name);
      WriteCheckpoint(name, grid, control_state_);
    }
  };

  fub::OutputFactory<fub::amrex::GriddingAlgorithm> factory{};
  factory.RegisterOutput<MakeCheckpoint>("Checkpoint", control_state);
  factory.RegisterOutput<GT::ControlOutput>("ControlOutput", control_state);
  factory.RegisterOutput<fub::amrex::WriteHdf5>(
      "HDF5", fub::VarNames<fub::Complete<fub::PerfectGasMix<1>>,
                            std::vector<std::string>>(tube_equation));
  using CounterOutput = fub::CounterOutput<fub::amrex::GriddingAlgorithm,
                                           std::chrono::nanoseconds>;
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> outputs(
      std::move(factory), fub::GetOptions(options, "Output"));

  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);

  outputs(*solver.GetGriddingAlgorithm());
  fub::RunSimulation(solver, run_options, wall_time_reference, outputs);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
  pybind11::scoped_interpreter interpreter{};
  std::optional<fub::ProgramOptions> opts = fub::ParseCommandLine(argc, argv);
  if (opts) {
    fub::InitializeLogging(MPI_COMM_WORLD,
                           fub::GetOptions(*opts, "LogOptions"));
    MyMain(*opts);
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}