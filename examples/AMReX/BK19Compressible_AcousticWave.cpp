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
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/Solver.hpp"
#include <AMReX_MLMG.H>

#include "fub/AMReX/CompressibleAdvectionIntegratorContext.hpp"
#include "fub/AMReX/MLMG/MLNodeHelmDualCstVel.hpp"
#include "fub/AMReX/solver/BK19LevelIntegrator.hpp"
#include "fub/equations/CompressibleAdvection.hpp"

struct AcousticWaveInitialData : fub::amrex::BK19PhysicalParameters {
  using Complete = fub::CompressibleAdvection<2>::Complete;
  AcousticWaveInitialData() {}

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration /*time*/) const {
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& mf = patch_level.data;
    fub::amrex::ForEachFab(mf, [&](const amrex::MFIter& mfi) {
      fub::CompressibleAdvection<2> equation{};
      amrex::FArrayBox& fab = mf[mfi];
      const amrex::Box& box = mfi.tilebox();
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(fab, equation, box);
      fub::ForEachIndex(fub::Box<0>(states), [&](int i, int j) {
        const double x = geom.CellCenter(i, 0);

        const double p = std::pow(1.0 + del0 * std::sin(wn * x),
                                  2.0 * gamma / (gamma - 1.0));
        const double rho = std::pow(p, 1.0 / gamma);
        const double c = std::sqrt(gamma * p / rho);
        const double Ma = std::sqrt(Msq);

        states.density(i, j) = rho;
        const double velocity0 = U0[0] + (p - 1.0) / (rho * c) / Ma;
        const double velocity1 = U0[1];
        states.PTdensity(i, j) = rho;

        states.momentum(i, j, 0) =
            states.density(i, j) * velocity0;
        states.momentum(i, j, 1) =
            states.density(i, j) * velocity1;
        states.PTinverse(i, j) = states.density(i, j) / states.PTdensity(i, j);
      });
    });

    // set initial values of pi
    amrex::MultiFab& pi = *patch_level.nodes;
    const double Gamma = (gamma - 1.0) / gamma;
    fub::amrex::ForEachFab(pi, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& fab = pi[mfi];
      fub::amrex::ForEachIndex(fab.box(), [&](auto... is) {
        ::amrex::IntVect i{int(is)...};

        ::amrex::Vector<double> coor(2);
        geom.LoNode(i, coor);
        const double x = coor[0];

        const double p = std::pow(1.0 + del0 * std::sin(wn * x),
                                  2.0 * gamma / (gamma - 1.0));
        fab(i, 0) = (pow(p, Gamma) - 1.0) / Msq;
      });
    });
  }

  const double del0 = 0.05;
  const double wn = 2.0 * M_PI;
  std::array<double, 2> U0{1.0, 0.0};
};

void MyMain(const fub::ProgramOptions& options) {
  using namespace fub::amrex;
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  ScopeGuard amrex_scope_guard{};

  const double h_ref{1.0};
  const double t_ref{1.0};
  const double T_ref{353.048780488};
  const double u_ref{h_ref / t_ref};

  // Here, some things are dimensional and others non-dimensionalized. Adjust???
  AcousticWaveInitialData inidat;
  inidat.R_gas = 287.0;
  inidat.gamma = 2.0;
  inidat.Msq = u_ref * u_ref / (inidat.R_gas * T_ref);
  inidat.c_p = inidat.gamma / (inidat.gamma - 1.0);
  inidat.alpha_p = 1.0;

  DataDescription desc{};
  desc.n_state_components = 5;
  desc.n_cons_components = 4;
  desc.n_node_components = 1;

  CartesianGridGeometry grid_geometry =
      fub::GetOptions(options, "GridGeometry");
  PatchHierarchyOptions hierarchy_options =
      fub::GetOptions(options, "PatchHierarchy");

  fub::SeverityLogger info = fub::GetInfoLogger();
  BOOST_LOG(info) << "GridGeometry:";
  grid_geometry.Print(info);
  BOOST_LOG(info) << "PatchHierarchy:";
  hierarchy_options.Print(info);

  PatchHierarchy hierarchy(desc, grid_geometry, hierarchy_options);

  using Complete = fub::CompressibleAdvection<2>::Complete;
  fub::CompressibleAdvection<2> equation{};

  fub::IndexMapping<fub::CompressibleAdvection<2>> index(equation);

  GradientDetector gradient(equation, std::pair{&Complete::PTinverse, 1.0e-2});

  std::shared_ptr grid = std::make_shared<GriddingAlgorithm>(
      std::move(hierarchy), inidat, TagAllOf{gradient, TagBuffer(2)});
  grid->InitializeHierarchy(0.0);

  // set number of MG levels to 1 (effectively no MG)
  amrex::LPInfo lp_info;
  lp_info.setMaxCoarseningLevel(0);

  auto box_array = grid->GetPatchHierarchy().GetPatchLevel(0).box_array;
  auto dmap = grid->GetPatchHierarchy().GetPatchLevel(0).distribution_mapping;
  auto linop = std::make_shared<amrex::MLNodeHelmDualCstVel>(
      amrex::Vector<amrex::Geometry>{grid->GetPatchHierarchy().GetGeometry(0)},
      amrex::Vector<amrex::BoxArray>{box_array},
      amrex::Vector<amrex::DistributionMapping>{dmap}, lp_info);

  linop->setDomainBC(
      {AMREX_D_DECL(amrex::LinOpBCType::Periodic, amrex::LinOpBCType::Periodic,
                    amrex::LinOpBCType::Periodic)},
      {AMREX_D_DECL(amrex::LinOpBCType::Periodic, amrex::LinOpBCType::Periodic,
                    amrex::LinOpBCType::Periodic)});

  fub::CompressibleAdvectionFluxMethod<2> flux_method{};

  HyperbolicMethod method{flux_method, EulerForwardTimeIntegrator(),
                          Reconstruction(fub::execution::seq, equation)};

  //   CompressibleAdvectionIntegratorContext simulation_data(grid, method, 2,
  //   0);
  CompressibleAdvectionIntegratorContext simulation_data(grid, method, 4, 2);

  fub::DimensionalSplitLevelIntegrator advection(
      //       fub::int_c<2>, std::move(simulation_data),
      //       fub::GodunovSplitting());
      fub::int_c<2>, std::move(simulation_data), fub::StrangSplitting());

  BK19LevelIntegratorOptions integrator_options =
      fub::GetOptions(options, "BK19LevelIntegrator");
  BOOST_LOG(info) << "BK19LevelIntegrator:";
  integrator_options.Print(info);
  BK19LevelIntegrator level_integrator(equation, std::move(advection), linop,
                                       inidat, integrator_options);
  fub::NoSubcycleSolver solver(std::move(level_integrator));

  CompressibleAdvectionAdvectiveFluxes& Pv =
      solver.GetContext().GetAdvectiveFluxes(0);
  RecomputeAdvectiveFluxes(index, Pv.on_faces, Pv.on_cells,
                           solver.GetContext().GetScratch(0),
                           solver.GetContext().GetGeometry(0).periodicity());

  using namespace std::literals::chrono_literals;
  std::string base_name = "BK19_CompAcousticWave/";

  fub::OutputFactory<GriddingAlgorithm> factory;
  factory.RegisterOutput<fub::AnyOutput<GriddingAlgorithm>>(
      "Plotfile", WriteBK19Plotfile{base_name});
  factory.RegisterOutput<fub::amrex::DebugOutput>(
      "DebugOutput",
      solver.GetGriddingAlgorithm()->GetPatchHierarchy().GetDebugStorage());
  fub::MultipleOutputs<GriddingAlgorithm> output{
      std::move(factory), fub::GetOptions(options, "Output")};

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(info) << "RunOptions:";
  run_options.Print(info);

  if (integrator_options.do_initial_projection) {
    solver.GetLevelIntegrator().InitialProjection(0);
  }

  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
  pybind11::scoped_interpreter interpreter{};
  {
    std::optional<fub::ProgramOptions> vm = fub::ParseCommandLine(argc, argv);
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
