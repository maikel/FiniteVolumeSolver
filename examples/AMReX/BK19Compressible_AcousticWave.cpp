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

#include "fub/AMReX/MLMG/MLNodeHelmDualCstVel.hpp"
#include "fub/AMReX/bk19/BK19IntegratorContext.hpp"
#include "fub/AMReX/bk19/BK19LevelIntegrator.hpp"
#include "fub/equations/CompressibleAdvection.hpp"

struct AcousticWaveInitialData {
  using Complete = fub::CompressibleAdvection<2>::Complete;
  AcousticWaveInitialData() {
  }

  void InitializeData(amrex::MultiFab& mf, const amrex::Geometry& geom) const {
    fub::amrex::ForEachFab(mf, [&](const amrex::MFIter& mfi) {
      fub::CompressibleAdvection<2> equation{};
      amrex::FArrayBox& fab = mf[mfi];
      const amrex::Box& box = mfi.tilebox();
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(fab, equation, box);
      fub::ForEachIndex(fub::Box<0>(states), [&](int i, int j) {
        const double x = geom.CellCenter(i, 0);
        const double y = geom.CellCenter(j, 1);
        const double dx = x;
        const double dy = y;

        states.density(i, j)     = rho0;
        states.velocity(i, j, 0) = U0[0];
        states.velocity(i, j, 1) = U0[1];
        states.PTdensity(i, j)   = 1.0;

        states.momentum(i, j, 0) =
            states.density(i, j) * states.velocity(i, j, 0);
        states.momentum(i, j, 1) =
            states.density(i, j) * states.velocity(i, j, 1);
        states.PTinverse(i, j) = states.density(i, j) / states.PTdensity(i, j);
      });
    });
  }

  const double rho0{1.0};
  const double R_gas{287.4};
  const double gamma{1.4};
  const double h_ref{100.0};
  const double t_ref{100.0};
  const double T_ref{300.0};
  const double u_ref{h_ref / t_ref};
  const double Msq{u_ref * u_ref / (R_gas * T_ref)};
  std::array<double, 2> U0{0.0, 0.0};
};

void MyMain(const fub::ProgramOptions& options) {
  using namespace fub::amrex;
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard amrex_scope_guard{};

  fub::amrex::DataDescription desc{};
  desc.n_state_components = 7;
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

  const AcousticWaveInitialData inidat{};
  const double Gamma = (inidat.gamma - 1.0) / inidat.gamma;

  using Complete = fub::CompressibleAdvection<2>::Complete;
  fub::CompressibleAdvection<2> equation{};
  // Here, this is not c_p but (gamma - 1) / gamma! Rename/Adjust???
  equation.c_p = Gamma;
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

  BK19IntegratorContext simulation_data(grid, method, 4, 2);
  const int nlevel = simulation_data.GetPatchHierarchy().GetNumberOfLevels();

  // set initial values of pi
  for (int level = 0; level < nlevel; ++level) {
    ::amrex::MultiFab& pi = simulation_data.GetPi(level);
    const ::amrex::Geometry& geom =
        grid->GetPatchHierarchy().GetGeometry(level);
    ForEachFab(pi, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& fab = pi[mfi];
      ForEachIndex(fab.box(), [&](auto... is) {
        ::amrex::IntVect i{int(is)...};

        ::amrex::Vector<double> coor(2);
        geom.LoNode(i, coor);
        const double dx = coor[0];
        const double dy = coor[1];

        fab(i, 0) = 0.0;
      });
    });
  }

  fub::DimensionalSplitLevelIntegrator advection(
      fub::int_c<2>, std::move(simulation_data), fub::StrangSplitting());

  BK19LevelIntegratorOptions integrator_options =
      fub::GetOptions(options, "BK19LevelIntegrator");
  BOOST_LOG(info) << "BK19LevelIntegrator:";
  integrator_options.Print(info);
  BK19LevelIntegrator level_integrator(equation, std::move(advection), linop,
                                       integrator_options);

  BK19AdvectiveFluxes& Pv = level_integrator.GetContext().GetAdvectiveFluxes(0);
  // Pv.on_faces[0].setVal(0.0);
  // Pv.on_faces[1].setVal(0.0);
  RecomputeAdvectiveFluxes(
      index, Pv.on_faces, Pv.on_cells,
      level_integrator.GetContext().GetScratch(0),
      level_integrator.GetContext().GetGeometry(0).periodicity());

  fub::NoSubcycleSolver solver(std::move(level_integrator));

  using namespace std::literals::chrono_literals;
  std::string base_name = "BK19_CompAcousticWave/";

  fub::OutputFactory<GriddingAlgorithm> factory;
  factory.RegisterOutput<fub::AnyOutput<GriddingAlgorithm>>("Plotfile", WriteBK19Plotfile{base_name});
  fub::MultipleOutputs<GriddingAlgorithm> output{std::move(factory), fub::GetOptions(options, "Output")};

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options = fub::GetOptions(options, "RunOptions");
  BOOST_LOG(info) << "RunOptions:";
  run_options.Print(info);
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
