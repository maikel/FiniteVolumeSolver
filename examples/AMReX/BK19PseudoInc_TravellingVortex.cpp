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

double p_coeff(double r, const std::vector<double>& coefficients) {
  if (r >= 1.0) {
    return 0.0;
  }

  double result = 0.0;
  int exponent = 12;
  for (double c : coefficients) {
    result += c * (std::pow(r, exponent) - 1.0);
    exponent += 1;
  }
  return result;
}

struct InitialData {
  using Complete = fub::CompressibleAdvection<2>::Complete;
  InitialData() {
    coefficients.resize(25);
    coefficients[0] = 1.0 / 12.0;
    coefficients[1] = -12.0 / 13.0;
    coefficients[2] = 9.0 / 2.0;
    coefficients[3] = -184.0 / 15.0;
    coefficients[4] = 609.0 / 32.0;
    coefficients[5] = -222.0 / 17.0;
    coefficients[6] = -38.0 / 9.0;
    coefficients[7] = 54.0 / 19.0;
    coefficients[8] = 783.0 / 20.0;
    coefficients[9] = -558.0 / 7.0;
    coefficients[10] = 1053.0 / 22.0;
    coefficients[11] = 1014.0 / 23.0;
    coefficients[12] = -1473.0 / 16.0;
    coefficients[13] = 204.0 / 5.0;
    coefficients[14] = 510.0 / 13.0;
    coefficients[15] = -1564.0 / 27.0;
    coefficients[16] = 153.0 / 8.0;
    coefficients[17] = 450.0 / 29.0;
    coefficients[18] = -269.0 / 15.0;
    coefficients[19] = 174.0 / 31.0;
    coefficients[20] = 57.0 / 32.0;
    coefficients[21] = -74.0 / 33.0;
    coefficients[22] = 15.0 / 17.0;
    coefficients[23] = -6.0 / 35.0;
    coefficients[24] = 1.0 / 72.0;
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
        const double dx = x - center[0];
        const double dy = y - center[1];
        const double r = std::sqrt(dx * dx + dy * dy);

        states.PTdensity(i, j) = 1.0;

        if (r < R0) {
          const double r_over_R0 = r / R0;
          const double uth =
              fac * std::pow(1.0 - r_over_R0, 6) * std::pow(r_over_R0, 6);

          states.density(i, j) =
              rho0 + del_rho * std::pow(1.0 - r_over_R0 * r_over_R0, 6);
          states.velocity(i, j, 0) = U0[0] - uth * (dy / r);
          states.velocity(i, j, 1) = U0[1] + uth * (dx / r);

          // pi(i, j) = Gamma * fac*fac * p_coeff(r_over_R0, coefficients);

        } else {
          states.density(i, j) = rho0;
          states.velocity(i, j, 0) = U0[0];
          states.velocity(i, j, 1) = U0[1];

          //           pi(i, j) = 0.0
        }
        states.momentum(i, j, 0) =
            states.density(i, j) * states.velocity(i, j, 0);
        states.momentum(i, j, 1) =
            states.density(i, j) * states.velocity(i, j, 1);
        states.PTinverse(i, j) = states.density(i, j) / states.PTdensity(i, j);
      });
    });
  }

  std::vector<double> coefficients;
  const double a_rho{1.0};
  const double rho0{a_rho * 0.5};
  const double del_rho{a_rho * 0.5};
  const double R0{0.4};
  const double fac{1024.0};
  std::array<double, 2> center{0.5, 0.5};
  std::array<double, 2> U0{1.0, 1.0};
};

void MyMain(const fub::ProgramOptions& options) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::DataDescription desc{};
  desc.n_state_components = 7;
  desc.n_cons_components = 4;
  desc.n_node_components = 1;

  fub::ProgramOptions grid_options =
      fub::ToMap(fub::GetOptionOr(options, "GridGeometry", pybind11::dict()));
  fub::amrex::CartesianGridGeometry geometry{grid_options};

  fub::ProgramOptions hierarchy_options = fub::ToMap(
      fub::GetOptionOr(options, "PatchHierarchyOptions", pybind11::dict()));
  fub::amrex::PatchHierarchyOptions hier_opts(hierarchy_options);

  fub::amrex::PatchHierarchy hierarchy(desc, geometry, hier_opts);

  using Complete = fub::CompressibleAdvection<2>::Complete;
  fub::CompressibleAdvection<2> equation{};
  fub::IndexMapping<fub::CompressibleAdvection<2>> index(equation);

  fub::amrex::GradientDetector gradient(
      equation, std::pair{&Complete::PTinverse, 1.0e-2});

  const InitialData inidat{};

  std::shared_ptr grid = std::make_shared<fub::amrex::GriddingAlgorithm>(
      std::move(hierarchy), inidat,
      fub::amrex::TagAllOf{gradient, fub::amrex::TagBuffer(2)});
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

  using namespace fub;
  CompressibleAdvectionFluxMethod<2> flux_method{};

  auto tag = fub::execution::seq;
  fub::amrex::HyperbolicMethod method{
      flux_method, fub::amrex::ForwardIntegrator(tag),
      fub::amrex::Reconstruction(tag, equation)};

  fub::amrex::BK19IntegratorContext simulation_data(grid, method, 2, 0);
  const int nlevel = simulation_data.GetPatchHierarchy().GetNumberOfLevels();

  // set initial values of pi
  const double gamma{1.4};
  const double Gamma = (gamma - 1.0) / gamma;
  // simulation_data.GetPi(0).setVal(1.0);

  for (int level = 0; level < nlevel; ++level) {
    ::amrex::MultiFab& pi = simulation_data.GetPi(level);
    const ::amrex::Geometry& geom =
        grid->GetPatchHierarchy().GetGeometry(level);
    fub::amrex::ForEachFab(pi, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& fab = pi[mfi];
      fub::amrex::ForEachIndex(fab.box(), [&](auto... is) {
        ::amrex::IntVect i{int(is)...};

        ::amrex::Vector<double> coor(2);
        geom.LoNode(i, coor);
        const double dx = coor[0] - inidat.center[0];
        const double dy = coor[1] - inidat.center[1];
        const double r = std::sqrt(dx * dx + dy * dy);

        if (r < inidat.R0) {
          const double r_over_R0 = r / inidat.R0;
          fab(i, 0) = Gamma * inidat.fac * inidat.fac *
                      p_coeff(r_over_R0, inidat.coefficients);

        } else {
          fab(i, 0) = 0.0;
        }
      });
    });
  }

  fub::DimensionalSplitLevelIntegrator advection(
      fub::int_c<2>, std::move(simulation_data), fub::GodunovSplitting());
//       fub::int_c<2>, std::move(simulation_data), fub::StrangSplitting());

  fub::ProgramOptions integrator_options = fub::ToMap(
      fub::GetOptionOr(options, "BK19LevelIntegrator", pybind11::dict()));
  fub::amrex::BK19LevelIntegrator level_integrator(
      equation, std::move(advection), linop, integrator_options);

  fub::amrex::BK19AdvectiveFluxes& Pv =
      level_integrator.GetContext().GetAdvectiveFluxes(0);
  // Pv.on_faces[0].setVal(0.0);
  // Pv.on_faces[1].setVal(0.0);
  fub::amrex::RecomputeAdvectiveFluxes(
      index, Pv.on_faces, Pv.on_cells,
      level_integrator.GetContext().GetScratch(0),
      level_integrator.GetContext().GetGeometry(0).periodicity());

  fub::NoSubcycleSolver solver(std::move(level_integrator));

  using namespace std::literals::chrono_literals;
  std::string base_name = "BK19_Pseudo_Incompressible/";

  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> output{};
  output.AddOutput(fub::MakeOutput<fub::amrex::GriddingAlgorithm>(
      {1}, {0.01s}, fub::amrex::WriteBK19Plotfile{base_name}));

  output(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.final_time = 1.0s;
  run_options.cfl = 0.45;
  // run_options.max_cycles = 10;
  fub::RunSimulation(solver, run_options, wall_time_reference, output);
}

std::optional<fub::ProgramOptions> ParseCommandLine(int argc, char** argv) {
  namespace po = boost::program_options;
  po::options_description desc{};
  std::string config_path{};
  desc.add_options()("config", po::value<std::string>(&config_path),
                     "Path to the config file which can be parsed.");
  po::variables_map vm;
  fub::ProgramOptions options{};
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("config")) {
      config_path = vm["config"].as<std::string>();
      options = fub::ParsePythonScript(config_path, MPI_COMM_WORLD);
    }
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
  return options;
}

int main(int argc, char** argv) {
  MPI_Init(nullptr, nullptr);
  fub::InitializeLogging(MPI_COMM_WORLD);
  pybind11::scoped_interpreter interpreter{};
  {
    fub::amrex::ScopeGuard _{};
    std::optional<fub::ProgramOptions> vm = ParseCommandLine(argc, argv);
    boost::log::sources::severity_logger<boost::log::trivial::severity_level>
        log(boost::log::keywords::severity = boost::log::trivial::info);
    if (vm) {
      fub::ProgramOptions integrator_options = fub::ToMap(
          fub::GetOptionOr(*vm, "BK19LevelIntegrator", pybind11::dict()));
      BOOST_LOG(log) << "BK19LevelIntegrator Options:";
      fub::amrex::BK19LevelIntegratorOptions(integrator_options).Print(log);

      fub::ProgramOptions hierarchy_options = fub::ToMap(
          fub::GetOptionOr(*vm, "PatchHierarchyOptions", pybind11::dict()));
      BOOST_LOG(log) << "PatchHierarchy Options:";
      fub::amrex::PatchHierarchyOptions(hierarchy_options).Print(log);

      fub::ProgramOptions grid_options =
          fub::ToMap(fub::GetOptionOr(*vm, "GridGeometry", pybind11::dict()));
      BOOST_LOG(log) << "CartesianGridGeometry Options:";
      fub::amrex::CartesianGridGeometry(grid_options).Print(log);
      MyMain(*vm);
    }
  }
  int flag = -1;
  MPI_Finalized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}
