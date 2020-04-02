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

struct CircleData {
  using Complete = fub::Complete<fub::Advection2d>;

  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration /*time*/) const {
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& data = patch_level.data;
    fub::amrex::ForEachFab(data, [&](const amrex::MFIter& mfi) {
      const ::amrex::Box& box = mfi.tilebox();
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], fub::Advection2d{{}}, box);
      fub::ForEachIndex(
          fub::Box<0>(states), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
            double x[AMREX_SPACEDIM] = {};
            geom.CellCenter(amrex::IntVect{AMREX_D_DECL(int(i), int(j), 0)}, x);
            const double norm = std::sqrt(x[0] * x[0] + x[1] * x[1]);
            if (norm < 0.25) {
              states.mass(i, j) = 3.0;
            } else {
              states.mass(i, j) = 1.0;
            }
          });
    });
  }
};

void MyMain(const fub::ProgramOptions& opts) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard guard{};

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM >= 2);

  fub::Advection2d equation{{1.0, 1.0}};
  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::amrex::CartesianGridGeometry geometry =
      fub::GetOptions(opts, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  geometry.Print(log);

  fub::amrex::PatchHierarchyOptions hier_opts =
      fub::GetOptions(opts, "PatchHierarchy");
  BOOST_LOG(log) << "PatchHierarchy:";
  hier_opts.Print(log);

  using State = fub::Advection2d::Complete;
  fub::amrex::GradientDetector gradient{equation,
                                        std::pair{&State::mass, 1e-3}};

  std::shared_ptr grid = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts), CircleData{},
      fub::amrex::TagAllOf(gradient, fub::amrex::TagBuffer(2)));
  grid->InitializeHierarchy(0.0);

  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethodAdapterAdapter(fub::GodunovMethod{equation}),
      fub::amrex::EulerForwardTimeIntegrator(), fub::amrex::NoReconstruction{}};

  const int scratch_gcw = 2;
  const int flux_gcw = 1;

  fub::DimensionalSplitLevelIntegrator level_integrator(
      fub::int_c<Dim>,
      fub::amrex::IntegratorContext(grid, method, scratch_gcw, flux_gcw),
      fub::GodunovSplitting());

  // fub::NoSubcycleSolver solver(std::move(level_integrator));
  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  using namespace std::literals::chrono_literals;
  using Plotfile = fub::amrex::PlotfileOutput<fub::Advection2d>;
  using CounterOutput = fub::CounterOutput<fub::amrex::GriddingAlgorithm,
                                           std::chrono::milliseconds>;
  fub::OutputFactory<fub::amrex::GriddingAlgorithm> factory{};
  factory.RegisterOutput<Plotfile>("Plotfile", equation);
  factory.RegisterOutput<CounterOutput>("CounterOutput", wall_time_reference);
  factory.RegisterOutput<fub::amrex::WriteHdf5>("HDF5");
  fub::MultipleOutputs<fub::amrex::GriddingAlgorithm> output(
      std::move(factory), fub::GetOptions(opts, "Output"));

  output(*grid);
  fub::RunOptions run_options = fub::GetOptions(opts, "RunOptions");
  BOOST_LOG(log) << "RunOptions:";
  run_options.Print(log);
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
