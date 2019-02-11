#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/PatchHierarchy.hpp"

#include "fub/equations/Advection.hpp"
#include "fub/geometry/Ball.hpp"
#include "fub/initial_data/RiemannProblem.hpp"

#include <AMReX.H>
#include <AMReX_Print.H>

#include "boost/iostreams/device/null.hpp"
#include "boost/iostreams/stream.hpp"

void RunSimulation() {
  constexpr int kDimension = AMREX_SPACEDIM;

  const std::array<int, kDimension> n_cells{32, 32, 32};
  const std::array<double, kDimension> lower{0.0, 0.0, 0.0};
  const std::array<double, kDimension> upper{1.0, 1.0, 1.0};
  // const std::chrono::duration<double> final_time = 1s;

  const std::array<double, kDimension> velocity{1.0, 0.0, 0.0};
  const fub::Advection equation{velocity};
  const fub::amrex::CartesianGridGeometry geometry{n_cells, {lower, upper}};

  fub::amrex::PatchHierarchy hierarchy(
      equation, geometry, {.ghost_layer_width = 1, .max_refinement_level = 1});

  // fub::GradientDetector tagging(&State::mass, 1e-2);
  // fub::amrex::GriddingAlgorithm gridding(hierarchy, tagging);
  fub::amrex::GriddingAlgorithm gridding(hierarchy);

  using State = fub::SingleState<fub::Advection<kDimension>::State>;
  fub::amrex::InitializeHierarchy(
      gridding,
      fub::RiemannProblem(State{1.0}, State{2.0}, fub::Ball(lower, 0.1)));

  // fub::HyperbolicSplitTimeIntegrator time_integrator(equation);
  // fub::MusclHancockMethod flux_method(equation,
  // fub::GodunovMethod(equation)); fub::StrangSplitting splitting{};

  // fub::amrex::HyperbolicSplitSystemSolver solver(splitting, time_integrator,
  //                                                flux_method);
  // solver.AdvanceHierarchy(gridding, final_time);

  // fub::amrex::AmrexOutput writer(equation);
  // writer.plot("advection", hierarchy);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  boost::iostreams::stream<boost::iostreams::null_sink> null{
      boost::iostreams::null_sink()};
  amrex::Initialize(MPI_COMM_WORLD, null, null);
  // amrex::Initialize(MPI_COMM_WORLD);
  try {
    RunSimulation();
  } catch (std::exception& e) {
    std::cerr << "An unexpected error happened: " << e.what() << '\n';
  }
  amrex::Finalize();
  MPI_Finalize();
}