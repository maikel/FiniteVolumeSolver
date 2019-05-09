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

#include "fub/equations/IdealGasMix.hpp"
#include "fub/CartesianCoordinates.hpp"
#include "fub/equations/ideal_gas_mix/KineticSourceTerm.hpp"
#include "fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp"

#include "fub/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/HyperbolicSplitSystemSolver.hpp"
#include "fub/SplitSystemSourceSolver.hpp"
#include "fub/boundary_condition/TransmissiveBoundary.hpp"

#include "fub/ext/Eigen.hpp"
#include "fub/flux_method/HllMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"

#include "fub/AMReX/FluxMethod.hpp"
#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/HyperbolicSplitIntegratorContext.hpp"
#include "fub/AMReX/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/AMReX/Reconstruction.hpp"
#include "fub/AMReX/ScopeGuard.hpp"

#include "fub/split_method/StrangSplitting.hpp"
#include "fub/tagging/GradientDetector.hpp"
#include "fub/tagging/TagBuffer.hpp"

#include "fub/RunSimulation.hpp"

#include <fmt/format.h>
#include <iostream>

#include <xmmintrin.h>

struct RiemannProblem {
  fub::IdealGasMix<2> equation_;
  void
  InitializeData(const fub::View<fub::Complete<fub::IdealGasMix<2>>>& states,
                 const fub::amrex::PatchHierarchy& hierarchy,
                 fub::amrex::PatchHandle patch) {
    const ::amrex::Geometry& geom = hierarchy.GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->tilebox();
    fub::CartesianCoordinates x =
        fub::amrex::GetCartesianCoordinates(geom, box);
    fub::FlameMasterReactor& reactor = equation_.GetReactor();
    reactor.SetMoleFractions("N2:79,O2:21,H2:42");
    const double high_temp = 1100.0;
    const double low_temp = 300.0;
    fub::Complete<fub::IdealGasMix<2>> complete(equation_);
    fub::ForEachIndex(fub::Box<0>(states), [&](auto... is) {
      if (x(is...)[0] < 0.01) {
        const double x0 = x(is...)[0];
        const double d = x0 / 0.01;
        reactor.SetTemperature(d * low_temp + (1.0 - d) * high_temp);
        reactor.SetPressure(101325.0);
      } else {
        reactor.SetTemperature(low_temp);
        reactor.SetPressure(101325.0);
      }
      equation_.CompleteFromReactor(complete);
      Store(states, complete, {is...});
    });
  }
};

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  const fub::amrex::ScopeGuard guard(argc, argv);

  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO |
                         _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW |
                         _MM_MASK_INVALID);

  constexpr int Dim = AMREX_SPACEDIM;

  const std::array<int, Dim> n_cells{AMREX_D_DECL(32, 32, 1)};
  const std::array<double, Dim> xlower{AMREX_D_DECL(0.0, 0.0, 0.0)};
  const std::array<double, Dim> xupper{AMREX_D_DECL(+0.1, +0.1, +0.1)};

  fub::Burke2012 mechanism{};
  fub::IdealGasMix<2> equation{fub::FlameMasterReactor(mechanism)};

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 1;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 2, 1)};

  using Complete = fub::IdealGasMix<2>::Complete;
  fub::GradientDetector gradient{equation,
                                 std::make_pair(&Complete::density, 5e-3),
                                 std::make_pair(&Complete::pressure, 5e-2)};

  RiemannProblem initial_data{equation};

  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(desc, geometry, hier_opts),
      fub::amrex::AdaptInitialData(initial_data, equation),
      fub::amrex::AdaptTagging(equation, gradient),
      fub::TransmissiveBoundary(equation));
  gridding->InitializeHierarchy(0.0);

  fub::HyperbolicSplitPatchIntegrator patch_integrator{equation};
  fub::HllMethod base_method(
      equation, fub::EinfeldtSignalVelocities<fub::IdealGasMix<2>>);
  fub::MusclHancockMethod flux_method{equation, base_method};

  const int gcw = base_method.GetStencilWidth();
  fub::HyperbolicSplitSystemSolver system_solver(
      fub::HyperbolicSplitLevelIntegrator(
          fub::amrex::HyperbolicSplitIntegratorContext(gridding, gcw),
          fub::amrex::HyperbolicSplitPatchIntegrator(patch_integrator),
          fub::amrex::FluxMethod(flux_method),
          fub::amrex::Reconstruction(equation)));

  fub::ideal_gas::KineticSourceTerm<2> source_term(equation, gridding);

  fub::SplitSystemSourceSolver solver(system_solver, source_term);

  std::string base_name = "IdealGasMix_1d_embed_in_Nd/";

  auto output = [&](const fub::amrex::PatchHierarchy& hierarchy,
                    std::ptrdiff_t cycle, fub::Duration) {
    std::string name = fmt::format("{}{:05}", base_name, cycle);
    amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::WritePlotFile(name, hierarchy, equation);
    amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  auto print_msg = [](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(solver.GetPatchHierarchy(), 0, 0.0s);
  fub::RunOptions run_options{};
  run_options.cfl = 0.8;
  run_options.final_time = 0.002s;
  run_options.output_interval = 0.00001s;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}
