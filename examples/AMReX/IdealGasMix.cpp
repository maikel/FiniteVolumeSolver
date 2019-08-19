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

#include <fmt/format.h>
#include <iostream>

#include <xmmintrin.h>

namespace fub {
struct TemperatureRamp {
  using Complete = IdealGasMix<1>::Complete;

  IdealGasMix<1> equation_;

  void InitializeData(::amrex::MultiFab& data, const ::amrex::Geometry& geom) {
    FlameMasterReactor& reactor = equation_.GetReactor();
    reactor.SetMoleFractions("N2:79,O2:21,H2:42");
    const double high_temp = 1450.0;
    const double low_temp = 300.0;
    Complete complete(equation_);

    amrex::ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
      View<Complete> states =
          amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      ForEachIndex(Box<0>(states), [&](std::ptrdiff_t i) {
        double x[AMREX_SPACEDIM] = {};
        geom.CellCenter(::amrex::IntVect{int(i)}, x);
        if (x[0] < 0.03) {
          reactor.SetMoleFractions("N2:79,O2:21,H2:42");
          reactor.SetTemperature(high_temp);
          reactor.SetPressure(101325.0);
        } else if (x[0] < 0.7) {
          reactor.SetMoleFractions("N2:79,O2:21,H2:42");
          reactor.SetTemperature(low_temp);
          reactor.SetPressure(101325.0);
        } else {
          reactor.SetMoleFractions("N2:79,O2:21");
          reactor.SetTemperature(low_temp);
          reactor.SetPressure(101325.0);
        }
        equation_.CompleteFromReactor(complete);
        Store(states, complete, {i});
      });
    });
  }
};

void WriteMatlabFile(std::ostream& out,
                     const amrex::PatchHierarchy& hierarchy) {
  const amrex::PatchLevel& level = hierarchy.GetPatchLevel(0);
  const ::amrex::Box domain = hierarchy.GetGeometry(0).Domain();
  const ::amrex::BoxArray ba(domain);
  const ::amrex::DistributionMapping dm(ba);
  ::amrex::MultiFab mf(ba, dm,
                       hierarchy.GetDataDescription().n_state_components, 0);
  mf.ParallelCopy(level.data);
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    out.write(reinterpret_cast<char*>(mf[0].dataPtr()),
              static_cast<std::streamsize>(mf[0].size()) *
                  static_cast<std::streamsize>(sizeof(double)));
  }
}

} // namespace fub

int main(int argc, char** argv) {
  // Store a reference timepoint to measure the wall time duration
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  // Initialize the AMReX library with default parameters.
  const fub::amrex::ScopeGuard guard(argc, argv);

  // Enable floating point exceptions.
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_DIV_ZERO |
                         _MM_MASK_OVERFLOW | _MM_MASK_UNDERFLOW |
                         _MM_MASK_INVALID);

  constexpr int Dim = AMREX_SPACEDIM;

  // Setup the domain parameters
  const std::array<int, Dim> n_cells{AMREX_D_DECL(200, 1, 1)};
  const std::array<double, Dim> xlower{AMREX_D_DECL(0.0, 0.0, 0.0)};
  const std::array<double, Dim> xupper{AMREX_D_DECL(+1.0, +0.03, +0.03)};

  // Define the equation which will be solved
  fub::Burke2012 mechanism{};
  fub::IdealGasMix<1> equation{fub::FlameMasterReactor(mechanism)};

  // Define the GriddingAlgorithm for this simulation and initialize data.
  // {{{
  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  fub::TemperatureRamp initial_data{equation};

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
  boundary.conditions.push_back(TransmissiveBoundary{fub::Direction::X, 1});
  //  boundary.conditions.push_back(IsentropicBoundary{equation, 101325.0,
  //  fub::Direction::X, 1});

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 4;
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
  fub::MusclHancockMethod flux_method{equation, hll_method};

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

  std::string base_name = "Tube_AMR/";

  auto output =
      [&](const std::shared_ptr<fub::amrex::GriddingAlgorithm>& gridding,
          std::ptrdiff_t cycle, fub::Duration) {
        std::string name = fmt::format("{}plt{:05}", base_name, cycle);
        std::ofstream out(name + ".dat", std::ofstream::binary);
        amrex::Print() << "Start output to '" << name << "'.\n";
        fub::amrex::WritePlotFile(name, gridding->GetPatchHierarchy(),
                                  equation);
        fub::WriteMatlabFile(out, gridding->GetPatchHierarchy());
        amrex::Print() << "Finished output to '" << name << "'.\n";
      };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), 0, 0.0s);
  fub::RunOptions run_options{};
  run_options.cfl = 0.8;
  run_options.final_time = 0.0001s;
  run_options.output_interval = 4e-6s;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
