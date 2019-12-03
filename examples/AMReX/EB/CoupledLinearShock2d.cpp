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
#include "fub/AMReX_CutCell.hpp"
#include "fub/Solver.hpp"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_AllRegular.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <boost/filesystem.hpp>

#include <cmath>
#include <iostream>

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 2;

struct TemperatureRamp {
  using Complete = fub::IdealGasMix<Tube_Rank>::Complete;

  fub::IdealGasMix<Tube_Rank> equation_;

  void InitializeData(::amrex::MultiFab& data, const ::amrex::Geometry& geom) {
    fub::FlameMasterReactor& reactor = equation_.GetReactor();
    reactor.SetMoleFractions("N2:79,O2:21,H2:42");
    const double high_temp = 2000.0;
    const double low_temp = 300.0;
    Complete complete(equation_);
    Complete right(equation_);
    reactor.SetMoleFractions("N2:79,O2:21,H2:42");
    reactor.SetTemperature(low_temp);
    reactor.SetPressure(101325.0);
    equation_.CompleteFromReactor(right);

    fub::amrex::ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
        double x[AMREX_SPACEDIM] = {};
        geom.CellCenter(::amrex::IntVect{int(i)}, x);
        const double x0 = x[0] + 1.0;
        if (0.0 < x0 && x0 < 0.05) {
          const double d = std::clamp(x0 / 0.05, 0.0, 1.0);
          reactor.SetMoleFractions("N2:79,O2:21,H2:42");
          reactor.SetTemperature(d * low_temp + (1.0 - d) * high_temp);
          reactor.SetPressure(101325.0);
          equation_.CompleteFromReactor(complete);
          fub::Store(states, complete, {i});
        } else if (-0.05 < x0 && x0 < 0.0) {
          const double d = std::clamp(-x0 / 0.05, 0.0, 1.0);
          reactor.SetMoleFractions("N2:79,O2:21,H2:42");
          reactor.SetTemperature(d * low_temp + (1.0 - d) * high_temp);
          reactor.SetPressure(101325.0);
          equation_.CompleteFromReactor(complete);
          fub::Store(states, complete, {i});
        } else {
          fub::Store(states, right, {i});
        }
      });
    });
  }
};

auto MakeTubeSolver(int num_cells, int n_level, fub::Burke2012& mechanism) {
  const std::array<int, 2> n_cells{num_cells, 1};
  const std::array<double, 2> xlower{-1.5, -0.015};
  const std::array<double, 2> xupper{-0.03, +0.015};
  const std::array<int, 2> periodicity{0, 0};

  amrex::RealBox xbox(xlower, xupper);
  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};

  using namespace fub::amrex;

  CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  DataDescription desc = MakeDataDescription(equation);

  PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = n_level;
  hier_opts.refine_ratio = amrex::IntVect{2, 1};

  amrex::Geometry geom(amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1}}, &xbox,
                       -1, periodicity.data());
  geom.refine(hier_opts.refine_ratio);

  using Complete = fub::IdealGasMix<1>::Complete;
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 5e-3),
                            std::make_pair(&Complete::pressure, 5e-2)};

  ::amrex::Box refine_box{{num_cells - 5, 0}, {num_cells - 1, 0}};
  ConstantBox constant_box{refine_box};

  TemperatureRamp initial_data{equation};

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1}}};

  //  PatchHierarchy hierarchy =
  //  ReadCheckpointFile("/Volumes/Maikel_Intenso/FiniteVolumeSolver_Build/MultiBlock_2d/Checkpoint/Tube_00897",
  //  desc, geometry, hier_opts);
  PatchHierarchy hierarchy(desc, geometry, hier_opts);
  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      std::move(hierarchy), initial_data,
      TagAllOf(gradient, constant_box, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

   fub::EinfeldtSignalVelocities<fub::IdealGasMix<1>> signals{};
   fub::HllMethod hll_method(equation, signals);
  //   fub::MusclHancockMethod flux_method{equation, hll_method};
  // fub::ideal_gas::MusclHancockPrimMethod<1> flux_method(equation);

  HyperbolicMethod method{
      FluxMethod(fub::execution::simd, hll_method),
      ForwardIntegrator(fub::execution::simd),
      Reconstruction(fub::execution::simd, equation)};

  return fub::amrex::IntegratorContext(gridding, method);
}

auto Rectangle(const std::array<double, 2>& lower,
               const std::array<double, 2>& upper) {
  amrex::EB2::PlaneIF lower_x({lower[0], lower[1]}, {0, +1});
  amrex::EB2::PlaneIF lower_y({lower[0], lower[1]}, {+1, 0});
  amrex::EB2::PlaneIF upper_x({upper[0], upper[1]}, {0, -1});
  amrex::EB2::PlaneIF upper_y({upper[0], upper[1]}, {-1, 0});
  return amrex::EB2::makeIntersection(lower_x, lower_y, upper_x, upper_y);
}

::amrex::Box BoxWhichContains(const ::amrex::RealVect& xlower,
                              const ::amrex::RealVect& xupper,
                              const ::amrex::Geometry& geom) {
  ::amrex::Box domain = geom.Domain();
  ::amrex::IntVect lo = domain.smallEnd();
  ::amrex::IntVect up = domain.bigEnd();
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    for (int i = domain.smallEnd(d); i < domain.bigEnd(d); ++i) {
      const double x = geom.CellCenter(i, d);
      if (x < xlower[d]) {
        lo[d] = std::max(lo[d], i);
      }
      if (x > xupper[d]) {
        up[d] = std::min(up[d], i);
      }
    }
  }
  return ::amrex::Box{lo, up};
}

auto MakePlenumSolver(int num_cells, int n_level, fub::Burke2012& mechanism) {
  const std::array<int, Plenum_Rank> n_cells{num_cells, num_cells};
  const std::array<double, Plenum_Rank> xlower{-0.03, -0.50};
  const std::array<double, Plenum_Rank> xupper{+0.97, +0.50};
  const std::array<int, Plenum_Rank> periodicity{0, 0};

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1}},
                              &xbox, -1, periodicity.data());

  auto embedded_boundary =
      amrex::EB2::makeUnion(Rectangle({-1.0, +0.015}, {0.0, 1.0}),
                            Rectangle({-1.0, -1.0}, {0.0, -0.015}));
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};

  // Make Gridding Algorithm

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  equation.GetReactor().SetMoleFractions("N2:79,O2:21");
  equation.GetReactor().SetTemperature(300.0);
  equation.GetReactor().SetPressure(101325.0);
  fub::Complete<fub::IdealGasMix<Plenum_Rank>> right(equation);
  equation.CompleteFromReactor(right);

  using namespace fub::amrex::cutcell;

  fub::amrex::cutcell::RiemannProblem initial_data(
      equation, fub::Halfspace({+1.0, 0.0, 0.0}, -0.04), right, right);

  PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces = MakeIndexSpaces(shop, coarse_geom, n_level);

  using State = fub::Complete<fub::IdealGasMix<Plenum_Rank>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.1},
                             std::pair{&State::density, 0.1}};

  const ::amrex::Box refine_box =
      BoxWhichContains({-0.1, -0.015}, {0.02, +0.015}, coarse_geom);
  ConstantBox constant_box{refine_box};

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1}}};

  auto desc = fub::amrex::MakeDataDescription(equation);
  //  PatchHierarchy hierarchy =
  //  ReadCheckpointFile("/Volumes/Maikel_Intenso/FiniteVolumeSolver_Build/MultiBlock_2d/Checkpoint/Plenum_00897",
  //  desc, geometry, options);
  PatchHierarchy hierarchy(desc, geometry, options);
  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      std::move(hierarchy), initial_data,
      TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(4)),
      boundary_condition);
  gridding->InitializeHierarchy(0.0);

  // Make Solver

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::ideal_gas::MusclHancockPrimMethod<Plenum_Rank> flux_method(equation);
  //  fub::MusclHancockMethod flux_method{equation, hll_method};
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{
      FluxMethod{fub::execution::openmp_simd, cutcell_method},
      fub::amrex::cutcell::TimeIntegrator{},
      Reconstruction{fub::execution::openmp_simd, equation}};

  return fub::amrex::cutcell::IntegratorContext(gridding, method);
}

int main() {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard _{};

  fub::InitializeLogging(MPI_COMM_WORLD);

  fub::Burke2012 mechanism{};

  const int n_level = 2;

  const int num_cells = 128;

  auto plenum = MakePlenumSolver(num_cells, n_level, mechanism);
  auto tube = MakeTubeSolver(3 * num_cells / 2 - ((3 * num_cells / 2) % 8),
                             n_level, mechanism);

  fub::amrex::BlockConnection connection;
  connection.direction = fub::Direction::X;
  connection.side = 0;
  connection.plenum.id = 0;
  connection.plenum.mirror_box = plenum.GetGriddingAlgorithm()
                                     ->GetPatchHierarchy()
                                     .GetGeometry(0)
                                     .Domain();
  connection.tube.id = 0;
  connection.tube.mirror_box =
      tube.GetGriddingAlgorithm()->GetPatchHierarchy().GetGeometry(0).Domain();

  fub::amrex::MultiBlockIntegratorContext context(
      fub::FlameMasterReactor(mechanism), {tube}, {plenum}, {connection});

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};
  fub::DimensionalSplitLevelIntegrator system_solver(fub::int_c<Plenum_Rank>,
                                                     std::move(context));

  fub::amrex::MultiBlockKineticSouceTerm source_term{
      fub::IdealGasMix<Tube_Rank>{mechanism},
      system_solver.GetGriddingAlgorithm()};
  fub::SplitSystemSourceLevelIntegrator level_integrator{std::move(system_solver), source_term};

  fub::SubcycleFineFirstSolver solver(std::move(level_integrator));

  std::string base_name = "MultiBlock_2d";
  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};
  using namespace std::literals::chrono_literals;
  auto output = fub::MakeOutput<fub::amrex::MultiBlockGriddingAlgorithm>({}, {0.001s / 30.0},
      [&](const fub::amrex::MultiBlockGriddingAlgorithm& gridding) {
    std::ptrdiff_t cycle = gridding.GetCycles();
    ::amrex::Print() << "Checkpointing.\n";
    fub::amrex::WriteCheckpointFile(
        fmt::format("{}/Checkpoint/Tube_{:05}", base_name, cycle),
        gridding.GetTubes()[0]->GetPatchHierarchy());
    fub::amrex::cutcell::WriteCheckpointFile(
        fmt::format("{}/Checkpoint/Plenum_{:05}", base_name, cycle),
        gridding.GetPlena()[0]->GetPatchHierarchy());
    std::string name = fmt::format("{}/Tube/plt{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::WritePlotFile(
        name, gridding.GetTubes()[0]->GetPatchHierarchy(), tube_equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
    name = fmt::format("{}/Plenum/plt{:05}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::cutcell::WritePlotFile(
        name, gridding.GetPlena()[0]->GetPatchHierarchy(), equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
      });
  (*output)(*solver.GetGriddingAlgorithm());
  fub::RunOptions run_options{};
  run_options.final_time = 0.004s;
  run_options.cfl = 0.4;
  fub::RunSimulation(solver, run_options, wall_time_reference, *output);
}
