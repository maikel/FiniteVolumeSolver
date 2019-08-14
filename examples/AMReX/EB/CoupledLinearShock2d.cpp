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
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>
#include <AMReX_EB2_IF_AllRegular.H>

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
    const double high_temp = 1350.0;
    const double low_temp = 300.0;
    Complete complete(equation_);

    fub::amrex::ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
        double x[AMREX_SPACEDIM] = {};
        geom.CellCenter(::amrex::IntVect{int(i)}, x);
        if (x[0] < 0.1) {
          const double d = x[0] / 0.1;
          reactor.SetMoleFractions("N2:79,O2:21,H2:42");
          reactor.SetTemperature(d * low_temp + (1.0 - d) * high_temp);
          reactor.SetPressure(101325.0);
        } else if (x[0] < 0.8) {
          reactor.SetMoleFractions("N2:79,O2:21,H2:42");
          reactor.SetTemperature(low_temp);
          reactor.SetPressure(101325.0);
        } else {
          reactor.SetMoleFractions("N2:79,O2:21");
          reactor.SetTemperature(low_temp);
          reactor.SetPressure(101325.0);
        }
        equation_.CompleteFromReactor(complete);
        fub::Store(states, complete, {i});
      });
    });
  }
};

auto MakeTubeSolver(int num_cells, fub::Burke2012& mechanism) {
  const std::array<int, 2> n_cells{num_cells, 1};
  const std::array<double, 2> xlower{0.0, 0.0};
  const std::array<double, 2> xupper{+1.0, +0.03};
  const std::array<int, 2> periodicity{0, 0};

  amrex::RealBox xbox(xlower, xupper);
  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};

  using namespace fub::amrex;

  CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

  DataDescription desc = MakeDataDescription(equation);

  PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 1;
  hier_opts.refine_ratio = amrex::IntVect{2, 1};

  amrex::Geometry geom(amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1}}, &xbox, -1, periodicity.data());
  geom.refine(hier_opts.refine_ratio);
  ::amrex::EB2::Build(::amrex::EB2::makeShop(::amrex::EB2::AllRegularIF()), geom, hier_opts.refine_ratio, 1, 1);

  using Complete = fub::IdealGasMix<1>::Complete;
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 5e-3),
                            std::make_pair(&Complete::pressure, 5e-2)};

  ::amrex::Box refine_box{{num_cells - 5, 0}, {num_cells - 1, 0}};
  ConstantBox constant_box{refine_box};

  TemperatureRamp initial_data{equation};

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(PatchHierarchy(desc, geometry, hier_opts), initial_data, TagAllOf(gradient, constant_box),
      boundary_condition);
  gridding->InitializeHierarchy(0.0);

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<1>> signals{};
  fub::HllMethod hll_method(equation, signals);
  //fub::MusclHancockMethod flux_method{equation, hll_method};

  HyperbolicMethod method{FluxMethod(fub::execution::seq, hll_method),
                          ForwardIntegrator(fub::execution::seq),
                          Reconstruction(fub::execution::seq, equation)};

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

auto MakePlenumSolver(int num_cells, fub::Burke2012& mechanism) {
  const std::array<int, Plenum_Rank> n_cells{num_cells, num_cells};
  const std::array<double, Plenum_Rank> xlower{-0.10, -0.15};
  const std::array<double, Plenum_Rank> xupper{+0.20, +0.15};
  const std::array<int, Plenum_Rank> periodicity{0, 0};

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1}},
                              &xbox, -1, periodicity.data());

  const int n_level = 1;

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
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

  ::amrex::Box refine_box{{0, 0}, {5, num_cells - 1}};
  ConstantBox constant_box{refine_box};

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, options), initial_data,
      TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  // Make Solver

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::MusclHancockMethod flux_method(equation, hll_method);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{fub::execution::seq, cutcell_method},
                          fub::amrex::cutcell::TimeIntegrator{},
                          Reconstruction{fub::execution::seq, equation}};

  return fub::amrex::cutcell::IntegratorContext(gridding, method);
}

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard _(argc, argv);
  fub::Burke2012 mechanism{};
  auto plenum = MakePlenumSolver(256, mechanism);
  auto tube = MakeTubeSolver(1600, mechanism);

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
  fub::HyperbolicSplitSystemSolver system_solver(
      fub::HyperbolicSplitLevelIntegrator(equation, std::move(context)));
  fub::amrex::MultiBlockKineticSouceTerm source_term{
      fub::IdealGasMix<1>{mechanism}, system_solver.GetGriddingAlgorithm()};
  fub::SplitSystemSourceSolver solver{system_solver, source_term};

  std::string base_name = "MultiBlock_2d";
  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};
  auto output =
      [&](std::shared_ptr<fub::amrex::MultiBlockGriddingAlgorithm> gridding,
          auto cycle, auto) {
        std::string name = fmt::format("{}/Tube/plt{:05}", base_name, cycle);
        ::amrex::Print() << "Start output to '" << name << "'.\n";
        fub::amrex::WritePlotFile(
            name, gridding->GetTubes()[0]->GetPatchHierarchy(), tube_equation);
        ::amrex::Print() << "Finished output to '" << name << "'.\n";
        name = fmt::format("{}/Plenum/plt{:05}", base_name, cycle);
        ::amrex::Print() << "Start output to '" << name << "'.\n";
        fub::amrex::cutcell::WritePlotFile(
            name, gridding->GetPlena()[0]->GetPatchHierarchy(), equation);
        ::amrex::Print() << "Finished output to '" << name << "'.\n";
      };
  auto print_msg = [&](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.final_time = 0.004s;
  run_options.output_interval = run_options.final_time / (60 * 4);
  run_options.cfl = 0.8;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}
