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

#include <boost/filesystem.hpp>

#include <cmath>
#include <iostream>

#include <xmmintrin.h>

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 3;

static constexpr double r_tube = 0.015;
static constexpr double r_inner = 0.5 * 0.130;
static constexpr double r_outer = 0.5 * 0.385;
static constexpr double r_tube_center = 0.5 * r_inner + 0.5 * r_outer;
static constexpr double alpha = 2. * M_PI / 5.;

auto Center(double x, double phi) -> ::amrex::RealArray {
  using std::cos;
  using std::sin;
  return {x, r_tube_center * sin(phi), r_tube_center * cos(phi)};
};

auto DomainAroundCenter(const ::amrex::RealArray& x, double rx)
    -> ::amrex::RealBox {
  return ::amrex::RealBox{{x[0] - rx, x[1] - r_tube, x[2] - r_tube},
                          {x[0] + rx, x[1] + r_tube, x[2] + r_tube}};
}

struct TemperatureRamp {
  using Complete = fub::IdealGasMix<Tube_Rank>::Complete;

  fub::IdealGasMix<Tube_Rank> equation_;

  void InitializeData(::amrex::MultiFab& data, const ::amrex::Geometry& geom) {
    fub::FlameMasterReactor& reactor = equation_.GetReactor();
    reactor.SetMoleFractions("N2:79,O2:21,H2:42");
    const double high_temp = 1250.0;
    const double low_temp = 300.0;
    Complete complete(equation_);

    fub::amrex::ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
        double x[AMREX_SPACEDIM] = {};
        geom.CellCenter(::amrex::IntVect{int(i)}, x);
        if (x[0] < -1.4) {
          const double d = std::clamp((x[0] + 1.5) / 0.15, 0.0, 1.0);
          reactor.SetMoleFractions("N2:79,O2:21,H2:42");
          reactor.SetTemperature(d * low_temp + (1.0 - d) * high_temp);
          reactor.SetPressure(101325.0);
        } else {
          reactor.SetMoleFractions("N2:79,O2:21,H2:42");
          reactor.SetTemperature(low_temp);
          reactor.SetPressure(101325.0);
        }
        equation_.CompleteFromReactor(complete);
        fub::Store(states, complete, {i});
      });
    });
  }
};

auto MakeTubeSolver(int num_cells, int n_level, double phi,
                    fub::Burke2012& mechanism) {
  const std::array<int, AMREX_SPACEDIM> n_cells{num_cells, 1, 1};
  constexpr double x_len = -0.03 + 1.5;
  constexpr double r_len = 0.5 * x_len;
  constexpr double x_mid = -0.5 * 0.03 - 0.5 * 1.5;
  amrex::RealBox xbox = DomainAroundCenter(Center(x_mid, phi), r_len);
  const std::array<int, AMREX_SPACEDIM> periodicity{0, 0, 0};

  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};

  using namespace fub::amrex;

  CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = xbox;

  DataDescription desc = MakeDataDescription(equation);

  PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = n_level;
  hier_opts.refine_ratio = amrex::IntVect{2, 1, 1};

  amrex::Geometry geom(
      amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1}}, &xbox,
      -1, periodicity.data());
  geom.refine(hier_opts.refine_ratio);
  ::amrex::EB2::Build(::amrex::EB2::makeShop(::amrex::EB2::AllRegularIF()),
                      geom, hier_opts.refine_ratio, 1, 1);

  using Complete = fub::IdealGasMix<1>::Complete;
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 1e-3),
                            std::make_pair(&Complete::pressure, 1e-2),
                            std::make_pair(&Complete::temperature, 1e-1)};

  ::amrex::Box refine_box{{num_cells - 5, 0, 0}, {num_cells - 1, 0, 0}};
  ConstantBox constant_box{refine_box};

  TemperatureRamp initial_data{equation};

  BoundarySet boundaries{
      {ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::X, 0},
       TransmissiveBoundary{fub::Direction::X, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(desc, geometry, hier_opts), initial_data,
      TagAllOf(gradient, constant_box), boundaries);
  gridding->InitializeHierarchy(0.0);

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<1>> signals{};
  fub::HllMethod hll_method(equation, signals);
  // fub::MusclHancockMethod flux_method{equation, hll_method};

  HyperbolicMethod method{FluxMethod(fub::execution::openmp, hll_method),
                          ForwardIntegrator(fub::execution::openmp),
                          Reconstruction(fub::execution::openmp, equation)};

  return fub::amrex::IntegratorContext(gridding, method);
}

::amrex::Box BoxWhichContains(const ::amrex::RealBox& xbox,
                              const ::amrex::Geometry& geom) {
  ::amrex::Box domain = geom.Domain();
  ::amrex::IntVect lo = domain.smallEnd();
  ::amrex::IntVect up = domain.bigEnd();
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    for (int i = domain.smallEnd(d); i < domain.bigEnd(d); ++i) {
      const double x = geom.CellCenter(i, d);
      if (x < xbox.lo(d)) {
        lo[d] = std::max(lo[d], i);
      }
      if (x > xbox.hi(d)) {
        up[d] = std::min(up[d], i);
      }
    }
  }
  return ::amrex::Box{lo, up};
}

auto MakePlenumSolver(int num_cells, int n_level, fub::Burke2012& mechanism) {
  const std::array<int, Plenum_Rank> n_cells{num_cells, num_cells, num_cells};
  const std::array<double, Plenum_Rank> xlower{-0.03, -0.5 * 0.56, -0.5 * 0.56};
  const std::array<double, Plenum_Rank> xupper{+0.53, +0.5 * 0.56, +0.5 * 0.56};
  const std::array<int, Plenum_Rank> periodicity{0, 0, 0};

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(amrex::Box{{}, ::amrex::IntVect(num_cells - 1)},
                              &xbox, -1, periodicity.data());

  auto embedded_boundary = amrex::EB2::makeUnion(
      amrex::EB2::makeIntersection(
          amrex::EB2::CylinderIF(r_outer, 0.5, 0, {0.25, 0.0, 0.0}, true),
          amrex::EB2::CylinderIF(r_tube_center, 1.0, 0,
                                 {1.0 - 1.0e-6, 0.0, 0.0}, true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 0.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 1.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 2.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 3.0 * alpha),
                                 true),
          amrex::EB2::CylinderIF(r_tube, 0.3, 0, Center(-0.1, 4.0 * alpha),
                                 true)),
      amrex::EB2::CylinderIF(r_inner, 1.0, 0, {0.25, 0.0, 0.0}, false));
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

  ::amrex::RealBox inlet{{-0.1, -0.015, -0.015}, {0.05, +0.015, +0.015}};
  const ::amrex::Box refine_box = BoxWhichContains(inlet, coarse_geom);
  ConstantBox constant_box{refine_box};

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1},
                                  TransmissiveBoundary{fub::Direction::Z, 0},
                                  TransmissiveBoundary{fub::Direction::Z, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, options), initial_data,
      TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
      boundary_condition);
  gridding->InitializeHierarchy(0.0);

  // Make Solver

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::MusclHancockMethod flux_method(equation, hll_method);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{fub::execution::openmp, cutcell_method},
                          fub::amrex::cutcell::TimeIntegrator{},
                          Reconstruction{fub::execution::openmp, equation}};

  return fub::amrex::cutcell::IntegratorContext(gridding, method);
}

int main(int /* argc */, char** /* argv */) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard _{};
  fub::Burke2012 mechanism{};

  const int n_level = 2;
  auto plenum = MakePlenumSolver(64, n_level, mechanism);
  auto tube = MakeTubeSolver(200, n_level, 0.0, mechanism);

  //  ::amrex::RealBox inlet{{-0.1, -0.015, -0.015}, {0.05, +0.015, +0.015}};

  fub::amrex::BlockConnection connection;
  connection.direction = fub::Direction::X;
  connection.side = 0;
  connection.plenum.id = 0;
  connection.plenum.mirror_box = BoxWhichContains(
      DomainAroundCenter(Center(0.0, 0.0), 0.03), plenum.GetGeometry(0));
  connection.tube.id = 0;
  connection.tube.mirror_box =
      tube.GetGriddingAlgorithm()->GetPatchHierarchy().GetGeometry(0).Domain();

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};

  fub::amrex::MultiBlockIntegratorContext context(
      fub::FlameMasterReactor(mechanism), {std::move(tube)},
      {std::move(plenum)}, {connection});

  fub::DimensionalSplitLevelIntegrator system_solver(fub::int_c<Plenum_Rank>,
                                                     std::move(context));
  fub::amrex::MultiBlockKineticSouceTerm source_term{
      fub::IdealGasMix<Tube_Rank>{mechanism},
      system_solver.GetGriddingAlgorithm()};

  fub::DimensionalSplitSystemSourceSolver solver{system_solver, source_term};

  std::string base_name = "MultiTube";
  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};

  // Write Checkpoints 25min + every 30min
  const fub::Duration checkpoint_offest = std::chrono::minutes(30);
  fub::Duration next_checkpoint = std::chrono::minutes(25);
  auto output =
      [&](std::shared_ptr<fub::amrex::MultiBlockGriddingAlgorithm> gridding,
          std::ptrdiff_t cycle, fub::Duration, int = 0) {
        std::chrono::steady_clock::time_point now =
            std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<fub::Duration>(
            now - wall_time_reference);
        if (duration > next_checkpoint) {
          fub::amrex::WriteCheckpointFile(
              fmt::format("{}/Checkpoint_Tube/{:05}", base_name, cycle),
              gridding->GetTubes()[0]->GetPatchHierarchy());
          fub::amrex::cutcell::WriteCheckpointFile(
              fmt::format("{}/Checkpoint_Plenum/{:05}", base_name, cycle),
              gridding->GetPlena()[0]->GetPatchHierarchy());
          next_checkpoint += checkpoint_offest;
        }
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

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.final_time = 0.020s;
  run_options.output_interval = {0.1e-3s};
  run_options.cfl = 0.8;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
