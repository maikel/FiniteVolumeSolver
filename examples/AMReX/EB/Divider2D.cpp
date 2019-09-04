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

#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>

#include <iostream>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <boost/container/pmr/vector.hpp>
#include <boost/filesystem.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

fub::Polygon ReadPolygonData(std::istream& input) {
  std::string line{};
  namespace pmr = boost::container::pmr;
  pmr::vector<double> xs{};
  pmr::vector<double> ys{};
  while (std::getline(input, line)) {
    double x{}, y{};
    std::istringstream linestream(line);
    linestream >> x >> y;
    if (linestream) {
      xs.push_back(x);
      ys.push_back(y);
    }
  }
  const double x0 = xs.front();
  const double x1 = xs.back();
  const double y0 = ys.front();
  const double y1 = ys.back();
  if (x0 != x1 || y0 != y1) {
    throw std::invalid_argument{
        "Invalid Input File: First and last entries are not the same point."};
  }
  return fub::Polygon(std::move(xs), std::move(ys));
}

void WriteMatlabData(std::ostream& out, const amrex::MultiFab& data,
                     const amrex::Geometry& geom, fub::Duration time_point,
                     std::ptrdiff_t cycle_number) {
  amrex::Box box = geom.Domain();
  amrex::BoxArray ba{box};
  amrex::DistributionMapping dm{ba, 1};
  amrex::MultiFab local_copy(ba, dm, data.nComp(), 0, amrex::MFInfo(),
                             data.Factory());
  local_copy.ParallelCopy(data, 0, 0, data.nComp());

  if (amrex::ParallelDescriptor::MyProc() == 0) {
    out << fmt::format("nx = {}\n", box.length(0));
    out << fmt::format("ny = {}\n", box.length(1));
    out << fmt::format("t = {}\n", time_point.count());
    out << fmt::format("cycle = {}\n", cycle_number);
    out << fmt::format("X Y Density VelocityX VelocityY Pressure\n");
    const amrex::FArrayBox& fab = local_copy[0];
    fub::ForEachIndex(
        fub::amrex::AsIndexBox<2>(fab.box()),
        [&](std::ptrdiff_t i, std::ptrdiff_t j) {
          double x[2] = {0.0, 0.0};
          amrex::IntVect iv{int(i), int(j)};
          geom.CellCenter(iv, x);
          const double density = fab(iv, 0) > 0.0 ? fab(iv, 0) : 0.0;
          const double velocity_x = density > 0.0 ? fab(iv, 1) / density : 0.0;
          const double velocity_y = density > 0.0 ? fab(iv, 2) / density : 0.0;
          const double pressure = density > 0.0 ? fab(iv, 4) : 0.0;
          out << fmt::format("{} {} {} {} {} {}\n", x[0], x[1], density,
                             velocity_x, velocity_y, pressure);
        });
    out.flush();
  }
}

void WriteCheckpoint(const std::string& base_name,
                     const fub::amrex::cutcell::PatchHierarchy& hierarchy,
                     std::ptrdiff_t cycle) {
  std::string path = fmt::format("{}/Checkpoint_{:05}", base_name, cycle);
  amrex::Print() << "Write Checkpoint File to '" << path << "'.\n";
  fub::amrex::cutcell::WriteCheckpointFile(path, hierarchy);
}

int main(int, char** argv) {
  static_assert(AMREX_SPACEDIM == 2);

  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  char* my_argv[] = {argv[0]};
  const fub::amrex::ScopeGuard _(1, my_argv);

  const std::array<int, 2> n_cells{8 * 150, 8 * 42};
  const std::array<double, 2> xlower{0.005, -0.016};
  const std::array<double, 2> xupper{0.155, +0.026};
  const std::array<int, 2> periodicity{0, 0};

  fub::PerfectGas<2> equation;

  using namespace fub::amrex::cutcell;
  using State = fub::Complete<fub::PerfectGas<2>>;
  GradientDetector gradients{equation, std::pair{&State::pressure, 0.05},
                             std::pair{&State::density, 0.005}};

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.22;
  cons.momentum << 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> post_shock_state;
  fub::CompleteFromCons(equation, post_shock_state, cons);

  const double shock_mach_number = 5.8;
  const fub::Array<double, 2, 1> normal{1.0, 0.0};

  ShockMachnumber initial_data(equation, fub::Halfspace({+1.0, 0.0, 0.0}, 0.01),
                               post_shock_state, shock_mach_number, normal);

  const State& pre_shock_state = initial_data.GetRiemannProblem().left;
  amrex::Print() << "Post-Shock-State:\n"
                 << "\tdensity: " << post_shock_state.density << " kg / m^3\n"
                 << "\tvelocity: "
                 << equation.Velocity(post_shock_state).transpose()
                 << " m / s\n"
                 << "\tpressure: " << post_shock_state.pressure << " Pa\n";

  amrex::Print() << "Calculated Pre-Shock-State:\n"
                 << "\tdensity: " << pre_shock_state.density << " kg / m^3\n"
                 << "\tvelocity: "
                 << equation.Velocity(pre_shock_state).transpose() << " m / s\n"
                 << "\tpressure: " << pre_shock_state.pressure << " Pa\n";

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1}},
                              &xbox, -1, periodicity.data());

  const int n_level = 1;

  std::ifstream input("wall_1.txt");
  fub::amrex::Geometry lower{ReadPolygonData(input)};
  input = std::ifstream("wall_2.txt");
  fub::amrex::Geometry mid{ReadPolygonData(input)};
  input = std::ifstream("wall_4.txt");
  fub::amrex::Geometry upper{ReadPolygonData(input)};

  auto embedded_boundary = amrex::EB2::makeUnion(lower, mid, upper);
  auto shop = amrex::EB2::makeShop(embedded_boundary);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  PatchHierarchyOptions options{};
  options.max_number_of_levels = n_level;
  options.index_spaces =
      fub::amrex::cutcell::MakeIndexSpaces(shop, coarse_geom, n_level);

  BoundarySet boundary_condition{{TransmissiveBoundary{fub::Direction::X, 0},
                                  TransmissiveBoundary{fub::Direction::X, 1},
                                  TransmissiveBoundary{fub::Direction::Y, 0},
                                  TransmissiveBoundary{fub::Direction::Y, 1}}};

  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      PatchHierarchy(equation, geometry, options), initial_data,
      TagAllOf(TagCutCells(), gradients, TagBuffer(4)), boundary_condition);
  gridding->InitializeHierarchy(0.0);

  fub::EinfeldtSignalVelocities<fub::PerfectGas<2>> signals{};
  fub::HllMethod hll_method{equation, signals};
  fub::MusclHancockMethod flux_method(equation, hll_method);
  fub::KbnCutCellMethod cutcell_method(std::move(flux_method), hll_method);
  HyperbolicMethod method{FluxMethod{fub::execution::seq, cutcell_method},
                          TimeIntegrator{},
                          Reconstruction{fub::execution::seq, equation}};

  fub::DimensionalSplitLevelIntegrator solver(fub::int_c<2>,
      IntegratorContext(gridding, method));

  std::string base_name = "Divider_58/";

  auto output = [&, count =
                        0LL](const std::shared_ptr<GriddingAlgorithm>& gridding,
                             std::ptrdiff_t cycle, fub::Duration) mutable {
    PatchHierarchy& hierarchy = gridding->GetPatchHierarchy();
    WriteCheckpoint(base_name, hierarchy, cycle);
    std::string name = fmt::format("{}{:05}.dat", base_name, count++);
    std::ofstream file(name);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    WriteMatlabData(file, hierarchy.GetPatchLevel(0).data,
                    hierarchy.GetGeometry(0), hierarchy.GetTimePoint(), cycle);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };

  auto print_msg = [&](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint());
  fub::RunOptions run_options{};
  run_options.final_time = 0.0005s;
  run_options.output_interval = 0.125 * 0.0000125s;
  run_options.cfl = 0.25 * 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}
