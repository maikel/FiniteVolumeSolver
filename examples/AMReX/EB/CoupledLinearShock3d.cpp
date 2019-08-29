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
#include <boost/program_options.hpp>

#include <cmath>
#include <iostream>

#include <xmmintrin.h>

void WriteMatlabData(std::ostream& out, const amrex::FArrayBox& fab,
                     const fub::IdealGasMix<3>& eq,
                     const amrex::Geometry& geom) {
  using namespace fub;
  auto view = fub::amrex::MakeView<const Complete<IdealGasMix<3>>>(fab, eq);
  out << fmt::format(
      "X Y Density VelocityX VelocityY VelocityZ Temperature Pressure\n");
  ForEachIndex(
      Box<0>(view), [&](std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k) {
        double x[3] = {0.0, 0.0, 0.0};
        ::amrex::IntVect iv{int(i), int(j), int(k)};
        geom.CellCenter(iv, x);
        const double density = view.density(i, j, k);
        const double velocity_x =
            density > 0.0 ? view.momentum(i, j, k, 0) / density : 0.0;
        const double velocity_y =
            density > 0.0 ? view.momentum(i, j, k, 1) / density : 0.0;
        const double velocity_z =
            density > 0.0 ? view.momentum(i, j, k, 2) / density : 0.0;
        const double temperature = view.temperature(i, j, k);
        const double pressure = view.pressure(i, j, k);
        out << fmt::format("{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e} "
                           "{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}\n",
                           x[0], x[1], density, velocity_x, velocity_y,
                           velocity_z, temperature, pressure);
      });
}

namespace fub::amrex::cutcell {

void Write2Dfrom3D(std::ostream& out, const PatchHierarchy& hierarchy,
                   const IdealGasMix<3>& eq, fub::Duration time_point,
                   std::ptrdiff_t cycle_number, MPI_Comm comm) {
  const std::size_t n_level =
      static_cast<std::size_t>(hierarchy.GetNumberOfLevels());
  std::vector<::amrex::Geometry> geoms{};
  geoms.reserve(n_level);
  std::vector<::amrex::MultiFab> data{};
  std::vector<::amrex::FArrayBox> fabs{};
  data.reserve(n_level);
  for (std::size_t level = 0; level < n_level; ++level) {
    const int ilvl = static_cast<int>(level);
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(ilvl);
    ::amrex::Box domain = level_geom.Domain();
    const int k =
        (level_geom.Domain().smallEnd(2) + level_geom.Domain().bigEnd(2)) / 2;
    domain.setSmall(2, k);
    domain.setBig(2, k);
    ::amrex::RealBox probDomain = level_geom.ProbDomain();
    //    const double dz = level_geom.CellSize(2);
    probDomain.setLo(2, 0.0);
    probDomain.setHi(2, 0.0);
    ::amrex::Geometry& geom = geoms.emplace_back(domain, &probDomain);
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(ilvl).data;
    ::amrex::FArrayBox local_fab(domain, level_data.nComp());
    local_fab.setVal(0.0);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& patch_data = level_data[mfi];
      const ::amrex::Box box = mfi.tilebox() & domain;
      local_fab.copy(patch_data, box);
    });
    int rank = -1;
    ::MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      ::amrex::FArrayBox& fab = fabs.emplace_back(domain, level_data.nComp());
      if (level > 0) {
        for (int comp = 0; comp < level_data.nComp(); ++comp) {
          for (int j = domain.smallEnd(1); j < domain.bigEnd(1); ++j) {
            for (int i = domain.smallEnd(0); i < domain.bigEnd(0); ++i) {
              ::amrex::IntVect fine_i{i, j, domain.smallEnd(2)};
              ::amrex::IntVect coarse_i = fine_i;
              coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
              if (fab(fine_i, 0) == 0.0) {
                fab(fine_i, comp) = fabs[level - 1](coarse_i, comp);
              }
            }
          }
        }
      }
      ::MPI_Reduce(local_fab.dataPtr(), fab.dataPtr(), local_fab.size(),
                   MPI_DOUBLE, MPI_SUM, 0, comm);
      if (level == n_level - 1) {
        out << fmt::format("nx = {}\n", domain.length(0));
        out << fmt::format("ny = {}\n", domain.length(1));
        out << fmt::format("t = {}\n", time_point.count());
        out << fmt::format("cycle = {}\n", cycle_number);
        WriteMatlabData(out, fab, eq, geom);
        out.flush();
      }
    } else {
      ::MPI_Reduce(local_fab.dataPtr(), nullptr, local_fab.size(), MPI_DOUBLE,
                   MPI_SUM, 0, comm);
    }
  }
}

std::vector<double>
GatherStates(const PatchHierarchy& hierarchy,
             basic_mdspan<const double, extents<3, dynamic_extent>> xs,
             MPI_Comm comm) {
  const int nlevel = hierarchy.GetNumberOfLevels();
  const int finest_level = nlevel - 1;
  const int ncomp = hierarchy.GetDataDescription().n_state_components;
  std::vector<double> buffer(xs.extent(1) * ncomp * nlevel);
  mdspan<double, 3> states(buffer.data(), xs.extent(1), ncomp, nlevel);
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(level).data;
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(level);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      ForEachIndex(mfi.tilebox(), [&](auto... is) {
        double lo[3]{};
        double hi[3]{};
        const ::amrex::IntVect iv{int(is)...};
        level_geom.LoNode(iv, lo);
        level_geom.HiNode(iv, hi);
        for (int k = 0; k < xs.extent(1); ++k) {
          if (lo[0] <= xs(0, k) && xs(0, k) < hi[0] && lo[1] <= xs(1, k) &&
              xs(1, k) < hi[1] && lo[2] <= xs(2, k) && xs(2, k) < hi[2]) {
            for (int comp = 0; comp < level_data.nComp(); ++comp) {
              states(k, comp, level) = level_data[mfi](iv, comp);
            }
          }
        }
      });
    });
  }
  std::vector<double> global_buffer(buffer.size());
  ::MPI_Allreduce(buffer.data(), global_buffer.data(), global_buffer.size(),
                  MPI_DOUBLE, MPI_SUM, comm);
  states = mdspan<double, 3>(global_buffer.data(), xs.extent(1), ncomp, nlevel);
  for (int level = 1; level < nlevel; ++level) {
    for (int comp = 0; comp < ncomp; ++comp) {
      for (int i = 0; i < xs.extent(1); ++i) {
        if (states(i, comp, level) == 0.0) {
          states(i, comp, level) = states(i, comp, level - 1);
        }
      }
    }
  }
  std::vector<double> result(&states(0, 0, finest_level),
                             &states(0, 0, finest_level) +
                                 xs.extent(1) * ncomp);
  return result;
}

} // namespace fub::amrex::cutcell

namespace fub::amrex {
std::vector<double>
GatherStates(const PatchHierarchy& hierarchy,
             basic_mdspan<const double, extents<3, dynamic_extent>> xs,
             MPI_Comm comm) {
  const int nlevel = hierarchy.GetNumberOfLevels();
  const int finest_level = nlevel - 1;
  const int ncomp = hierarchy.GetDataDescription().n_state_components;
  std::vector<double> buffer(xs.extent(1) * ncomp * nlevel);
  mdspan<double, 3> states(buffer.data(), xs.extent(1), ncomp, nlevel);
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(level).data;
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(level);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      ForEachIndex(mfi.tilebox(), [&](auto... is) {
        double lo[3]{};
        double hi[3]{};
        const ::amrex::IntVect iv{int(is)...};
        level_geom.LoNode(iv, lo);
        level_geom.HiNode(iv, hi);
        for (int k = 0; k < xs.extent(1); ++k) {
          if (lo[0] <= xs(0, k) && xs(0, k) < hi[0] && lo[1] <= xs(1, k) &&
              xs(1, k) < hi[1] && lo[2] <= xs(2, k) && xs(2, k) < hi[2]) {
            for (int comp = 0; comp < level_data.nComp(); ++comp) {
              states(k, comp, level) = level_data[mfi](iv, comp);
            }
          }
        }
      });
    });
  }
  std::vector<double> global_buffer(buffer.size());
  ::MPI_Allreduce(buffer.data(), global_buffer.data(), global_buffer.size(),
                  MPI_DOUBLE, MPI_SUM, comm);
  states = mdspan<double, 3>(global_buffer.data(), xs.extent(1), ncomp, nlevel);
  for (int level = 1; level < nlevel; ++level) {
    for (int comp = 0; comp < ncomp; ++comp) {
      for (int i = 0; i < xs.extent(1); ++i) {
        if (states(i, comp, level) == 0.0) {
          states(i, comp, level) = states(i, comp, level - 1);
        }
      }
    }
  }
  std::vector<double> result(&states(0, 0, finest_level),
                             &states(0, 0, finest_level) +
                                 xs.extent(1) * ncomp);
  return result;
}

} // namespace fub::amrex

static constexpr int Tube_Rank = 1;
static constexpr int Plenum_Rank = 3;

struct TemperatureRamp {
  using Complete = fub::IdealGasMix<Tube_Rank>::Complete;

  fub::IdealGasMix<Tube_Rank> equation_;
  double ignition_pos_;

  void InitializeData(::amrex::MultiFab& data, const ::amrex::Geometry& geom) {
    fub::FlameMasterReactor& reactor = equation_.GetReactor();
    reactor.SetMoleFractions("N2:79,O2:21,H2:42");
    const double high_temp = 1450.0;
    const double low_temp = 300.0;
    Complete complete(equation_);
    Complete cold(equation_);
    reactor.SetMoleFractions("N2:79,O2:21,H2:42");
    reactor.SetTemperature(low_temp);
    reactor.SetPressure(101325.0);
    equation_.CompleteFromReactor(cold);

    fub::amrex::ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
      fub::View<Complete> states =
          fub::amrex::MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
      fub::ForEachIndex(fub::Box<0>(states), [&](std::ptrdiff_t i) {
        double x[AMREX_SPACEDIM] = {};
        geom.CellCenter(::amrex::IntVect{int(i)}, x);
        const double x0 = ignition_pos_ - x[0];
        if (0.0 < x0 && x0 < 0.05) {
          const double d = std::clamp(x0 / 0.05, 0.0, 1.0);
          reactor.SetMoleFractions("N2:79,O2:21,H2:42");
          reactor.SetTemperature(d * high_temp + (1.0 - d) * low_temp);
          reactor.SetPressure(101325.0);
          equation_.CompleteFromReactor(complete);
          fub::Store(states, complete, {i});
        } else if (-0.05 < x0 && x0 < 0.0) {
          const double d = std::clamp(-x0 / 0.05, 0.0, 1.0);
          reactor.SetMoleFractions("N2:79,O2:21,H2:42");
          reactor.SetTemperature(d * high_temp + (1.0 - d) * low_temp);
          reactor.SetPressure(101325.0);
          equation_.CompleteFromReactor(complete);
          fub::Store(states, complete, {i});
        } else {
          fub::Store(states, cold, {i});
        }
      });
    });
  }
};

auto MakeTubeSolver(int num_cells, double length, double ignition_pos,
                    int n_level, fub::Burke2012& mechanism) {
  const std::array<int, AMREX_SPACEDIM> n_cells{num_cells, 1, 1};
  const double tube_radius = 0.015;
  const std::array<double, AMREX_SPACEDIM> xlower{-length - 0.03, -tube_radius,
                                                  -tube_radius};
  const std::array<double, AMREX_SPACEDIM> xupper{-0.03, tube_radius,
                                                  tube_radius};
  const std::array<int, AMREX_SPACEDIM> periodicity{0, 0, 0};

  amrex::RealBox xbox(xlower, xupper);
  fub::IdealGasMix<Tube_Rank> equation{fub::FlameMasterReactor(mechanism)};

  using namespace fub::amrex;

  CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);

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
  GradientDetector gradient{equation, std::make_pair(&Complete::density, 5e-2),
                            std::make_pair(&Complete::pressure, 5e-2)};

  ::amrex::Box refine_box{{num_cells - 5, 0, 0}, {num_cells - 1, 0, 0}};
  ConstantBox constant_box{refine_box};

  TemperatureRamp initial_data{equation, ignition_pos};

  BoundarySet boundaries{
      {ReflectiveBoundary{fub::execution::seq, equation, fub::Direction::X, 0},
       TransmissiveBoundary{fub::Direction::X, 1}}};

  //  PatchHierarchy hierarchy =
  //  ReadCheckpointFile("/Volumes/Maikel_Intenso/FiniteVolumeSolver_Build_3d/LongLinearShock_3d/Checkpoint/Tube_00108",
  //  MakeDataDescription(equation), geometry, hier_opts);
  PatchHierarchy hierarchy(equation, geometry, hier_opts);
  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      hierarchy, initial_data, TagAllOf(gradient, constant_box, TagBuffer(2)),
      boundaries);
  gridding->InitializeHierarchy(0.0);

  //  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Tube_Rank>> signals{};
  //  fub::HllMethod hll_method(equation, signals);
  // fub::MusclHancockMethod flux_method{equation, hll_method};
  fub::ideal_gas::MusclHancockPrimMethod<Tube_Rank> flux_method(equation);

  HyperbolicMethod method{FluxMethod(fub::execution::seq, flux_method),
                          ForwardIntegrator(fub::execution::seq),
                          Reconstruction(fub::execution::seq, equation)};

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

auto MakePlenumSolver(int num_cells, double length, int n_level,
                      fub::Burke2012& mechanism) {
  const std::array<int, Plenum_Rank> n_cells{num_cells, num_cells, num_cells};
  const double half = 0.5 * length;

  const std::array<double, Plenum_Rank> xlower{-0.03, -half, -half};
  const std::array<double, Plenum_Rank> xupper{length - 0.03, +0.15, +0.15};
  const std::array<int, Plenum_Rank> periodicity{0, 0, 0};

  amrex::RealBox xbox(xlower, xupper);
  amrex::Geometry coarse_geom(
      amrex::Box{{}, {n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1}}, &xbox,
      -1, periodicity.data());

  auto embedded_boundary = amrex::EB2::makeIntersection(
      amrex::EB2::PlaneIF({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, false),
      amrex::EB2::CylinderIF(0.015, -1.0, 0, {1e6, 0.0, 0.0}, true));
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
                             std::pair{&State::density, 0.05}};

  ::amrex::RealBox inlet{{-0.1, -0.015, -0.015}, {0.01, +0.015, +0.015}};
  const ::amrex::Box refine_box = BoxWhichContains(inlet, coarse_geom);
  ConstantBox constant_box{refine_box};

  const double p0 = 101325.0;
  BoundarySet boundary_condition{
      {IsentropicPressureBoundary{equation, p0, fub::Direction::X, 0},
       IsentropicPressureBoundary{equation, p0, fub::Direction::X, 1},
       IsentropicPressureBoundary{equation, p0, fub::Direction::Y, 0},
       IsentropicPressureBoundary{equation, p0, fub::Direction::Y, 1},
       IsentropicPressureBoundary{equation, p0, fub::Direction::Z, 0},
       IsentropicPressureBoundary{equation, p0, fub::Direction::Z, 1}}};

  //  PatchHierarchy hierarchy =
  //  ReadCheckpointFile("/Volumes/Maikel_Intenso/FiniteVolumeSolver_Build_3d/LongLinearShock_3d/Checkpoint/Plenum_00108",
  //  fub::amrex::MakeDataDescription(equation), geometry, options);
  PatchHierarchy hierarchy(equation, geometry, options);
  std::shared_ptr gridding = std::make_shared<GriddingAlgorithm>(
      std::move(hierarchy), initial_data,
      TagAllOf(TagCutCells(), gradients, constant_box, TagBuffer(2)),
      boundary_condition);
  gridding->InitializeHierarchy(0.0);

  // Make Solver

  fub::EinfeldtSignalVelocities<fub::IdealGasMix<Plenum_Rank>> signals{};
  fub::HllMethod hll_method{equation, signals};
  //  fub::MusclHancockMethod flux_method(equation, hll_method);
  fub::ideal_gas::MusclHancockPrimMethod<Plenum_Rank> flux_method(equation);
  fub::KbnCutCellMethod cutcell_method(flux_method, hll_method);

  HyperbolicMethod method{FluxMethod{fub::execution::simd, cutcell_method},
                          fub::amrex::cutcell::TimeIntegrator{},
                          Reconstruction{fub::execution::simd, equation}};

  return fub::amrex::cutcell::IntegratorContext(gridding, method);
}

struct ProgramOptions {
  double final_time{0.20};
  double cfl{0.8};
  double plenum_domain_length{1.0};
  int plenum_n_cells{128};
  int max_refinement_level{1};
  double tube_domain_length{1.47};
  double tube_ignition_position{0.6 - 1.47};
  double output_interval{1.0E-4};
};

void MyMain(const ProgramOptions& po);

std::optional<ProgramOptions> ParseCommandLine(int argc, char** argv) {
  namespace po = boost::program_options;
  ProgramOptions opts;
  po::options_description desc("Program Options");
  desc.add_options()("help", "Print help messages")(
      "cfl", po::value<double>(&opts.cfl)->default_value(opts.cfl),
      "Set the CFL condition")("max_refinement_level,r",
                               po::value<int>(&opts.max_refinement_level)
                                   ->default_value(opts.plenum_domain_length),
                               "Set the maximal refinement level")(
      "ncells,n",
      po::value<int>(&opts.plenum_n_cells)->default_value(opts.plenum_n_cells),
      "Set number of cells in each direction for the plenum")(
      "final_time,t",
      po::value<double>(&opts.final_time)->default_value(opts.final_time),
      "Set the final simulation time")(
      "plenum_len,p",
      po::value<double>(&opts.plenum_domain_length)
          ->default_value(opts.plenum_domain_length),
      "Set the base length for the plenum")(
      "tube_len,l",
      po::value<double>(&opts.tube_domain_length)
          ->default_value(opts.tube_domain_length),
      "Set the base length for the tube")(
      "tube_ignition,i",
      po::value<double>(&opts.tube_ignition_position)
          ->default_value(opts.tube_ignition_position),
      "Set the position for the ignition of the detonation inside the tube")(
      "output_interval,o",
      po::value<double>(&opts.output_interval)
          ->default_value(opts.output_interval),
      "Sets the output interval");
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
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

  amrex::Print() << fmt::format(
      "[Info] Simulation will run with the following options:\n[Info]\n");
  amrex::Print() << fmt::format("[Info] final_time = {}s\n", opts.final_time);
  amrex::Print() << fmt::format("[Info] output_interval = {}s\n",
                                opts.output_interval);
  amrex::Print() << fmt::format("[Info] cfl = {}\n", opts.cfl);
  amrex::Print() << fmt::format("[Info] max_refinement_level = {}\n",
                                opts.max_refinement_level);
  std::array<double, 3> plo{-0.03, -0.5 * opts.plenum_domain_length,
                            -0.5 * opts.plenum_domain_length};
  std::array<double, 3> pup{opts.plenum_domain_length - 0.03,
                            +0.5 * opts.plenum_domain_length,
                            +0.5 * opts.plenum_domain_length};
  amrex::Print() << fmt::format(
      "[Info] plenum_domain = [{}, {}] x [{}, {}] x [{}, {}]\n", plo[0], pup[0],
      plo[1], pup[1], plo[2], pup[2]);
  amrex::Print() << fmt::format(
      "[Info] plenum_n_cells = {}\n[Info] plenum_dx = {}\n",
      opts.plenum_n_cells, opts.plenum_domain_length / opts.plenum_n_cells);
  std::array<double, 3> tlo{-opts.tube_domain_length - 0.03, -0.015, -0.015};
  std::array<double, 3> tup{-0.03, +0.015, +0.015};
  amrex::Print() << fmt::format(
      "[Info] tube_domain = [{}, {}] x [{}, {}] x [{}, {}]\n", tlo[0], tup[0],
      tlo[1], tup[1], tlo[2], tup[2]);
  const double ratio = opts.tube_domain_length / opts.plenum_domain_length;
  const int tube_n_cells =
      static_cast<int>(ratio * static_cast<double>(opts.plenum_n_cells));
  amrex::Print() << fmt::format(
      "[Info] tube_n_cells = {}\n[Info] tube_dx = {}\n", tube_n_cells,
      opts.tube_domain_length / tube_n_cells);

  return opts;
}

int main(int argc, char** argv) {
  fub::amrex::ScopeGuard _{};
  std::optional<ProgramOptions> po = ParseCommandLine(argc, argv);
  if (po) {
    MyMain(*po);
  }
}

void MyMain(const ProgramOptions& po) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::Burke2012 mechanism{};

  const int n_level = po.max_refinement_level;
  const int n_plenum_cells = po.plenum_n_cells;
  const int n_tube_cells = [&] {
    const double len_ratio = po.tube_domain_length / po.plenum_domain_length;
    int cells =
        static_cast<int>(len_ratio * static_cast<double>(n_plenum_cells));
    cells -= cells % 8;
    return cells;
  }();
  auto plenum = MakePlenumSolver(n_plenum_cells, po.plenum_domain_length,
                                 n_level, mechanism);
  auto tube = MakeTubeSolver(n_tube_cells, po.tube_domain_length,
                             po.tube_ignition_position, n_level, mechanism);

  ::amrex::RealBox inlet{{-0.1, -0.015, -0.015}, {0.05, +0.015, +0.015}};

  fub::amrex::BlockConnection connection;
  connection.direction = fub::Direction::X;
  connection.side = 0;
  connection.plenum.id = 0;
  connection.plenum.mirror_box = BoxWhichContains(inlet, plenum.GetGeometry(0));
  connection.tube.id = 0;
  connection.tube.mirror_box =
      tube.GetGriddingAlgorithm()->GetPatchHierarchy().GetGeometry(0).Domain();

  fub::IdealGasMix<Plenum_Rank> equation{mechanism};
  fub::IdealGasMix<Tube_Rank> tube_equation{mechanism};

  fub::amrex::MultiBlockIntegratorContext context(
      fub::FlameMasterReactor(mechanism), {std::move(tube)},
      {std::move(plenum)}, {connection});

  fub::DimensionalSplitLevelIntegrator system_solver(fub::int_c<Plenum_Rank>,
                                                     std::move(context));

  fub::amrex::MultiBlockKineticSouceTerm source_term{
      tube_equation, system_solver.GetGriddingAlgorithm()};

  fub::DimensionalSplitSystemSourceSolver solver{system_solver, source_term};

  std::string base_name = "LongLinearShock_3d";

  std::vector<double> probes_buffer(6 * 3);
  fub::basic_mdspan<double, fub::extents<3, fub::dynamic_extent>> probes(
      probes_buffer.data(), 6);
  probes(0, 0) = 1.0E-6 * 0.03;
  probes(0, 1) = 4.0 * 0.03;
  probes(0, 2) = 9.0 * 0.03;
  probes(0, 3) = 14.0 * 0.03;
  probes(0, 4) = 19.0 * 0.03;
  probes(0, 5) = 24.0 * 0.03;

  std::vector<double> tube_probes_buffer(5 * 3);
  fub::basic_mdspan<double, fub::extents<3, fub::dynamic_extent>> tube_probes(
      tube_probes_buffer.data(), 5);
  tube_probes(0, 0) = -4.0 * 0.03;
  tube_probes(0, 1) = -9.0 * 0.03;
  tube_probes(0, 2) = -14.0 * 0.03;
  tube_probes(0, 3) = -19.0 * 0.03;
  tube_probes(0, 4) = -24.0 * 0.03;

  // Write Checkpoints 0min + every 5min
  //  const fub::Duration checkpoint_offest = std::chrono::minutes(5);
  //  fub::Duration next_checkpoint = std::chrono::minutes(0);
  auto output =
      [&](std::shared_ptr<fub::amrex::MultiBlockGriddingAlgorithm> gridding,
          auto cycle, auto timepoint, int output_num) {
        if (output_num == 0) {
          ::amrex::Print() << "Start Checkpointing.\n";
          fub::amrex::WriteCheckpointFile(
              fmt::format("{}/Checkpoint/Tube_{:05}", base_name, cycle),
              gridding->GetTubes()[0]->GetPatchHierarchy());
          fub::amrex::cutcell::WriteCheckpointFile(
              fmt::format("{}/Checkpoint/Plenum_{:05}", base_name, cycle),
              gridding->GetPlena()[0]->GetPatchHierarchy());
          ::amrex::Print() << "Finish Checkpointing.\n";
        }
        if (output_num >= 0) {
          std::string name = fmt::format("{}/Tube/plt{:05}", base_name, cycle);
          ::amrex::Print() << "Start output to '" << name << "'.\n";
          fub::amrex::WritePlotFile(
              name, gridding->GetTubes()[0]->GetPatchHierarchy(),
              tube_equation);
          ::amrex::Print() << "Finished output to '" << name << "'.\n";
          name = fmt::format("{}/Plenum/plt{:05}", base_name, cycle);
          ::amrex::Print() << "Start output to '" << name << "'.\n";
          fub::amrex::cutcell::WritePlotFile(
              name, gridding->GetPlena()[0]->GetPatchHierarchy(), equation);
          ::amrex::Print() << "Finished output to '" << name << "'.\n";

          ::amrex::Print() << "Begin Matlab Output.\n";
          std::ofstream out(
              fmt::format("{}/Plenum_{:05}.dat", base_name, cycle),
              std::ios::trunc);
          fub::amrex::cutcell::Write2Dfrom3D(
              out, gridding->GetPlena()[0]->GetPatchHierarchy(), equation,
              timepoint, cycle, context.GetMpiCommunicator());
          ::amrex::Print() << "End Matlab Output.\n";

          ::amrex::Print() << "Start Output for Probes.\n";
          {
            std::vector<double> buffer =
                GatherStates(gridding->GetTubes()[0]->GetPatchHierarchy(),
                             tube_probes, context.GetMpiCommunicator());
            fub::mdspan<const double, 2> states(
                buffer.data(), tube_probes.extent(1),
                buffer.size() / tube_probes.extent(1));
            for (int i = tube_probes.extent(1) - 1; i >= 0; --i) {
              const double rho = states(i, 0);
              const double u = states(i, 1) / rho;
              const double T = states(i, 16);
              const double p = states(i, 14);
              const double t = timepoint.count();
              const double x = tube_probes(0, i);
              ::amrex::Print()
                  << fmt::format("{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< "
                                 "24.15g}{:< 24.15g}{:< 24.15g}\n",
                                 x, t, rho, u, T, p);
            }
          }

          {
            std::vector<double> buffer =
                GatherStates(gridding->GetPlena()[0]->GetPatchHierarchy(),
                             probes, context.GetMpiCommunicator());
            fub::mdspan<const double, 2> states(buffer.data(), probes.extent(1),
                                                buffer.size() /
                                                    probes.extent(1));
            for (int i = 0; i < probes.extent(1); ++i) {
              const double rho = states(i, 0);
              const double u = states(i, 1) / rho;
              const double v = states(i, 2) / rho;
              const double w = states(i, 3) / rho;
              const double a = states(i, 17);
              const double Ma = std::sqrt(u * u + v * v + w * w) / a;
              const double T = states(i, 18);
              const double p = states(i, 16);
              const double t = timepoint.count();
              const double x = probes(0, i);
              ::amrex::Print()
                  << fmt::format("{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< "
                                 "24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}\n",
                                 x, t, rho, Ma, a, T, p);
            }
          }
          ::amrex::Print() << "End Output for Probes.\n";
        }
      };

  using namespace std::literals::chrono_literals;
  output(solver.GetGriddingAlgorithm(), solver.GetCycles(),
         solver.GetTimePoint(), 0);
  fub::RunOptions run_options{};
  run_options.final_time = fub::Duration(po.final_time);
  run_options.output_interval =
      std::vector<fub::Duration>{fub::Duration(po.output_interval), 0.0s};
  run_options.output_frequency = std::vector<int>{0, 1};
  run_options.cfl = po.cfl;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     fub::amrex::print);
}
