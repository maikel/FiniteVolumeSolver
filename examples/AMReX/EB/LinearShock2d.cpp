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

#include "fub/CartesianCoordinates.hpp"
#include "fub/equations/PerfectGas.hpp"

#include "fub/HyperbolicSplitCutCellPatchIntegrator.hpp"
#include "fub/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/HyperbolicSplitSystemSolver.hpp"
#include "fub/cutcell_method/KbnStabilisation.hpp"
#include "fub/grid/AMReX/FillCutCellData.hpp"
#include "fub/grid/AMReX/GriddingAlgorithm.hpp"
#include "fub/grid/AMReX/ScopeGuard.hpp"
#include "fub/grid/AMReX/cutcell/FluxMethod.hpp"
#include "fub/grid/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/grid/AMReX/cutcell/HyperbolicSplitIntegratorContext.hpp"
#include "fub/grid/AMReX/cutcell/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/grid/AMReX/cutcell/Reconstruction.hpp"
#include "fub/grid/AMReX/cutcell/Tagging.hpp"
#include "fub/tagging/GradientDetector.hpp"
#include "fub/tagging/TagBuffer.hpp"

#include "fub/flux_method/HllMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"

#include "fub/RunSimulation.hpp"

#include "fub/ForEach.hpp"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <iostream>
#include <xmmintrin.h>

constexpr std::ptrdiff_t ipow(int base, int exponent) {
  std::ptrdiff_t prod{1};
  while (exponent > 0) {
    prod *= base;
    exponent -= 1;
  }
  return prod;
}

amrex::Geometry MakeFinestGeometry(const std::array<int, 2>& n_cells,
                                   const std::array<double, 2>& lower,
                                   const std::array<double, 2>& upper,
                                   const std::array<int, 2>& periodicity,
                                   int n_levels) {
  const int n_fine_cells0 =
      static_cast<int>(ipow(2, n_levels - 1) * n_cells[0]);
  const int n_fine_cells1 =
      static_cast<int>(ipow(2, n_levels - 1) * n_cells[1]);
  amrex::Box box(amrex::IntVect{0, 0},
                 amrex::IntVect{n_fine_cells0 - 1, n_fine_cells1 - 1});
  amrex::RealBox realbox{lower, upper};
  return amrex::Geometry(box, realbox, -1, periodicity);
}

struct ShockData {
  using Equation = fub::PerfectGas<2>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  void InitializeData(const fub::View<Complete>& data,
                      fub::amrex::PatchHandle patch) {
    using namespace fub;
    const ::amrex::Geometry& geom = hierarchy->GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->tilebox();
    const auto& factory = hierarchy->GetPatchLevel(patch.level).factory;
    const auto& flags = factory->getMultiEBCellFlagFab()[*patch.iterator];
    ::amrex::FabType type = flags.getType(box);
    if (type == ::amrex::FabType::covered) {
      ForEachVariable<Complete>(
          [](const auto& var) {
            span<double> span = var.Span();
            std::fill(span.begin(), span.end(), 0.0);
          },
          data);
      return;
    }

    CartesianCoordinates x = fub::amrex::GetCartesianCoordinates(geom, box);
    if (type == ::amrex::FabType::regular) {
      ForEachIndex(Box<0>(data), [&](auto... is) {
        if (x(is...)[0] < -0.04) {
          Store(data, left, {is...});
        } else {
          Store(data, right, {is...});
        }
      });
    } else {
      CutCellData<2> eb = hierarchy->GetCutCellData(patch, Direction::X);
      const PatchDataView<const ::amrex::EBCellFlag, 2>& flags = eb.flags;
      ForEachIndex(Box<0>(data), [&](auto... is) {
        if (flags(is...).isCovered()) {
          Store(data, zero, {is...});
        } else if (x(is...)[0] < -0.04) {
          Store(data, left, {is...});
        } else {
          Store(data, right, {is...});
        }
      });
    }
  }

  std::shared_ptr<fub::amrex::cutcell::PatchHierarchy> hierarchy;
  Equation equation;
  Complete left;
  Complete right;
  Complete zero{equation};
};

struct TagCutCells {
  using Equation = fub::PerfectGas<2>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  void TagCellsForRefinement(const fub::PatchDataView<char, 2>& tags,
                             const fub::View<const Complete>& /* states */,
                             const fub::CartesianCoordinates& /* x */) {
    const fub::span<char> span = tags.Span();
    std::fill(span.begin(), span.end(), char(0));
  }

  void TagCellsForRefinement(const fub::PatchDataView<char, 2>& tags,
                             const fub::View<const Complete>& /* states */,
                             const fub::CutCellData<2>& cutcell_data,
                             const fub::CartesianCoordinates& /* x */) {
    const auto& flags = cutcell_data.flags;
    fub::ForEachIndex(Intersect(flags.Box(), tags.Box()), [&](auto... is) {
      const amrex::EBCellFlag& flag = flags(is...);
      if (flag.isSingleValued() && !tags(is...)) {
        tags(is...) = 1;
      }
    });
  }
};

struct TransmissiveBoundary {
  using Equation = fub::PerfectGas<2>;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  void operator()(const fub::PatchDataView<double, 3>& data,
                  fub::amrex::PatchHandle, fub::Location location,
                  int fill_width, fub::Duration) {
    fub::View<Complete> complete =
        fub::amrex::MakeView<fub::View<Complete>>(data, equation);
    const std::size_t dir = static_cast<std::size_t>(location.direction);
    const std::size_t other_dir = (dir + 1) % 2;
    std::array<std::ptrdiff_t, 3> origin = data.Origin();
    if (location.side == 0) {
      for (std::ptrdiff_t j = 0; j < data.Extent(other_dir); ++j) {
        std::array<std::ptrdiff_t, 2> dest_index{origin[0], origin[1]};
        std::array<std::ptrdiff_t, 2> source_index{origin[0], origin[1]};
        dest_index[other_dir] += j;
        source_index[other_dir] += j;
        source_index[dir] += fill_width;
        Load(state, complete, source_index);
        for (int i = 0; i < fill_width; ++i) {
          dest_index[dir] = origin[dir] + i;
          Store(complete, state, dest_index);
        }
      }
    } else if (location.side == 1) {
      const std::ptrdiff_t n = fub::Extents<0>(complete).extent(dir);
      for (std::ptrdiff_t j = 0; j < data.Extent(other_dir); ++j) {
        std::array<std::ptrdiff_t, 2> dest_index{origin[0], origin[1]};
        std::array<std::ptrdiff_t, 2> source_index{origin[0], origin[1]};
        dest_index[other_dir] += j;
        source_index[other_dir] += j;
        source_index[dir] += n - 1 - fill_width;
        Load(state, complete, source_index);
        for (int i = 0; i < fill_width; ++i) {
          dest_index[dir] = origin[dir] + n - fill_width + i;
          Store(complete, state, dest_index);
        }
      }
    }
  }

  std::shared_ptr<fub::amrex::cutcell::PatchHierarchy> hierarchy;
  Equation equation;
  Complete state{equation};
};

auto Rectangle(const std::array<double, 2>& lower,
               const std::array<double, 2>& upper) {
  amrex::EB2::PlaneIF lower_x({lower[0], lower[1]}, {0, +1});
  amrex::EB2::PlaneIF lower_y({lower[0], lower[1]}, {+1, 0});
  amrex::EB2::PlaneIF upper_x({upper[0], upper[1]}, {0, -1});
  amrex::EB2::PlaneIF upper_y({upper[0], upper[1]}, {-1, 0});
  return amrex::EB2::makeIntersection(lower_x, lower_y, upper_x, upper_y);
}

int main(int argc, char** argv) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();
  fub::amrex::ScopeGuard _(argc, argv);

  const std::array<int, 2> n_cells{128, 128};
  const std::array<double, 2> xlower{-0.10, -0.15};
  const std::array<double, 2> xupper{+0.20, +0.15};
  const std::array<int, 2> periodicity{0, 0};

  const int n_level = 3;

  amrex::Geometry finest_geom =
      MakeFinestGeometry(n_cells, xlower, xupper, periodicity, n_level);

  auto embedded_boundary =
      amrex::EB2::makeUnion(Rectangle({-1.0, +0.015}, {0.0, 1.0}),
                            Rectangle({-1.0, -1.0}, {0.0, -0.015}));

  auto shop = amrex::EB2::makeShop(embedded_boundary);
  amrex::EB2::Build(shop, finest_geom, n_level - 1, n_level - 1);

  fub::PerfectGas<2> equation;
  fub::amrex::DataDescription desc = fub::amrex::MakeDataDescription(equation);

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = periodicity;

  fub::amrex::PatchHierarchyOptions options;
  options.max_number_of_levels = n_level;

  auto hierarchy = std::make_shared<fub::amrex::cutcell::PatchHierarchy>(
      desc, geometry, options);

  amrex::Geometry geom = hierarchy->GetGeometry(0);
  amrex::Box domain{{}, {n_cells[0], n_cells[1]}};
  amrex::BoxArray ba{domain};
  amrex::DistributionMapping dm(ba);
  auto eb_factory =
      amrex::makeEBFabFactory(geom, ba, dm, {4, 4, 4}, amrex::EBSupport::full);

  amrex::WriteSingleLevelPlotfile("LinearShock_Geom", eb_factory->getVolFrac(),
                                  {"VolumeFraction"}, geom, 0.0, 0);

  fub::Conservative<fub::PerfectGas<2>> cons;
  cons.density = 1.0;
  cons.momentum << 0.0, 0.0;
  cons.energy = 101325.0 * equation.gamma_minus_1_inv;
  fub::Complete<fub::PerfectGas<2>> right;
  fub::CompleteFromCons(equation, right, cons);

  cons.energy *= 2.5;
  fub::Complete<fub::PerfectGas<2>> left;
  fub::CompleteFromCons(equation, left, cons);

  ShockData initial_data{hierarchy, equation, left, right};
  TagCutCells cutcells{};
  using State = fub::Complete<fub::PerfectGas<2>>;
  fub::GradientDetector gradients{std::pair{&State::pressure, 0.05},
                                  std::pair{&State::density, 0.005}};
  TransmissiveBoundary boundary{hierarchy, equation};

  fub::HyperbolicSplitCutCellPatchIntegrator patch_integrator{equation};
  fub::HllMethod base_method{equation,
                             fub::EinfeldtSignalVelocities<fub::PerfectGas<2>>};
  fub::MusclHancockMethod muscl_method{equation, base_method};
  fub::KbnCutCellMethod cutcell_method(muscl_method);

  auto gridding = std::make_shared<fub::amrex::cutcell::GriddingAlgorithm>(
      hierarchy, fub::amrex::AdaptInitialData(initial_data, equation),
      fub::amrex::cutcell::AdaptTagging(equation, hierarchy, cutcells,
                                        gradients, fub::TagBuffer(4)),
      boundary);

  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  gridding->InitializeHierarchy(0.0);

  fub::amrex::cutcell::WritePlotFile("LinearShock_0000", *hierarchy, equation);

  const int gcw = 2;
  fub::HyperbolicSplitSystemSolver solver(fub::HyperbolicSplitLevelIntegrator(
      fub::amrex::cutcell::HyperbolicSplitIntegratorContext(gridding, gcw),
      fub::amrex::cutcell::HyperbolicSplitPatchIntegrator(patch_integrator),
      fub::amrex::cutcell::FluxMethod(cutcell_method),
      fub::amrex::cutcell::Reconstruction(equation)));

  std::string base_name = "AMReX/LinearShock_";

  auto output = [&](auto& hierarchy, int cycle, fub::Duration) {
    std::string name = fmt::format("{}{:04}", base_name, cycle);
    ::amrex::Print() << "Start output to '" << name << "'.\n";
    fub::amrex::cutcell::WritePlotFile(name, *hierarchy, equation);
    ::amrex::Print() << "Finished output to '" << name << "'.\n";
  };
  auto print_msg = [](const std::string& msg) { ::amrex::Print() << msg; };

  using namespace std::literals::chrono_literals;
  output(hierarchy, 0, 0s);
  fub::RunOptions run_options{};
  run_options.final_time = 0.002s;
  run_options.output_interval = 0.0000125s;
  run_options.cfl = 0.5 * 0.9;
  fub::RunSimulation(solver, run_options, wall_time_reference, output,
                     print_msg);
}