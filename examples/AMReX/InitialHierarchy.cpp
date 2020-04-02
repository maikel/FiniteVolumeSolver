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
  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration /*time*/) const {
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& data = patch_level.data;
    fub::amrex::ForEachFab(data, [&](const amrex::MFIter& mfi) {
      const ::amrex::Box& box = mfi.tilebox();
      amrex::FArrayBox& fab = data[mfi];
      fub::amrex::ForEachIndex(box, [&](int i, int j) {
        const double x = geom.CellCenter(i, 0);
        const double y = geom.CellCenter(j, 1);
        const double norm2 = x * x + y * y;
        constexpr double r2 = 0.25 * 0.25;
        amrex::IntVect iv(i, j);
        if (norm2 < r2) {
          fab(iv, 0) = 3.0;
        } else {
          fab(iv, 0) = 1.0;
        }
      });
    });
  }
};

int main() {
  fub::amrex::ScopeGuard guard{};

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM == 2);

  const std::array<int, Dim> n_cells{128, 128};
  const std::array<double, Dim> xlower{-1.0, -1.0};
  const std::array<double, Dim> xupper{+1.0, +1.0};

  fub::Advection2d equation{{1.0, 1.0}};

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = std::array<int, Dim>{1, 1};

  fub::amrex::PatchHierarchyOptions hier_opts{};
  hier_opts.max_number_of_levels = 4;

  fub::amrex::PatchHierarchy hierarchy(equation, geometry, hier_opts);

  using State = fub::Advection2d::Complete;
  fub::amrex::GradientDetector gradient{equation,
                                        std::pair{&State::mass, 1e-2}};

  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      hierarchy, CircleData{}, gradient);
  gridding->InitializeHierarchy(0.0);

  fub::amrex::WritePlotFile("InitialHierarchy/plt00000",
                            gridding->GetPatchHierarchy(), equation);
}
