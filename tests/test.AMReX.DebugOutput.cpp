// Copyright (c) 2018 Maikel Nadolski
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
#include "fub/Equation.hpp"
#include "fub/AMReX/output/DebugOutput.hpp"

#include <iostream>

struct CircleData {
  void InitializeData(amrex::MultiFab& data,
                      const amrex::Geometry& geom) const {
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
  fub::amrex::ScopeGuard amrex_scope_guard{};

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

  auto grid = std::make_shared<fub::amrex::GriddingAlgorithm>(
      hierarchy, CircleData{}, gradient);
  grid->InitializeHierarchy(0.0);

  const int nlevels = grid->GetPatchHierarchy().GetNumberOfLevels();
  std::size_t size = static_cast<std::size_t>(nlevels);
  ::amrex::Vector<const ::amrex::MultiFab*> mf(size);
  ::amrex::Vector<const ::amrex::Geometry*> geoms(size);
  for (std::size_t i = 0; i < size; ++i) {
    mf[i] = &grid->GetPatchHierarchy().GetPatchLevel(static_cast<int>(i)).data;
    geoms[i] = &grid->GetPatchHierarchy().GetGeometry(static_cast<int>(i));
  }

  ::amrex::Vector<const ::amrex::MultiFab*> mf2(2);
  ::amrex::Vector<const ::amrex::Geometry*> geoms2(2);
  mf2[0] = &grid->GetPatchHierarchy().GetPatchLevel(1).data;
  mf2[1] = &grid->GetPatchHierarchy().GetPatchLevel(3).data;
  geoms2[0] = &grid->GetPatchHierarchy().GetGeometry(1);
  geoms2[1] = &grid->GetPatchHierarchy().GetGeometry(3);

  ::amrex::BoxArray ba3 = mf[3]->boxArray();
  ba3.surroundingNodes();
  ::amrex::DistributionMapping dm3 = mf[3]->DistributionMap();
  ::amrex::Geometry geom3 = *geoms[3];
  ::amrex::MultiFab mf3(ba3, dm3, 1, 0);
  mf3.setVal(1.0);

  using Traits = fub::StateTraits<fub::Complete<fub::Advection2d>>;
  constexpr auto names = Traits::names;

  fub::amrex::DebugStorage& debug = *grid->GetPatchHierarchy().GetDebugStorage();
  debug.Enable();
  fub::amrex::DebugSnapshotProxy dbg_snapshot = debug.AddSnapshot("test");

  dbg_snapshot.SaveData(mf, std::get<0>(names), geoms);
  dbg_snapshot.SaveData(mf2, "Mass_sub", geoms2);
  dbg_snapshot.SaveData(mf3, "const_val", geom3);

  debug.FlushData("DebugTest");

  dbg_snapshot = debug.AddSnapshot("test");
  dbg_snapshot.SaveData(mf, std::get<0>(names), geoms);

  debug.FlushData("DebugTest");

}
