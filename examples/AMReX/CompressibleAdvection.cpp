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

#include "fub/Solver.hpp"
#include "fub/AMReX.hpp"

#include "fub/equations/CompressibleAdvection.hpp"

struct InitialData {
  using Complete = fub::CompressibleAdvection<2>::Complete;
  void InitializeData(amrex::MultiFab& mf, const amrex::Geometry& geom) {
    fub::amrex::ForEachFab(mf, [&](const amrex::MFIter& mfi) {
      fub::CompressibleAdvection<2> equation{};
      amrex::FArrayBox& fab = mf[mfi];
      const amrex::Box& box = mfi.tilebox();
      fub::View<Complete> states = fub::amrex::MakeView<Complete>(fab, equation, box);
      fub::ForEachIndex(fub::Box<0>(states), [&](int i, int j) {
        states.density(i, j) = 1.0;
        states.momentum(i, j, 0) = 0.0;
        states.momentum(i, j, 1) = 0.0;
        states.PTdensity(i, j) = 1.0;
        states.velocity(i, j, 0) = 0.0;
        states.velocity(i, j, 1) = 0.0;
        states.PTinverse(i, j) = 1.0;
      });
    });
  }
};

int main()
{
  fub::amrex::ScopeGuard guard{};
  fub::amrex::DataDescription desc{};
  desc.n_state_components = 7;
  desc.n_node_components = 0;
  desc.n_face_components = 1;

  fub::amrex::CartesianGridGeometry geometry{};
  geometry.cell_dimensions = std::array<int, 2>{128, 128};
  geometry.coordinates = amrex::RealBox({-1.0, -1.0},
                                        {+1.0, +1.0});
  geometry.periodicity = std::array<int, 2>{1, 1};

  fub::amrex::PatchHierarchyOptions hier_opts;
  hier_opts.max_number_of_levels = 1;
  hier_opts.refine_ratio = amrex::IntVect{AMREX_D_DECL(2, 2, 1)};

  fub::amrex::PatchHierarchy hierarchy(desc, geometry, hier_opts);

  fub::amrex::GriddingAlgorithm grid(std::move(hierarchy), InitialData{}, fub::amrex::TagBuffer{0});
  grid.InitializeHierarchy(0.0);

  fub::amrex::FluxMethod fm(fub::execution::simd, fub::CompressibleAdvectionFluxMethod<2>());
}