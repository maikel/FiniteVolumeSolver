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

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB_LSCore.H>

#include <array>
#include <iostream>
#include <memory>

auto Tube() {
  return amrex::EB2::makeIntersection(
      amrex::EB2::PlaneIF({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, false),
      amrex::EB2::CylinderIF(0.015, -1.0, 0, {1e6, 0.0, 0.0}, true));
}

amrex::Geometry MakeGeometry(const std::array<int, 3>& n_cells,
                             const std::array<double, 3>& lower,
                             const std::array<double, 3>& upper,
                             const std::array<int, 3>& periodicity) {
  amrex::Box box({0, 0, 0}, {n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1});
  amrex::RealBox realbox{lower, upper};
  return amrex::Geometry(box, realbox, -1, periodicity);
}

void MyMain(int cells) {
  const std::array<int, 3> n_cells{cells, cells, cells};
  const std::array<double, 3> xlower{-0.10, -0.15, -0.15};
  const std::array<double, 3> xupper{+0.20, +0.15, +0.15};
  const std::array<int, 3> periodicity{0, 0, 0};
  const amrex::Geometry geom =
      MakeGeometry(n_cells, xlower, xupper, periodicity);
  auto embedded_boundary = Tube();
  auto shop = amrex::EB2::makeShop(embedded_boundary);
  amrex::EB2::Build(shop, geom, 0, 0);

  amrex::Box box({0, 0, 0}, {cells - 1, cells - 1, cells - 1});
  amrex::BoxArray box_array(box);
  amrex::DistributionMapping distribution_mapping(box_array);

  std::unique_ptr<amrex::EBFArrayBoxFactory> factory =
      amrex::makeEBFabFactory(geom, box_array, distribution_mapping, {4, 4, 4},
                              ::amrex::EBSupport::full);

  const auto& flags = factory->getMultiEBCellFlagFab();
  const amrex::MultiFab& volumes = factory->getVolFrac();
  for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
    const amrex::EBCellFlagFab& flag_fab = flags[mfi];
    const amrex::FArrayBox& vol_fab = volumes[mfi];
    const amrex::Box& box = mfi.validbox();
    amrex::IntVect len = box.length();
    for (int k = 0; k < len[2]; ++k) {
      for (int j = 0; j < len[1]; ++j) {
        for (int i = 0; i < len[0]; ++i) {
          amrex::EBCellFlag flag = flag_fab({i, j, k});
          double volume = vol_fab({i, j, k});
          if (flag.isSingleValued() && !(volume > 0.0)) {
            amrex::Print() << i << ", " << j << ", " << k << '\n';
          }
        }
      }
    }
  }
}

int main(int argc, char** argv) {
  ::amrex::Initialize(argc, argv);
  MyMain(80);
  ::amrex::Finalize();
}
