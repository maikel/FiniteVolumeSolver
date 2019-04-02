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

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <Eigen/Dense>

#ifdef __APPLE__
#include <xmmintrin.h>
#else
#include <fenv.h>
#endif

Eigen::Vector2d OrthogonalTo(const Eigen::Vector2d& x) {
  return Eigen::Vector2d{x[1], -x[0]};
}

auto Wedge(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
  Eigen::Vector2d p0{0.0, 0.1};
  Eigen::Vector2d norm1 = OrthogonalTo(p1 - p0).normalized();
  Eigen::Vector2d norm2 = OrthogonalTo(p2 - p0).normalized();
  amrex::EB2::PlaneIF plane1({p0[0], p0[1]}, {norm1[0], norm1[1]});
  amrex::EB2::PlaneIF plane2({p0[0], p0[1]}, {norm2[0], norm2[1]}, false);
  return amrex::EB2::makeComplement(
      amrex::EB2::makeIntersection(plane1, plane2));
}

amrex::Geometry MakeGeometry(const std::array<int, 2>& n_cells,
                             const std::array<double, 2>& lower,
                             const std::array<double, 2>& upper,
                             const std::array<int, 2>& periodicity) {
  amrex::Box box(amrex::IntVect{0, 0},
                 amrex::IntVect{n_cells[0] - 1, n_cells[1] - 1});
  amrex::RealBox realbox{lower, upper};
  return amrex::Geometry(box, realbox, -1, periodicity);
}

void MyMain(int cells) {
  const std::array<int, 2> n_cells{cells, cells};
  const std::array<double, 2> xlower{-0.5, 0.0};
  const std::array<double, 2> xupper{+0.5, +1.0};
  const std::array<int, 2> periodicity{0, 0};
  const amrex::Geometry finest_geom =
      MakeGeometry(n_cells, xlower, xupper, periodicity);
  auto embedded_boundary = Wedge({-1.0, +1.0}, {+1.0, +1.0});
  auto shop = amrex::EB2::makeShop(embedded_boundary);
  amrex::EB2::Build(shop, finest_geom, 0, 0);
}

int main(int argc, char** argv) {
  ::amrex::Initialize(argc, argv);
#ifdef __APPLE__
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
#else
  ::feenableexcept(~FE_INVALID);
#endif
  MyMain(64);
  ::amrex::Finalize();
}
