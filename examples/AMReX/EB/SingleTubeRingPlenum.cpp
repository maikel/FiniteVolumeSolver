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

#include <AMReX.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB_LSCore.H>
#include <AMReX_PlotFileUtil.H>

#include <algorithm>
#include <cmath>
#include <vector>
#include <utility>

/// Rotates a specified 2d geometry around the x axis and produce a 3d geometry
template <typename Base> class Extrude : private Base {
public:
  explicit Extrude(const Base& base) : Base(base) {}

  double operator()(std::array<double, 3> x) const noexcept {
    const double x0 = x[0];
    const double y0 = std::sqrt(x[1] * x[1] + x[2] * x[2]);
    return Base::operator()(x0, y0);
  }

  double operator()(double x, double y, double z) const noexcept {
    return this->operator()({x, y, z});
  }
};

/// Given arrays of x and y coordinates of a convext polygon vertices to define
/// a 2D geometry
class Polygon {
public:
  Polygon(std::vector<double> xs, std::vector<double> ys)
      : xs_{std::move(xs)}, ys_{std::move(ys)} {}

  double operator()(double x, double y) const noexcept;

  double operator()(std::array<double, 2> xs) const noexcept {
    return this->operator()(xs[0], xs[1]);
  }

private:
  std::vector<double> xs_;
  std::vector<double> ys_;
};

int main(int argc, char** argv) {
  ::amrex::Initialize(argc, argv);
  {
    // define some constants for lengths in the geometry

    const double r_tube = 0.015;
    const double d_tube = 2 * r_tube;
    const double plenum_length = 0.25;
    const double plenum_outlet_radius{r_tube};
    const double plenum_jump{2 * r_tube};
    const double r_inner = r_tube;
    const double r_tube_center =
        0.5 * (r_inner + r_inner + 2 * plenum_jump + d_tube);

    // define domain lengths and extents

    const double xlength = 0.36;
    const double ylength = 0.3;

    const double xlower = -0.03;
    const double xupper = xlower + xlength;
    const double ylower = -0.5 * ylength;
    const double yupper = ylower + ylength;

    const double ratio = ylength / xlength;

    auto nCellsY = [ratio](int nx) {
      int base = int(nx * ratio);
      int ny = base - base % 8;
      return ny;
    };

    const int nx = 128;

    std::array<int, AMREX_SPACEDIM> n_cells{
        AMREX_D_DECL(nx, nCellsY(nx), nCellsY(nx))};
    std::array<int, AMREX_SPACEDIM> periodicity{};

    amrex::RealBox xbox{{AMREX_D_DECL(xlower, ylower, ylower)},
                        {AMREX_D_DECL(xupper, yupper, yupper)}};

    // generate the EB geometry

    auto MakePlenum2D = [](auto... points) {
      std::vector<double> xs{}, ys{};

      (xs.push_back(std::get<0>(points)), ...);
      (ys.push_back(std::get<1>(points)), ...);

      return Polygon(std::move(xs), std::move(ys));
    };

    // span a 2D object from given polygon vertices
    auto plenum2D = MakePlenum2D(
        std::pair{+0.00, r_inner}, std::pair{+plenum_length, r_inner},
        std::pair{+plenum_length + 0.03,
                  r_inner + plenum_jump + r_tube - plenum_outlet_radius},
        std::pair{+1.00, r_inner + plenum_jump + r_tube - plenum_outlet_radius},
        std::pair{+1.00, r_inner + plenum_jump + r_tube + plenum_outlet_radius},
        std::pair{+plenum_length + 0.03,
                  r_inner + plenum_jump + r_tube + plenum_outlet_radius},
        std::pair{+plenum_length, r_inner + 2 * plenum_jump + d_tube},
        std::pair{+0.00, r_inner + 2 * plenum_jump + d_tube},
        std::pair{+0.00, r_inner});

#if AMREX_SPACEDIM == 3
    auto embedded_boundary = amrex::EB2::makeComplement(Extrude(plenum2D));
#elif AMREX_SPACEDIM == 2
    auto embedded_boundary = amrex::EB2::makeComplement(plenum2D);
#endif
    auto shop = amrex::EB2::makeShop(embedded_boundary);

    amrex::Geometry coarse_geom(
        amrex::Box{
            {}, {AMREX_D_DECL(n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1)}},
        &xbox, -1, periodicity.data());

    // amrex::EB2::Build(shop, coarse_geom, 0, 0);

    amrex::Geometry fine_geom = coarse_geom;
    fine_geom.refine(amrex::IntVect(2));
    amrex::EB2::Build(shop, fine_geom, 1, 1);

    // create some dummy variable to see some output
    amrex::BoxArray ba(coarse_geom.Domain());
    amrex::DistributionMapping dm(ba);
    std::unique_ptr<amrex::EBFArrayBoxFactory> factory =
        amrex::makeEBFabFactory(coarse_geom, ba, dm, {4, 4, 4},
                                amrex::EBSupport::full);
    amrex::MultiFab mf(ba, dm, 1, 0, amrex::MFInfo(), *factory);
    mf.setVal(0.0);

    int size = 1;

    amrex::Vector<std::string> varnames{};
    varnames.push_back("Dummy");
    amrex::Vector<const amrex::MultiFab*> mfs(size);
    amrex::Vector<amrex::Geometry> geoms(size);
    amrex::Vector<int> level_steps(size);
    amrex::Vector<amrex::IntVect> ref_ratio(size);
    const double time_point = 0.0;

    mfs[0] = &mf;
    geoms[0] = coarse_geom;
    level_steps[0] = 0;
    ref_ratio[0] = amrex::IntVect(1);
    amrex::EB_WriteMultiLevelPlotfile("Plotfile", 1, mfs, varnames, geoms,
                                      time_point, level_steps, ref_ratio);
  }
  amrex::Finalize();
}


// implementation of the polygon distance function

struct Point {
  double x;
  double y;
};

double Dot(Point p, Point q) noexcept { return p.x * q.x + p.y * q.y; }

double Distance(Point p, Point q) noexcept {
  Point diff{q.x - p.x, q.y - p.y};
  return std::sqrt(Dot(diff, diff));
}

double Distance(Point v, Point w, Point p) {
  const Point w_minus_v{w.x - v.x, w.y - v.y};
  const double length_squared = Dot(w_minus_v, w_minus_v);
  const Point p_minus_v{p.x - v.x, p.y - v.y};

  const double t =
      std::clamp(Dot(p_minus_v, w_minus_v) / length_squared, 0.0, 1.0);

  const Point proj{(1.0 - t) * v.x + t * w.x, (1.0 - t) * v.y + t * w.y};

  const double distance =
      (length_squared > 0.0) ? Distance(p, proj) : Distance(p, v);

  return distance;
}

int CountCrossedLines(Point v, Point w, Point p) {
  const double t = ((p.y - v.y) / (w.y - v.y));
  const double xIntercept = (1.0 - t) * v.x + t * w.x;
  return ((v.y < p.y && p.y <= w.y) || (w.y < p.y && p.y <= v.y)) &&
         (xIntercept < p.x);
}

double Polygon::operator()(double x0, double y0) const noexcept {
  const double* first_x = xs_.data();
  const double* first_y = ys_.data();
  const double* last_x = xs_.data() + xs_.size();

  double min_distance = std::numeric_limits<double>::max();

  int num_crossed_lines = 0;

  const Point p{x0, y0};

  std::ptrdiff_t count = last_x - first_x;
  while ((first_x + 1) != last_x) {
    const Point v{*first_x, *first_y};
    const Point w{*(first_x + 1), *(first_y + 1)};

    num_crossed_lines += CountCrossedLines(v, w, p);
    min_distance = std::min(min_distance, Distance(v, w, p));

    ++first_x;
    ++first_y;
  }
  const int sign = num_crossed_lines % 2 ? 1 : -1;
  return sign * min_distance;
}
