#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_LSCore.H>

#include <Eigen/Dense>
#include <memory>

Eigen::Vector2d OrthogonalTo(const Eigen::Vector2d& x) {
  return {x[1], -x[0]};
}

auto Ramp(const Eigen::Vector2d& p1) {
  Eigen::Vector2d p0{0.0, 0.0};
  Eigen::Vector2d norm1 = OrthogonalTo(p1 - p0).normalized();
  return amrex::EB2::PlaneIF({p0[0], p0[1]}, {norm1[0], norm1[1]}, false);
}

amrex::Geometry MakeGeometry(const std::array<int, 2>& n_cells,
                             const std::array<double, 2>& lower,
                             const std::array<double, 2>& upper) {
  amrex::Box box(amrex::IntVect{0, 0},
                 amrex::IntVect{n_cells[0] - 1, n_cells[1] - 1});
  amrex::RealBox realbox{lower, upper};
  return amrex::Geometry(box, realbox, -1, std::array<int, 2>{});
}

void MyMain() {
  const std::array<int, 2> n_cells{128, 128};
  const std::array<double, 2> xlower{-1.0, -1.0};
  const std::array<double, 2> xupper{+1.0, +1.0};
  const amrex::Geometry geom = MakeGeometry(n_cells, xlower, xupper);
  amrex::EB2::Build(amrex::EB2::makeShop(Ramp({-1.0, +1.0})), geom, 1, 1);

  amrex::BoxList list;
  list.push_back(amrex::Box({104, 0}, {127, 39}));
  list.push_back(amrex::Box({88, 8}, {103, 23}));
  list.push_back(amrex::Box({72, 24}, {103, 39}));
  list.push_back(amrex::Box({56, 40}, {71, 55}));
  list.push_back(amrex::Box({72, 40}, {127, 71}));
  list.push_back(amrex::Box({24, 72}, {39, 87}));
  list.push_back(amrex::Box({40, 56}, {71, 87}));
  list.push_back(amrex::Box({96, 72}, {127, 103}));
  list.push_back(amrex::Box({8, 88}, {55, 103}));
  list.push_back(amrex::Box({0, 104}, {39, 127}));
  amrex::BoxArray ba(list);
  amrex::DistributionMapping dm(ba);

  amrex::Vector<int> ngrow_per_support_level(4, 4);
  std::unique_ptr<amrex::EBFArrayBoxFactory> factory = amrex::makeEBFabFactory(geom, ba, dm, ngrow_per_support_level, amrex::EBSupport::full);

  const amrex::FabArray<amrex::EBCellFlagFab>& cellflags = factory->getMultiEBCellFlagFab();
  amrex::MultiCutFab mcf(ba, dm, 1, 4, cellflags);
  mcf.setVal(42.0);

  amrex::MultiCutFab copy(ba, dm, 1, 4, cellflags);
  copy.ParallelCopy(mcf, 0, 0, 1, 4, 4);
}

int main(int argc, char** argv) {
  ::amrex::Initialize(argc, argv);
  MyMain();
  ::amrex::Finalize();
}
