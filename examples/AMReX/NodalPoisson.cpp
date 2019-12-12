#include "fub/AMReX/CartesianGridGeometry.hpp"
#include "fub/AMReX/ScopeGuard.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/core/mdspan.hpp"

#include <AMReX.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

int main(int argc, char** argv) {
  fub::amrex::ScopeGuard guard(argc, argv);

  amrex::RealBox rb({0., 0.}, {1., 1.});
  amrex::Array<int, AMREX_SPACEDIM> is_periodic{1, 1};

  amrex::Box domain0({0, 0}, {63, 63});
  amrex::Box domain1 = domain0;
  domain1.refine(2);

  amrex::Box box1({32, 32}, {63, 63});

  amrex::Geometry geom0(domain0, &rb, -1, is_periodic.data());
  amrex::Geometry geom1(domain1, &rb, -1, is_periodic.data());

  amrex::BoxArray box_array0(domain0);
  box_array0.maxSize(8);

  amrex::BoxArray box_array1(box1);

  amrex::BoxArray node_boxes0 = amrex::convert(box_array0, {1, 1});
  amrex::BoxArray node_boxes1 = amrex::convert(box_array1, {1, 1});

  amrex::DistributionMapping dmap0(box_array0);
  amrex::DistributionMapping dmap1(box_array1);

  const int ncomp = 1;
  const int no_ghosts = 0;
  amrex::MultiFab solution(node_boxes0, dmap0, ncomp, no_ghosts);
  amrex::MultiFab rhs(node_boxes0, dmap0, ncomp, no_ghosts);
  amrex::MultiFab exact_solution(node_boxes0, dmap0, ncomp, no_ghosts);
  amrex::MultiFab sigma(box_array0, dmap0, ncomp, no_ghosts);

  amrex::MultiFab solution1(node_boxes1, dmap1, ncomp, no_ghosts);
  amrex::MultiFab rhs1(node_boxes1, dmap1, ncomp, no_ghosts);
  amrex::MultiFab exact_solution1(node_boxes1, dmap1, ncomp, no_ghosts);
  amrex::MultiFab sigma1(box_array1, dmap1, ncomp, no_ghosts);

  solution.setVal(0.0);
  solution1.setVal(0.0);
  sigma.setVal(1.0);
  sigma1.setVal(1.0);

  const double* dx = geom0.CellSize();
  const double* dx1 = geom1.CellSize();

  for (amrex::MFIter mfi(rhs, true); mfi.isValid(); ++mfi) {
    amrex::FArrayBox& fab = rhs[mfi];
    fub::mdspan<double, 2> fabv = fub::amrex::MakeMdSpan(fab, 0);

    amrex::Box box = rhs[mfi].box();
    int i0 = box.smallEnd(0);
    int j0 = box.smallEnd(1);
    // fub::CartesianCoordinates coords =
    // fub::amrex::GetCartesianCoordinates(geom0, box);
    for (int j = 0; j < fabv.extent(1); ++j) {
      for (int i = 0; i < fabv.extent(0); ++i) {
        const std::array<double, 2> x{(i + i0) * dx[0], (j + j0) * dx[1]};
        fabv(i, j) =
            -6 * M_PI * std::cos(2 * M_PI * x[0]) * std::cos(2 * M_PI * x[1]) -
            6 * M_PI * std::cos(4 * M_PI * x[0]) * std::cos(4 * M_PI * x[1]);
      }
    }

    // fab.setFormat(amrex::FABio::FAB_ASCII);
    // fab.writeOn(std::cout);
    // std::cout << '\n';
  }

  for (amrex::MFIter mfi(rhs1, true); mfi.isValid(); ++mfi) {
    amrex::FArrayBox& fab = rhs1[mfi];
    fub::mdspan<double, 2> fabv = fub::amrex::MakeMdSpan(fab, 0);

    amrex::Box box = rhs1[mfi].box();
    int i0 = box.smallEnd(0);
    int j0 = box.smallEnd(1);
    // fub::CartesianCoordinates coords =
    // fub::amrex::GetCartesianCoordinates(geom0, box);
    for (int j = 0; j < fabv.extent(1); ++j) {
      for (int i = 0; i < fabv.extent(0); ++i) {
        const std::array<double, 2> x{(i + i0) * dx1[0], (j + j0) * dx1[1]};
        fabv(i, j) =
            -6 * M_PI * std::cos(2 * M_PI * x[0]) * std::cos(2 * M_PI * x[1]) -
            6 * M_PI * std::cos(4 * M_PI * x[0]) * std::cos(4 * M_PI * x[1]);
      }
    }

    // fab.setFormat(amrex::FABio::FAB_ASCII);
    // fab.writeOn(std::cout);
    // std::cout << '\n';
  }

  // ::amrex::WriteMultiLevelPlotfile(plotfilename, nlevels, mf, varnames,
  // geoms,
  //                                    time_point, level_steps, ref_ratio);

  amrex::IntVect unit = amrex::IntVect::TheUnitVector();
  amrex::WriteMultiLevelPlotfile("rhs", 2, {&rhs, &rhs1}, {"rhs"},
                                 {geom0, geom1}, 0.0, {0, 0}, {unit, 2 * unit});

  amrex::MLNodeLaplacian linop({geom0, geom1}, {box_array0, box_array1},
                               {dmap0, dmap1});

  linop.setDomainBC(
      {AMREX_D_DECL(amrex::LinOpBCType::Periodic, amrex::LinOpBCType::Periodic,
                    amrex::LinOpBCType::Periodic)},
      {AMREX_D_DECL(amrex::LinOpBCType::Periodic, amrex::LinOpBCType::Periodic,
                    amrex::LinOpBCType::Periodic)});

  linop.setLevelBC(0, nullptr);

  linop.setSigma(0, sigma);

  amrex::MLMG mlmg(linop);
  mlmg.setVerbose(2);
  mlmg.setBottomVerbose(2);

  mlmg.solve({&solution, &solution1}, {&rhs, &rhs1}, 1e-6, 0.0);

  amrex::WriteMultiLevelPlotfile("solution", 2, {&solution, &solution1},
                                 {"solution"}, {geom0, geom1}, 0.0, {0, 0},
                                 {unit, 2 * unit});
}
