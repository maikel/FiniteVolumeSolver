#include "fub/AMReX/CartesianGridGeometry.hpp"
#include "fub/AMReX/ScopeGuard.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/MLMG/MLNodeHelmDualCstVel.hpp"
#include "fub/core/mdspan.hpp"

#include <AMReX.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

int main(int argc, char** argv) {
  fub::amrex::ScopeGuard guard(argc, argv);

  int mg_verbose     = 4;
  int bottom_verbose = 4;
  int max_iter       = 100;
  int n_cell         = 32;
  int max_grid_size  = 32;
  amrex::Real reltol = 1.e-10;
  amrex::Real abstol = 1.e-15;

        // Define the absolute tolerance; note that this argument is optional
//         Real abstol = 1.e-15;


  amrex::RealBox rb({AMREX_D_DECL(0., 0., 0.)}, {AMREX_D_DECL(1., 1., 1.)});
  amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

  amrex::Box domain(amrex::IntVect{AMREX_D_DECL(0, 0, 0)},
                    amrex::IntVect{AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)});

  amrex::Geometry geom(domain, &rb, amrex::CoordSys::cartesian, is_periodic.data());

  amrex::BoxArray box_array(domain);
  box_array.maxSize(max_grid_size);
  amrex::BoxArray node_boxes = amrex::convert(box_array, {AMREX_D_DECL(1, 1, 1)});

  amrex::DistributionMapping dmap(box_array);

  // Store plotfile variables
  amrex::MultiFab plotfile_mf(box_array, dmap, 2*AMREX_SPACEDIM+2, 0);

  const int no_ghosts = 0;

  //  Create the cell-centered velocity field we want to project
  amrex::MultiFab vel(box_array, dmap, AMREX_SPACEDIM, 1);

  amrex::MultiFab phi(node_boxes, dmap, 1, no_ghosts);
  amrex::MultiFab rhs(node_boxes, dmap, 1, no_ghosts);
  amrex::MultiFab div(node_boxes, dmap, 1, no_ghosts);
  amrex::MultiFab sigma(box_array, dmap, 1, no_ghosts);

  phi.setVal(0.0);
  div.setVal(0.0);
  sigma.setVal(1.0);

  const auto dx = geom.CellSizeArray();

  for (amrex::MFIter mfi(vel, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.tilebox();
    amrex::Array4<double> const v = vel.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      constexpr double pi = 3.1415926535897932;
      constexpr double tpi = 2.*pi;
      constexpr double fpi = 4.*pi;
      constexpr double fac = tpi*tpi*AMREX_SPACEDIM;

      double x = i*dx[0];
      double y = j*dx[1];
      double z = k*dx[2];

      v(i,j,k,0) = (std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z))
          + 0.25 * (std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z));
      v(i,j,k,1) = (std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z))
          + 0.25 * (std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z));

    });
  }

  // copy velocity into plotfile
  amrex::MultiFab::Copy(plotfile_mf, vel, 0, 0, AMREX_SPACEDIM, 0);

  // Setup linear operator, AKA the nodal Laplacian
  amrex::LPInfo lp_info;
  //set number of MG levels to 1 (effectively no MG)
  lp_info.setMaxCoarseningLevel(0);

  amrex::MLNodeHelmDualCstVel linop({geom}, {box_array}, {dmap}, lp_info);

  linop.setDomainBC(
      {AMREX_D_DECL(amrex::LinOpBCType::Periodic, amrex::LinOpBCType::Periodic,
                    amrex::LinOpBCType::Periodic)},
      {AMREX_D_DECL(amrex::LinOpBCType::Periodic, amrex::LinOpBCType::Periodic,
                    amrex::LinOpBCType::Periodic)});

  linop.setSigma(0, sigma);

  // cell-centered contributions to RHS
  amrex::MultiFab S_cc(box_array, dmap, 1, 1);
  S_cc.setVal(0.0); // Set it to zero for this example

  // node-centered contributions to RHS
  amrex::MultiFab S_nd(node_boxes, dmap, 1, 1);
  S_nd.setVal(0.0); // Set it to zero for this example

  // Compute RHS -- vel must be cell-centered
  linop.compRHS({&rhs}, {&vel}, {&S_nd}, {&S_cc});

  amrex::MLMG nodal_solver(linop);
  nodal_solver.setVerbose(mg_verbose);
  nodal_solver.setBottomVerbose(bottom_verbose);
  nodal_solver.setMaxIter(max_iter);

// if we want to use another bottom solver
//   nodal_solver.setBottomSolver(BottomSolver::smoother);
//   nodal_solver.setBottomSolver(BottomSolver::cg);

  nodal_solver.solve({&phi}, {&rhs}, reltol, abstol);

  // create cell-centered multifab to hold value of -sigma*grad(phi) at cell-centers
  amrex::MultiFab fluxes(box_array, dmap, AMREX_SPACEDIM, 1);
  fluxes.setVal(0.0);

  // get fluxes from solver
  nodal_solver.getFluxes( {&fluxes} );

  // apply projection explicitly --  vel = vel - sigma * grad(phi)
  amrex::MultiFab::Add(vel, fluxes, 0, 0, AMREX_SPACEDIM, 0);

  // copy velocity into plotfile
  amrex::average_node_to_cellcenter(plotfile_mf, AMREX_SPACEDIM, rhs, 0, 1);
  amrex::MultiFab::Copy(plotfile_mf, vel, 0, AMREX_SPACEDIM+1, AMREX_SPACEDIM, 0);

  // Compute divergence -- vel must be cell-centered
  linop.compDivergence({&div}, {&vel});
  amrex::average_node_to_cellcenter(plotfile_mf, 2*AMREX_SPACEDIM+1, div, 0, 1);


  std::string base_name = "NodalProjectionSingleLevel/";
  std::string name = base_name + "plot";

  amrex::Print() << "Start output to '" << name << "'.\n";
  WriteSingleLevelPlotfile(name, plotfile_mf,
    { "xvelold", "yvelold", "rhsold" ,"xvelnew", "yvelnew", "rhsnew" },
    geom, 0.0, 0);
  amrex::Print() << "Finished output to '" << name << "'.\n";

}
