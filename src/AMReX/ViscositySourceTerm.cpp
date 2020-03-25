// Copyright (c) 2020 Stefan Vater
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

#include "fub/AMReX/ViscositySourceTerm.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/equations/PerfectGas.hpp"

namespace fub::amrex {

template <int Rank>
ViscositySourceTerm<Rank>::ViscositySourceTerm(const PerfectGas<Rank>& eq,
                                 const std::shared_ptr<GriddingAlgorithm>& grid)
    : equation_(eq), grid_(grid) {
  ResetHierarchyConfiguration(grid);

  info_solve.setMaxCoarseningLevel(m_mg_max_coarsening_level);
  info_apply.setMaxCoarseningLevel(0);

}

// ViscositySourceTerm::ViscositySourceTerm(const ViscositySourceTerm& other)
//     : equation_(other.equation_) {
//   Ax_.reserve(other.Ax_.size());
//   for (const ::amrex::MultiFab& other_Ax : other.Ax_) {
//     ::amrex::BoxArray box_array = other_Ax.boxArray();
//     ::amrex::DistributionMapping distribution_map = other_Ax.DistributionMap();
//     const int ncomp = other_Ax.nComp();
//     const int ngrow = other_Ax.nGrow();
//     ::amrex::MultiFab& Ax =
//         Ax_.emplace_back(box_array, distribution_map, ncomp, ngrow);
//     Ax.copy(other_Ax);
//   }
// }
//
// ViscositySourceTerm& ViscositySourceTerm::operator=(const ViscositySourceTerm& other) {
//   ViscositySourceTerm tmp(other);
//   return (*this = std::move(tmp));
// }

namespace {
} // namespace

template <int Rank>
void ViscositySourceTerm<Rank>::PreAdvanceHierarchy() {
  ::amrex::Print() << "\nHello!\n\n";

  const fub::amrex::PatchHierarchy& hier = grid_->GetPatchHierarchy();
  const auto nlevels = hier.GetNumberOfLevels();
  const ::amrex::Geometry geom = hier.GetGeometry(0);
  ::amrex::Vector<::amrex::Array<::amrex::LinOpBCType,AMREX_SPACEDIM>> r(3);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (geom.isPeriodic(dir)) {
      r[0][dir] = ::amrex::LinOpBCType::Periodic;
      r[1][dir] = ::amrex::LinOpBCType::Periodic;
        r[2][dir] = ::amrex::LinOpBCType::Periodic;
    }
  }

  m_reg_solve_op.reset(new ::amrex::MLTensorOp(hier.GetGeometries(),
                                      hier.GetBoxArrays(),
                                      hier.GetDistributionMappings(),
                                      info_solve));
  m_reg_solve_op->setMaxOrder(m_mg_maxorder);
  // TODO: set correct bc
  m_reg_solve_op->setDomainBC(r, r);

//   if (m_incflo->need_divtau()) {
    m_reg_apply_op.reset(new ::amrex::MLTensorOp(hier.GetGeometries(),
                                      hier.GetBoxArrays(),
                                      hier.GetDistributionMappings(),
                                      info_apply));
    m_reg_apply_op->setMaxOrder(m_mg_maxorder);
    // TODO: set correct bc
    m_reg_apply_op->setDomainBC(r, r);
// }

    // TODO: set correct dt!
    const double dt = 1.0;
    const ::amrex::Vector<const ::amrex::MultiFab*> data = hier.GetData();

    m_reg_solve_op->setScalars(1.0, dt);
//     for (int lev = 0; lev < nlevels; ++lev) {
//       m_reg_solve_op->setACoeffs(lev, *density[lev]);
//       ::amrex::Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_velocity_eta_to_faces(lev, *eta[lev]);
//       m_reg_solve_op->setShearViscosity(lev, GetArrOfConstPtrs(b));
//     }

    ::amrex::Vector<::amrex::MultiFab> rhs(nlevels);
//     for (int lev = 0; lev < nlevels; ++lev) {
//         rhs[lev].define(velocity[lev]->boxArray(),
//                         velocity[lev]->DistributionMap(), AMREX_SPACEDIM, 0);
//
//         for (MFIter mfi(rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
//             Box const& bx = mfi.tilebox();
//             Array4<Real> const& rhs_a = rhs[lev].array(mfi);
//             Array4<Real const> const& vel_a = velocity[lev]->const_array(mfi);
//             Array4<Real const> const& rho_a = density[lev]->const_array(mfi);
//             amrex::ParallelFor(bx,AMREX_SPACEDIM,
//             [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//             {
//                 rhs_a(i,j,k,n) = rho_a(i,j,k) * vel_a(i,j,k,n);
//             });
//         }
//
//         m_reg_solve_op->setLevelBC(lev, velocity[lev]);
//     }

    ::amrex::MLMG mlmg(*m_reg_solve_op);

//     // The default bottom solver is BiCG
//     if (m_bottom_solver_type == "smoother")
//     {
//         mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
//     }
//     else if (m_bottom_solver_type == "hypre")
//     {
//         mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
//     }
//     // Maximum iterations for MultiGrid / ConjugateGradients
    mlmg.setMaxIter(m_mg_max_iter);
    mlmg.setMaxFmgIter(m_mg_max_fmg_iter);
    mlmg.setCGMaxIter(m_mg_cg_maxiter);

    // Verbosity for MultiGrid / ConjugateGradients
    mlmg.setVerbose(m_mg_verbose);
    mlmg.setCGVerbose(m_mg_cg_verbose);

//     mlmg.solve(velocity, ::amrex::GetVecOfConstPtrs(rhs), m_mg_rtol, m_mg_atol);
}

template <int Rank>
void ViscositySourceTerm<Rank>::ResetHierarchyConfiguration(
    const std::shared_ptr<amrex::GriddingAlgorithm>& grid) {
}

template <int Rank>
Duration ViscositySourceTerm<Rank>::ComputeStableDt(int /* level */) {
  return Duration(std::numeric_limits<double>::max());
}

template <int Rank>
Result<void, TimeStepTooLarge>
ViscositySourceTerm<Rank>::AdvanceLevel(IntegratorContext& simulation_data, int level,
                              Duration time_step_size) {
  ::amrex::MultiFab& data = simulation_data.GetScratch(level);
  const double dt = time_step_size.count();

  return boost::outcome_v2::success();
}

template class ViscositySourceTerm<1>;
template class ViscositySourceTerm<2>;
template class ViscositySourceTerm<3>;

} // namespace fub::amrex
