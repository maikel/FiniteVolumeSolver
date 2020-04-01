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

// We give names to some magic zeros and ones.
inline constexpr int no_ghosts = 0;
inline constexpr int one_ghost_cell_width = 1;
inline constexpr int one_component = 1;

namespace fub::amrex {

template <int Rank>
ViscositySourceTerm<Rank>::ViscositySourceTerm(const PerfectGas<Rank>& eq,
                  const std::shared_ptr<GriddingAlgorithm>& grid,
                  const IndexMapping<Equation>& index)
    : equation_(eq), grid_(std::move(grid)), index_(index) {
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
void ViscositySourceTerm<Rank>::PreAdvanceLevel(int level, Duration dt,
                                                std::pair<int, int> subcycle) {

  if (level == 0) {
    const fub::amrex::PatchHierarchy& hier = grid_->GetPatchHierarchy();
    const int nlevels = hier.GetNumberOfLevels();
    const ::amrex::Geometry geom = hier.GetGeometry(0);
    ::amrex::Vector<::amrex::Array<::amrex::LinOpBCType,AMREX_SPACEDIM>> r(2);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (geom.isPeriodic(dir)) {
        r[0][dir] = ::amrex::LinOpBCType::Periodic;
        r[1][dir] = ::amrex::LinOpBCType::Periodic;
//         r[2][dir] = ::amrex::LinOpBCType::Periodic;
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
      const ::amrex::Vector<const ::amrex::MultiFab*> data = hier.GetData();
      ::amrex::Vector<::amrex::MultiFab> rhs(nlevels);

      m_reg_solve_op->setScalars(1.0, dt.count());
      for (int lev = 0; lev < nlevels; ++lev) {
        const ::amrex::BoxArray& ba = data[lev]->boxArray();
        const ::amrex::DistributionMapping& dm = data[lev]->DistributionMap();
        ::amrex::MultiFab density(ba, dm, one_component, 0);
        ::amrex::MultiFab::Copy(density, *data[lev], index_.density, 0,
                       one_component, 0);

        m_reg_solve_op->setACoeffs(lev, density);
        // TODO: We might want to introduce a space dependent eta
//         ::amrex::Array<MultiFab,AMREX_SPACEDIM> b = m_incflo->average_velocity_eta_to_faces(lev, *eta[lev]);
//         m_reg_solve_op->setShearViscosity(lev, GetArrOfConstPtrs(b));
        m_reg_solve_op->setShearViscosity(lev, eta_);

        rhs[lev].define(ba, dm, AMREX_SPACEDIM, 0);
        ::amrex::MultiFab::Copy(rhs[lev], *data[lev], index_.momentum[0], 0,
                       index_.momentum.size(), 0);


        ::amrex::MultiFab velocity(ba, dm, AMREX_SPACEDIM, 1);
        ::amrex::MultiFab::Copy(velocity, *data[lev], index_.momentum[0], 0,
                       index_.momentum.size(), 0);

        for (std::size_t i = 0; i < index_.momentum.size(); ++i) {
          ::amrex::MultiFab::Divide(velocity, *data[lev], i, index_.density,
                       one_component, 0);
        }


        const ::amrex::MultiFab* vel = &velocity;
        m_reg_solve_op->setLevelBC(lev, vel);
    }

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

    ::amrex::Vector<::amrex::MultiFab*> vel_corr(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        const ::amrex::BoxArray& ba = data[lev]->boxArray();
        const ::amrex::DistributionMapping& dm = data[lev]->DistributionMap();
        vel_corr[lev] = (new ::amrex::MultiFab(ba, dm, AMREX_SPACEDIM, 0));
    }

    mlmg.solve(vel_corr, ::amrex::GetVecOfConstPtrs(rhs), m_mg_rtol, m_mg_atol);
  }
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
