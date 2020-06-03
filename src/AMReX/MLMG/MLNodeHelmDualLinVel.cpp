// Copyright (c) 2017, The Regents of the University of California,
// through Lawrence Berkeley National Laboratory and the Alliance for
// Sustainable Energy, LLC., through National Renewable Energy Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of
// Energy).  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// (1) Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence
// Berkeley National Laboratory, Alliance for Sustainable Energy, LLC.,
// National Renewable Energy Laboratory, U.S. Dept. of Energy nor the
// names of its contributors may be used to endorse or promote products
// derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Notice: This file is copied and modified from AMReX-Codes
// (https://amrex-codes.github.io/).

// clang-format off
#include <limits>

#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>

#include <fub/AMReX/MLMG/MLNodeHelmDualLinVel.hpp>
#include <src/AMReX/MLMG/MLNodeHelmDualLinVel_K.cpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-parameter"

namespace amrex {

MLNodeHelmDualLinVel::MLNodeHelmDualLinVel (const Vector<Geometry>& a_geom,
                                  const Vector<BoxArray>& a_grids,
                                  const Vector<DistributionMapping>& a_dmap,
                                  const LPInfo& a_info,
                                  const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

MLNodeHelmDualLinVel::~MLNodeHelmDualLinVel ()
{}

void
MLNodeHelmDualLinVel::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info,
                         const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLNodeHelmDualLinVel::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info, a_factory);

    m_sigma.resize(m_num_amr_levels);
    m_sigmacross.resize(m_num_amr_levels);
    m_alpha.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_sigma[amrlev].resize(m_num_mg_levels[amrlev]);
        m_sigmacross[amrlev].resize(m_num_mg_levels[amrlev]);
        m_alpha[amrlev].resize(m_num_mg_levels[amrlev]);
        const int mglev = 0;
        const int idim = 0;
        m_sigma[amrlev][mglev][idim].reset
            (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1,
                          MFInfo(), *m_factory[amrlev][0]));
        m_sigma[amrlev][mglev][idim]->setVal(0.0);
        m_sigmacross[amrlev][mglev][idim].reset
            (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], AMREX_SPACEDIM*(AMREX_SPACEDIM-1), 1));
        m_sigmacross[amrlev][mglev][idim]->setVal(0.0);
        const BoxArray& ba = amrex::convert(m_grids[amrlev][mglev],
                                            IntVect(AMREX_D_DECL(1,1,1)));
        m_alpha[amrlev][mglev].define(ba, m_dmap[amrlev][mglev], 1, 0);
        m_alpha[amrlev][mglev].setVal(0.0);
    }

}

void
MLNodeHelmDualLinVel::setSigma (int amrlev, const MultiFab& a_sigma)
{
    MultiFab::Copy(*m_sigma[amrlev][0][0], a_sigma, 0, 0, 1, 0);
}

void
MLNodeHelmDualLinVel::setSigmaCross (int amrlev, const MultiFab& a_sigmacross)
{
    MultiFab::Copy(*m_sigmacross[amrlev][0][0], a_sigmacross, 0, 0, AMREX_SPACEDIM*(AMREX_SPACEDIM-1), 0);
}

void
MLNodeHelmDualLinVel::setAlpha(int amrlev, const MultiFab& alpha)
{
    MultiFab::Copy(m_alpha[amrlev][0], alpha, 0, 0, 1, 0);
    m_is_bottom_singular = false;
}

void
MLNodeHelmDualLinVel::compDivergence (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel)
{
    compRHS(rhs, vel, Vector<const MultiFab*>(), Vector<MultiFab*>());
}

void
MLNodeHelmDualLinVel::compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
                          const Vector<const MultiFab*>& rhnd,
                          const Vector<MultiFab*>& a_rhcc)
{
    BL_PROFILE("MLNodeHelmDualLinVel::compRHS()");

    if (!m_masks_built) buildMasks();

    const auto lobc = LoBC();
    const auto hibc = HiBC();

    Vector<std::unique_ptr<MultiFab> > rhcc(m_num_amr_levels);
    Vector<std::unique_ptr<MultiFab> > rhs_cc(m_num_amr_levels);

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        const Geometry& geom = m_geom[ilev][0];
        AMREX_ASSERT(vel[ilev]->nComp() >= AMREX_SPACEDIM);
        AMREX_ASSERT(vel[ilev]->nGrow() >= 1);
        vel[ilev]->FillBoundary(0, AMREX_SPACEDIM, IntVect(1), geom.periodicity());

        if (ilev < a_rhcc.size() && a_rhcc[ilev])
        {
            rhcc[ilev].reset(new MultiFab(a_rhcc[ilev]->boxArray(),
                                          a_rhcc[ilev]->DistributionMap(), 1, 1));
            rhcc[ilev]->setVal(0.0);
            MultiFab::Copy(*rhcc[ilev], *a_rhcc[ilev], 0, 0, 1, 0);
            rhcc[ilev]->FillBoundary(geom.periodicity());

            rhs_cc[ilev].reset(new MultiFab(rhs[ilev]->boxArray(),
                                            rhs[ilev]->DistributionMap(), 1, 0));
        }

        const auto dxinvarr = geom.InvCellSizeArray();
        const Box& nddom = amrex::surroundingNodes(geom.Domain());

        const iMultiFab& dmsk = *m_dirichlet_mask[ilev][0];

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*rhs[ilev],mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& rhsarr = rhs[ilev]->array(mfi);
            Array4<Real const> const& velarr = vel[ilev]->const_array(mfi);
            Array4<int const> const& dmskarr = dmsk.const_array(mfi);

            {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                {

                    mlndhelm_divu(i,j,k,rhsarr,velarr,dmskarr,dxinvarr);
                });
            }

            mlndhelm_impose_neumann_bc(bx, rhsarr, nddom, lobc, hibc);

            if (rhcc[ilev])
            {
                Array4<Real> const& rhs_cc_a = rhs_cc[ilev]->array(mfi);
                Array4<Real const> const& rhccarr = rhcc[ilev]->const_array(mfi);
                {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        rhs_cc_a(i,j,k) = mlndhelm_rhcc(i, j, k, rhccarr, dmskarr);
                    });
                }

                mlndhelm_impose_neumann_bc(bx, rhs_cc_a, nddom, lobc, hibc);
            }
        }
    }

    Vector<std::unique_ptr<MultiFab> > frhs(m_num_amr_levels);

    for (int ilev = 0; ilev < m_num_amr_levels-1; ++ilev)
    {
        const Geometry& cgeom = m_geom[ilev  ][0];
        const Geometry& fgeom = m_geom[ilev+1][0];

        frhs[ilev].reset(new MultiFab(amrex::coarsen(rhs[ilev+1]->boxArray(),2),
                                      rhs[ilev+1]->DistributionMap(), 1, 0));
        frhs[ilev]->setVal(0.0);

        const Box& cccdom = cgeom.Domain();
        const auto fdxinv = fgeom.InvCellSizeArray();
        const iMultiFab& fdmsk = *m_dirichlet_mask[ilev+1][0];

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            FArrayBox vfab, rfab, rhccfab;
            for (MFIter mfi(*frhs[ilev],mfi_info); mfi.isValid(); ++mfi)
            {
                const Box& cvbx = mfi.validbox();
                const Box& fvbx = amrex::refine(cvbx,2);
                const Box& cbx = mfi.tilebox();
                const Box& fbx = amrex::refine(cbx,2);

                const Box& cc_fbx = amrex::enclosedCells(fbx);
                const Box& cc_fvbx = amrex::enclosedCells(fvbx);

                const Box& bx_vel = amrex::grow(cc_fbx,2) & amrex::grow(cc_fvbx,1);
                Box b = bx_vel & cc_fvbx;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    if (m_lobc[0][idim] == LinOpBCType::inflow)
                    {
                        if (b.smallEnd(idim) == cccdom.smallEnd(idim)) {
                            b.growLo(idim, 1);
                        }
                    }
                    if (m_hibc[0][idim] == LinOpBCType::inflow)
                    {
                        if (b.bigEnd(idim) == cccdom.bigEnd(idim)) {
                            b.growHi(idim, 1);
                        }
                    }
                }

                vfab.resize(bx_vel, AMREX_SPACEDIM);
                Elixir veli = vfab.elixir();
                Array4<Real> const& varr = vfab.array();

                const Box& bx_rhs = amrex::grow(fbx,1);
                const Box& b2 = bx_rhs & amrex::grow(fvbx,-1);
                rfab.resize(bx_rhs);
                Elixir reli = rfab.elixir();
                Array4<Real> const& rarr = rfab.array();

                Array4<Real const> const& varr_orig = vel[ilev+1]->const_array(mfi);
                AMREX_HOST_DEVICE_FOR_4D(bx_vel, AMREX_SPACEDIM, i, j, k, n,
                {
                    if (b.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                        varr(i,j,k,n) = varr_orig(i,j,k,n);
                    } else {
                        varr(i,j,k,n) = 0.0;
                    }
                });

                Array4<Real const> const& rarr_orig = rhs[ilev+1]->const_array(mfi);
                AMREX_HOST_DEVICE_FOR_3D(bx_rhs, i, j, k,
                {
                    if (b2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                        rarr(i,j,k) = rarr_orig(i,j,k);
                    } else {
                        rarr(i,j,k) = 0.0;
                    }
                    mlndhelm_divu_compute_fine_contrib(i,j,k,fvbx,rarr,varr,fdxinv);
                });

                Array4<Real> const& rhsarr = frhs[ilev]->array(mfi);
                Array4<int const> const& mskarr = fdmsk.const_array(mfi);
                AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                {
                    mlndhelm_divu_add_fine_contrib(i,j,k,fvbx,rhsarr,rarr,mskarr);
                });

                if (rhcc[ilev+1])
                {
                    const Box& bx_rhcc = amrex::grow(cc_fbx,2);
                    const Box& b3 = bx_rhcc & cc_fvbx;

                    rhccfab.resize(bx_rhcc);
                    Elixir rhcceli = rhccfab.elixir();
                    Array4<Real> const& rhccarr = rhccfab.array();

                    Array4<Real const> const& rhccarr_orig = rhcc[ilev+1]->const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(bx_rhcc, i, j, k,
                    {
                        if (b3.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                            rhccarr(i,j,k) = rhccarr_orig(i,j,k);
                        } else {
                            rhccarr(i,j,k) = 0.0;
                        }
                    });

                    AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                    {
                        mlndhelm_rhcc_fine_contrib(i,j,k,fvbx,rhsarr,rhccarr,mskarr);
                    });
                }
            }
        }
    }

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        if (rhs_cc[ilev]) {
            MultiFab::Add(*rhs[ilev], *rhs_cc[ilev], 0, 0, 1, 0);
        }
    }

    for (int ilev = m_num_amr_levels-2; ilev >= 0; --ilev)
    {
        const Geometry& cgeom = m_geom[ilev][0];

        MultiFab crhs(rhs[ilev]->boxArray(), rhs[ilev]->DistributionMap(), 1, 0);
        crhs.setVal(0.0);
        crhs.ParallelAdd(*frhs[ilev], cgeom.periodicity());

        const Box& cccdom = cgeom.Domain();
        const Box& cnddom = amrex::surroundingNodes(cccdom);
        const auto cdxinv = cgeom.InvCellSizeArray();
        const iMultiFab& cdmsk = *m_dirichlet_mask[ilev][0];
        const iMultiFab& c_nd_mask = *m_nd_fine_mask[ilev];
        const iMultiFab& c_cc_mask = *m_cc_fine_mask[ilev];
        const auto& has_fine_bndry = *m_has_fine_bndry[ilev];

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*rhs[ilev]); mfi.isValid(); ++mfi)
        {
            if (has_fine_bndry[mfi])
            {
                const Box& bx = mfi.tilebox();

                Array4<Real> const& rhsarr = rhs[ilev]->array(mfi);
                Array4<Real const> const& velarr = vel[ilev]->const_array(mfi);
                Array4<Real const> const& crhsarr = crhs.const_array(mfi);
                Array4<int const> const& cdmskarr = cdmsk.const_array(mfi);
                Array4<int const> const& ndmskarr = c_nd_mask.const_array(mfi);
                Array4<int const> const& ccmskarr = c_cc_mask.const_array(mfi);

                Array4<Real const> rhccarr = (rhcc[ilev])
                    ? rhcc[ilev]->const_array(mfi) : Array4<Real const>{};

                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    mlndhelm_divu_cf_contrib(i,j,k,rhsarr,velarr,crhsarr,rhccarr,
                                            cdmskarr,ndmskarr,ccmskarr,
                                            cdxinv,cnddom,lobc,hibc);
                });
            }
        }
    }

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        if (ilev < rhnd.size() && rhnd[ilev]) {
            MultiFab::Add(*rhs[ilev], *rhnd[ilev], 0, 0, 1, 0);
        }
    }

}

void
MLNodeHelmDualLinVel::updateVelocity (const Vector<MultiFab*>& vel, const Vector<MultiFab const*>& sol) const
{

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        const auto& sigma = *m_sigma[amrlev][0][0];
        const auto dxinv = m_geom[amrlev][0].InvCellSizeArray();
        for (MFIter mfi(*vel[amrlev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& varr = vel[amrlev]->array(mfi);
            Array4<Real const> const& solarr = sol[amrlev]->const_array(mfi);
            Array4<Real const> const& sigmaarr = sigma.const_array(mfi);
            {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                {
                    mlndhelm_mknewu(i,j,k,varr,solarr,sigmaarr,dxinv);
                });
            }
        }
    }
}

void
MLNodeHelmDualLinVel::getFluxes (const Vector<MultiFab*> & a_flux, const Vector<MultiFab*>& a_sol) const
{

    AMREX_ASSERT(a_flux[0]->nComp() >= AMREX_SPACEDIM);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        const auto& sigma = *m_sigma[amrlev][0][0];
        const auto dxinv = m_geom[amrlev][0].InvCellSizeArray();

        // Initialize to zero because we only want -(sigma * grad(phi))

        for (MFIter mfi(sigma, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& farr = a_flux[amrlev]->array(mfi);
            Array4<Real const> const& solarr = a_sol[amrlev]->const_array(mfi);
            Array4<Real const> const& sigmaarr = sigma.array(mfi);

            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, AMREX_SPACEDIM, i, j, k, n,
            {
                farr(i,j,k,n) = 0.0;
            });

            {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                {
                    mlndhelm_mknewu(i,j,k,farr,solarr,sigmaarr,dxinv);
                });
            }
        }
    }
}

void
MLNodeHelmDualLinVel::averageDownCoeffs ()
{
    BL_PROFILE("MLNodeHelmDualLinVel::averageDownCoeffs()");

        for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
        {
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                int ndims = (m_use_harmonic_average) ? AMREX_SPACEDIM : 1;
                for (int idim = 0; idim < ndims; ++idim)
                {
                    if (m_sigma[amrlev][mglev][idim] == nullptr) {
                        if (mglev == 0) {
                            m_sigma[amrlev][mglev][idim].reset
                                (new MultiFab(*m_sigma[amrlev][mglev][0], amrex::make_alias, 0, 1));
                        } else {
                            m_sigma[amrlev][mglev][idim].reset
                                (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
                            m_sigma[amrlev][mglev][idim]->setVal(0.0);
                        }
                    }
                    if (m_sigmacross[amrlev][mglev][idim] == nullptr) {
                        if (mglev == 0) {
                            m_sigmacross[amrlev][mglev][idim].reset
                                (new MultiFab(*m_sigmacross[amrlev][mglev][0], amrex::make_alias, 0, 1));
                        } else {
                            m_sigmacross[amrlev][mglev][idim].reset
                                (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
                            m_sigmacross[amrlev][mglev][idim]->setVal(0.0);
                        }
                    }
                }
                if (m_alpha[amrlev][mglev].nComp() == 0) {
                    if (mglev > 0) {
                        const BoxArray& ba = amrex::convert(m_grids[amrlev][mglev],
                                                            IntVect(AMREX_D_DECL(1,1,1)));
                        m_alpha[amrlev][mglev].define(ba, m_dmap[amrlev][mglev], 1, 0);
                        m_alpha[amrlev][mglev].setVal(0.0);
                    }
            }
        }
    }

    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        averageDownCoeffsSameAmrLevel(amrlev);
        averageDownCoeffsToCoarseAmrLevel(amrlev);
    }

    averageDownCoeffsSameAmrLevel(0);

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        if (m_use_harmonic_average) {
            int mglev = 0;
            FillBoundaryCoeff(*m_sigma[amrlev][mglev][0], m_geom[amrlev][mglev]);
            FillBoundaryCoeff(*m_sigmacross[amrlev][mglev][0], m_geom[amrlev][mglev]);
            for (mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    if (m_sigma[amrlev][mglev][idim]) {
                        FillBoundaryCoeff(*m_sigma[amrlev][mglev][idim], m_geom[amrlev][mglev]);
                    }

                    if (m_sigmacross[amrlev][mglev][idim]) {
                        FillBoundaryCoeff(*m_sigmacross[amrlev][mglev][idim], m_geom[amrlev][mglev]);
                    }
                }
            }
        } else {
            int idim = 0;
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                if (m_sigma[amrlev][mglev][idim]) {
                    FillBoundaryCoeff(*m_sigma[amrlev][mglev][idim], m_geom[amrlev][mglev]);
                }

                if (m_sigmacross[amrlev][mglev][idim]) {
                    FillBoundaryCoeff(*m_sigmacross[amrlev][mglev][idim], m_geom[amrlev][mglev]);
                }
            }
        }
    }
}

void
MLNodeHelmDualLinVel::averageDownCoeffsToCoarseAmrLevel (int flev)
{
    const int mglev = 0;
    const int idim = 0;  // other dimensions are just aliases
    amrex::average_down(*m_sigma[flev][mglev][idim], *m_sigma[flev-1][mglev][idim], 0, 1,
                        m_amr_ref_ratio[flev-1]);
    amrex::average_down(*m_sigmacross[flev][mglev][idim], *m_sigmacross[flev-1][mglev][idim], 0, 1,
                        m_amr_ref_ratio[flev-1]);

    // NOTE: The following line has not been tested yet!
    amrex::average_down_nodal(m_alpha[flev][mglev], m_alpha[flev-1][mglev],
                        amrex::IntVect(m_amr_ref_ratio[flev-1]));
}

void
MLNodeHelmDualLinVel::averageDownCoeffsSameAmrLevel (int amrlev)
{
    const int nsigma = (m_use_harmonic_average) ? AMREX_SPACEDIM : 1;

    for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
    {
        for (int idim = 0; idim < nsigma; ++idim)
        {
            const MultiFab& fine = *m_sigma[amrlev][mglev-1][idim];
            MultiFab& crse = *m_sigma[amrlev][mglev][idim];
            bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
            MultiFab cfine;
            if (need_parallel_copy) {
                const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
                cfine.define(ba, fine.DistributionMap(), 1, 0);
            }

            MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> const& cfab = pcrse->array(mfi);
                Array4<Real const> const& ffab = fine.const_array(mfi);
                if (idim == 0) {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                    {
                        mlndhelm_avgdown_coeff_x(i,j,k,cfab,ffab);
                    });
                } else if (idim == 1) {
#if (AMREX_SPACEDIM >= 2)
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                    {
                        mlndhelm_avgdown_coeff_y(i,j,k,cfab,ffab);
                    });
#endif
                } else {
#if (AMREX_SPACEDIM == 3)
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                    {
                        mlndhelm_avgdown_coeff_z(i,j,k,cfab,ffab);
                    });
#endif
                }
            }

            if (need_parallel_copy) {
                crse.ParallelCopy(cfine);
            }
        }
        for (int idim = 0; idim < nsigma; ++idim)
        {
            const MultiFab& fine = *m_sigmacross[amrlev][mglev-1][idim];
            MultiFab& crse = *m_sigmacross[amrlev][mglev][idim];
            bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
            MultiFab cfine;
            if (need_parallel_copy) {
                const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
                cfine.define(ba, fine.DistributionMap(), 1, 0);
            }

            MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> const& cfab = pcrse->array(mfi);
                Array4<Real const> const& ffab = fine.const_array(mfi);
                if (idim == 0) {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                    {
                        mlndhelm_avgdown_coeff_x(i,j,k,cfab,ffab);
                    });
                } else if (idim == 1) {
#if (AMREX_SPACEDIM >= 2)
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                    {
                        mlndhelm_avgdown_coeff_y(i,j,k,cfab,ffab);
                    });
#endif
                } else {
#if (AMREX_SPACEDIM == 3)
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                    {
                        mlndhelm_avgdown_coeff_z(i,j,k,cfab,ffab);
                    });
#endif
                }
            }

            if (need_parallel_copy) {
                crse.ParallelCopy(cfine);
            }
        }
        const MultiFab& fine = m_alpha[amrlev][mglev-1];
        MultiFab& crse = m_alpha[amrlev][mglev];
        average_down_nodal(fine, crse, IntVect(2));
    }
}

void
MLNodeHelmDualLinVel::FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom)
{
    BL_PROFILE("MLNodeHelmDualLinVel::FillBoundaryCoeff()");

    sigma.FillBoundary(geom.periodicity());

        const Box& domain = geom.Domain();
        const auto lobc = LoBC();
        const auto hibc = HiBC();

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sigma, mfi_info); mfi.isValid(); ++mfi)
        {
            Array4<Real> const& sfab = sigma.array(mfi);
            mlndhelm_fillbc_cc<Real>(mfi.validbox(),sfab,domain,lobc,hibc);
        }
}

void
MLNodeHelmDualLinVel::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
    if (!m_masks_built) buildMasks();

    const iMultiFab& cfmask = *m_nd_fine_mask[amrlev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(resmsk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<int> const& rmsk = resmsk.array(mfi);
        Array4<int const> const& fmsk = cfmask.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
        {
            if (fmsk(i,j,k) == crse_fine_node) rmsk(i,j,k) = 1;
        });
    }
}

void
MLNodeHelmDualLinVel::prepareForSolve ()
{
    BL_PROFILE("MLNodeHelmDualLinVel::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    buildMasks();

    averageDownCoeffs();
}

void
MLNodeHelmDualLinVel::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLNodeHelmDualLinVel::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
    }

    MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][cmglev-1];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> cfab = pcrse->array(mfi);
        Array4<Real const> const& ffab = fine.const_array(mfi);
        Array4<int const> const& mfab = dmsk.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndhelm_restriction(i,j,k,cfab,ffab,mfab);
            });
    }

    if (need_parallel_copy) {
        crse.ParallelCopy(cfine);
    }
}

void
MLNodeHelmDualLinVel::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
    BL_PROFILE("MLNodeHelmDualLinVel::interpolation()");

    const auto& sigma = m_sigma[amrlev][fmglev];

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    const MultiFab* cmf = &crse;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
        cfine.ParallelCopy(crse);
        cmf = &cfine;
    }

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][fmglev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const& ffab = fine.array(mfi);
        Array4<Real const> const& cfab = cmf->const_array(mfi);
        Array4<int const> const& mfab = dmsk.const_array(mfi);
        if (m_use_harmonic_average && fmglev > 0)
        {
            AMREX_D_TERM(Array4<Real const> const& sxfab = sigma[0]->const_array(mfi);,
                         Array4<Real const> const& syfab = sigma[1]->const_array(mfi);,
                         Array4<Real const> const& szfab = sigma[2]->const_array(mfi););
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndhelm_interpadd_ha(i,j,k,ffab,cfab,AMREX_D_DECL(sxfab,syfab,szfab),mfab);
            });
        }
        else
        {
            Array4<Real const> const& sfab = sigma[0]->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndhelm_interpadd_aa(i,j,k,ffab,cfab,sfab,mfab);
            });
        }
    }
}

void
MLNodeHelmDualLinVel::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
    const auto& amrrr = AMRRefRatio(camrlev);
    amrex::average_down(fine_sol, crse_sol, 0, 1, amrrr);

    if (isSingular(0))
    {
        MultiFab frhs(fine_rhs.boxArray(), fine_rhs.DistributionMap(), 1, 1);
        MultiFab::Copy(frhs, fine_rhs, 0, 0, 1, 0);
        restrictInteriorNodes(camrlev, crse_rhs, frhs);
    }
}

void
MLNodeHelmDualLinVel::restrictInteriorNodes (int camrlev, MultiFab& crhs, MultiFab& a_frhs) const
{
    const BoxArray& fba = a_frhs.boxArray();
    const DistributionMapping& fdm = a_frhs.DistributionMap();

    MultiFab* frhs = nullptr;
    std::unique_ptr<MultiFab> mf;
    if (a_frhs.nGrow() == 1)
    {
        frhs = &a_frhs;
    }
    else
    {
        mf.reset(new MultiFab(fba, fdm, 1, 1));
        frhs = mf.get();
        MultiFab::Copy(*frhs, a_frhs, 0, 0, 1, 0);
    }

    const Geometry& cgeom = m_geom[camrlev  ][0];

    const iMultiFab& fdmsk = *m_dirichlet_mask[camrlev+1][0];

    MultiFab cfine(amrex::coarsen(fba, 2), fdm, 1, 0);

    applyBC(camrlev+1, 0, *frhs, BCMode::Inhomogeneous, StateMode::Solution);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cfine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& cfab = cfine.array(mfi);
        Array4<Real const> const& ffab = frhs->const_array(mfi);
        Array4<int const> const& mfab = fdmsk.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndhelm_restriction(i,j,k,cfab,ffab,mfab);
            });
    }

    MultiFab tmp_crhs(crhs.boxArray(), crhs.DistributionMap(), 1, 0);
    tmp_crhs.setVal(0.0);
    tmp_crhs.ParallelCopy(cfine, cgeom.periodicity());

    const iMultiFab& c_nd_mask = *m_nd_fine_mask[camrlev];
    const auto& has_fine_bndry = *m_has_fine_bndry[camrlev];

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(crhs, mfi_info); mfi.isValid(); ++mfi)
    {
        if (has_fine_bndry[mfi])
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& dfab = crhs.array(mfi);
            Array4<Real const> const& sfab = tmp_crhs.const_array(mfi);
            Array4<int const> const& mfab = c_nd_mask.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                if (mfab(i,j,k) == fine_node) dfab(i,j,k) = sfab(i,j,k);
            });
        }
    }
}

void
MLNodeHelmDualLinVel::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLNodeHelmDualLinVel::Fapply()");

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& sigmacross = m_sigmacross[amrlev][mglev];
    const auto& alpha = m_alpha[amrlev][mglev];
    const auto dxinvarr = m_geom[amrlev][mglev].InvCellSizeArray();

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real const> const& xarr = in.const_array(mfi);
        Array4<Real> const& yarr = out.array(mfi);
        Array4<int const> const& dmskarr = dmsk.const_array(mfi);

        if (m_use_harmonic_average && mglev > 0)
        {
            AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                         Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                         Array4<Real const> const& szarr = sigma[2]->const_array(mfi););
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                mlndhelm_adotx_ha(i,j,k,yarr,xarr,AMREX_D_DECL(sxarr,syarr,szarr), dmskarr,
                                 dxinvarr);
            });
        }
        else
        {
            Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
            Array4<Real const> const& scarr = sigmacross[0]->const_array(mfi);
            Array4<Real const> const& aarr = alpha.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                mlndhelm_adotx_aa(i,j,k,yarr,xarr,sarr,scarr,aarr,dmskarr,
                                 dxinvarr);
            });
        }
    }
}

void
MLNodeHelmDualLinVel::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
    BL_PROFILE("MLNodeHelmDualLinVel::Fsmooth()");

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& sigmacross = m_sigmacross[amrlev][mglev];
    const auto& alpha = m_alpha[amrlev][mglev];
    const auto dxinvarr = m_geom[amrlev][mglev].InvCellSizeArray();

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        constexpr int nsweeps = 4;
        MultiFab Ax(sol.boxArray(), sol.DistributionMap(), 1, 0);

        if (m_use_harmonic_average && mglev > 0)
        {
            for (MFIter mfi(sol); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                             Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                             Array4<Real const> const& szarr = sigma[2]->const_array(mfi););
                Array4<Real> const& solarr = sol.array(mfi);
                Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                Array4<int const> const& dmskarr = dmsk.const_array(mfi);
                Array4<Real> const& Axarr = Ax.array(mfi);

                for (int ns = 0; ns < nsweeps; ++ns) {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        mlndhelm_adotx_ha(i,j,k,Axarr,solarr,AMREX_D_DECL(sxarr,syarr,szarr), dmskarr,
                                         dxinvarr);
                    });
                    AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
                    {
                        mlndhelm_jacobi_ha (tbx, solarr, Axarr, rhsarr, AMREX_D_DECL(sxarr,syarr,szarr),
                                           dmskarr, dxinvarr);
                    });
                }
            }
        }
        else
        {
            for (MFIter mfi(sol); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
                Array4<Real const> const& scarr = sigmacross[0]->const_array(mfi);
                Array4<Real const> const& aarr = alpha.const_array(mfi);
                Array4<Real> const& solarr = sol.array(mfi);
                Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                Array4<int const> const& dmskarr = dmsk.const_array(mfi);
                Array4<Real> const& Axarr = Ax.array(mfi);

                for (int ns = 0; ns < nsweeps; ++ns) {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        mlndhelm_adotx_aa(i,j,k,Axarr,solarr,sarr,scarr, aarr,dmskarr,
                                         dxinvarr);
                    });
                    AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
                    {
                        mlndhelm_jacobi_aa (tbx, solarr, Axarr, rhsarr, sarr,
                                           scarr, aarr, dmskarr, dxinvarr);
                    });
                }
            }
        }

        if (nsweeps > 1) nodalSync(amrlev, mglev, sol);
    }
    else // cpu
#endif
    {
        constexpr int nsweeps = 2;
        if (m_use_gauss_seidel)
        {
            if (m_use_harmonic_average && mglev > 0)
            {
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.validbox();
                    AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                                 Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                                 Array4<Real const> const& szarr = sigma[2]->const_array(mfi););
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    for (int ns = 0; ns < nsweeps; ++ns) {
                        mlndhelm_gauss_seidel_ha(bx, solarr, rhsarr,
                                                AMREX_D_DECL(sxarr,syarr,szarr),
                                                dmskarr, dxinvarr);
                    }
                }
            }
            else
            {
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.validbox();
                    Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
                    Array4<Real const> const& scarr = sigmacross[0]->const_array(mfi);
                    Array4<Real const> const& aarr = alpha.const_array(mfi);
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    for (int ns = 0; ns < nsweeps; ++ns) {
                        mlndhelm_gauss_seidel_aa(bx, solarr, rhsarr,
                                                sarr, scarr, aarr, dmskarr, dxinvarr);
                    }
                }
            }

            nodalSync(amrlev, mglev, sol);
        }
        else
        {
            MultiFab Ax(sol.boxArray(), sol.DistributionMap(), 1, 0);
            Fapply(amrlev, mglev, Ax, sol);

            if (m_use_harmonic_average && mglev > 0)
            {
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol,true); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                                 Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                                 Array4<Real const> const& szarr = sigma[2]->const_array(mfi););
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& Axarr = Ax.const_array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    mlndhelm_jacobi_ha (bx, solarr, Axarr, rhsarr, AMREX_D_DECL(sxarr,syarr,szarr),
                                       dmskarr, dxinvarr);
                }
            }
            else
            {
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol,true); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
                    Array4<Real const> const& scarr = sigmacross[0]->const_array(mfi);
                    Array4<Real const> const& aarr = alpha.const_array(mfi);
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& Axarr = Ax.const_array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    mlndhelm_jacobi_aa (bx, solarr, Axarr, rhsarr, sarr, scarr, aarr,
                                       dmskarr, dxinvarr);
                }
            }
        }
    }
}

void
MLNodeHelmDualLinVel::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    BL_PROFILE("MLNodeHelmDualLinVel::normalize()");

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& sigmacross = m_sigmacross[amrlev][mglev];
    const auto& alpha = m_alpha[amrlev][mglev];
    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& arr = mf.array(mfi);
        Array4<int const> const& dmskarr = dmsk.const_array(mfi);
        if (m_use_harmonic_average && mglev > 0)
        {
            AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                         Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                         Array4<Real const> const& szarr = sigma[2]->const_array(mfi););

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlndhelm_normalize_ha(tbx,arr,AMREX_D_DECL(sxarr,syarr,szarr),dmskarr,dxinv);
            });
        }
        else
        {
            Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
            Array4<Real const> const& scarr = sigmacross[0]->const_array(mfi);
            Array4<Real const> const& aarr = alpha.const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlndhelm_normalize_aa(tbx,arr,sarr,scarr,aarr,dmskarr,dxinv);
            });
        }
    }
}

// void
// MLNodeHelmDualLinVel::compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& a_phi,
//                                          const MultiFab& vold, const MultiFab* rhcc,
//                                          const BoxArray& fine_grids, const IntVect& ref_ratio)
// {
//     BL_PROFILE("MLNodeHelmDualLinVel::SyncResCrse()");
//
//     sync_resid.setVal(0.0);
//
//     const Geometry& geom = m_geom[0][0];
//     const DistributionMapping& dmap = m_dmap[0][0];
//     const BoxArray& ccba = m_grids[0][0];
//     const BoxArray& ndba = amrex::convert(ccba, IntVect::TheNodeVector());
//     const BoxArray& ccfba = amrex::convert(fine_grids, IntVect::TheZeroVector());
//     const auto lobc = LoBC();
//     const auto hibc = HiBC();
//
//     // cell-center, 1: coarse; 0: covered by fine
//     const int owner = 1;
//     const int nonowner = 0;
//     iMultiFab crse_cc_mask = amrex::makeFineMask(ccba, dmap, IntVect(1), ccfba, ref_ratio,
//                                                  geom.periodicity(), owner, nonowner);
//
//     const Box& ccdom = geom.Domain();
// #ifdef _OPENMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//     for (MFIter mfi(crse_cc_mask); mfi.isValid(); ++mfi)
//     {
//         Array4<int> const& fab = crse_cc_mask.array(mfi);
//         mlndhelm_fillbc_cc<int>(mfi.validbox(),fab,ccdom,lobc,hibc);
//     }
//
//     MultiFab phi(ndba, dmap, 1, 1);
// #ifdef _OPENMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//     for (MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
//     {
//         const Box& bx = mfi.tilebox();
//         const Box& gbx = mfi.growntilebox();
//         Array4<Real> const& fab = phi.array(mfi);
//         Array4<Real const> const& a_fab = a_phi.const_array(mfi);
//         Array4<int const> const& msk = crse_cc_mask.const_array(mfi);
//         AMREX_HOST_DEVICE_FOR_3D(gbx, i, j, k,
//         {
//             if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
//                 fab(i,j,k) = a_fab(i,j,k);
//                 mlndhelm_zero_fine(i,j,k,fab,msk,nonowner);
//             } else {
//                 fab(i,j,k) = 0.0;
//             }
//         });
//     }
//
//     const auto& nddom = amrex::surroundingNodes(ccdom);
//
//     const auto dxinv = geom.InvCellSizeArray();
//
//     const MultiFab& sigma_orig = *m_sigma[0][0][0];
//     const iMultiFab& dmsk = *m_dirichlet_mask[0][0];
//
//     MFItInfo mfi_info;
//     if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
// #ifdef _OPENMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//     {
//         FArrayBox rhs, u;
//         for (MFIter mfi(sync_resid,mfi_info); mfi.isValid(); ++mfi)
//         {
//             const Box& bx = mfi.tilebox();
//
//             auto typ = FabType::regular;
//             if (typ != FabType::covered)
//             {
//                 const Box& bxg1 = amrex::grow(bx,1);
//                 const Box& ccbxg1 = amrex::enclosedCells(bxg1);
//                 Array4<int const> const& cccmsk = crse_cc_mask.const_array(mfi);
//
//                 bool has_fine;
//                 if (Gpu::inLaunchRegion()) {
//                     AMREX_ASSERT(ccbxg1 == crse_cc_mask[mfi].box());
//                     has_fine = Reduce::AnyOf(ccbxg1,
//                                              [cccmsk,nonowner] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> bool
//                     {
//                         return cccmsk(i,j,k) == nonowner;
//                     });
//                 } else {
//                     has_fine = mlndhelm_any_fine_sync_cells(ccbxg1,cccmsk,nonowner);
//                 }
//
//                 if (has_fine)
//                 {
//                     const Box& ccvbx = amrex::enclosedCells(mfi.validbox());
//
//                     u.resize(ccbxg1, AMREX_SPACEDIM);
//                     Elixir ueli = u.elixir();
//                     Array4<Real> const& uarr = u.array();
//
//                     Box b = ccbxg1 & ccvbx;
//                     for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
//                     {
//                         if (m_lobc[0][idim] == LinOpBCType::inflow)
//                         {
//                             if (b.smallEnd(idim) == ccdom.smallEnd(idim)) {
//                                 b.growLo(idim, 1);
//                             }
//                         }
//                         if (m_hibc[0][idim] == LinOpBCType::inflow)
//                         {
//                             if (b.bigEnd(idim) == ccdom.bigEnd(idim)) {
//                                 b.growHi(idim, 1);
//                             }
//                         }
//                     }
//
//                     Array4<Real const> const& voarr = vold.const_array(mfi);
//                     AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
//                     {
// 		        if (b.contains(IntVect(AMREX_D_DECL(i,j,k))) and cccmsk(i,j,k)){
//                             AMREX_D_TERM(uarr(i,j,k,0) = voarr(i,j,k,0);,
//                                          uarr(i,j,k,1) = voarr(i,j,k,1);,
//                                          uarr(i,j,k,2) = voarr(i,j,k,2););
//                         } else {
//                             AMREX_D_TERM(uarr(i,j,k,0) = 0.0;,
//                                          uarr(i,j,k,1) = 0.0;,
//                                          uarr(i,j,k,2) = 0.0;);
//                         }
//                     });
//
//                     rhs.resize(bx);
//                     Elixir rhseli = rhs.elixir();
//                     Array4<Real> const& rhsarr = rhs.array();
//                     Array4<int const> const& dmskarr = dmsk.const_array(mfi);
//
//                     {
//                         AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
//                         {
//                             mlndhelm_divu(i,j,k,rhsarr,uarr,dmskarr,dxinv);
//                         });
//                     }
//
//                     if (rhcc)
//                     {
//                         Array4<Real> rhccarr = uarr;
//                         Array4<Real const> const& rhccarr_orig = rhcc->const_array(mfi);
//                         const Box& b2 = ccbxg1 & ccvbx;
//                         AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
//                         {
//  			    if (b2.contains(IntVect(AMREX_D_DECL(i,j,k))) and cccmsk(i,j,k)){
//                                 rhccarr(i,j,k) = rhccarr_orig(i,j,k);
//                             } else {
//                                 rhccarr(i,j,k) = 0.0;
//                             }
//                         });
//
//                         {
//                             AMREX_HOST_DEVICE_FOR_3D (bx, i, j, k,
//                             {
//                                 Real rhs2 = mlndhelm_rhcc(i,j,k,rhccarr,dmskarr);
//                                 rhsarr(i,j,k) += rhs2;
//                             });
//                         }
//                     }
//
//                     Array4<Real> const& sync_resid_a = sync_resid.array(mfi);
//                     Array4<Real const> const& phiarr = phi.const_array(mfi);
//                     Array4<Real const> const& sigmaarr_orig = sigma_orig.const_array(mfi);
//                     {
//                         Array4<Real> sigmaarr = uarr;
//                         const Box& ibx = ccbxg1 & amrex::enclosedCells(mfi.validbox());
//                         AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
//                         {
//                             if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k))) and cccmsk(i,j,k)) {
//                                 sigmaarr(i,j,k) = sigmaarr_orig(i,j,k);
//                             } else {
//                                 sigmaarr(i,j,k) = 0.0;
//                             }
//                         });
//
//                         AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
//                         {
//                             mlndhelm_adotx_aa(i, j, k, sync_resid_a, phiarr, sigmaarr, dmskarr,
//                                              dxinv);
//                             mlndhelm_crse_resid(i, j, k, sync_resid_a, rhsarr, cccmsk, nddom, lobc, hibc);
//                         });
//                     }
//                 }
//             }
//         }
//     }
// }
//
// void
// MLNodeHelmDualLinVel::compSyncResidualFine (MultiFab& sync_resid, const MultiFab& phi, const MultiFab& vold,
//                                        const MultiFab* rhcc)
// {
//     BL_PROFILE("MLNodeHelmDualLinVel::SyncResFine()");
//
//     const MultiFab& sigma_orig = *m_sigma[0][0][0];
//     const iMultiFab& dmsk = *m_dirichlet_mask[0][0];
//
//     const Geometry& geom = m_geom[0][0];
//     const Box& ccdom = geom.Domain();
//     const auto dxinv = geom.InvCellSizeArray();
//
//     MFItInfo mfi_info;
//     if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
// #ifdef _OPENMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//     {
//         FArrayBox rhs, u;
//         IArrayBox tmpmask;
//         for (MFIter mfi(sync_resid,mfi_info); mfi.isValid(); ++mfi)
//         {
//             const Box& bx = mfi.tilebox();
//
//             auto typ = FabType::regular;
//             if (typ != FabType::covered)
//             {
//                 const Box& gbx = mfi.growntilebox();
//                 const Box& vbx = mfi.validbox();
//                 const Box& ccvbx = amrex::enclosedCells(vbx);
//                 const Box& bxg1 = amrex::grow(bx,1);
//                 const Box& ccbxg1 = amrex::enclosedCells(bxg1);
//
//                 u.resize(ccbxg1, AMREX_SPACEDIM);
//                 Elixir ueli = u.elixir();
//                 Array4<Real> const& uarr = u.array();
//
//                 Box ovlp = ccvbx & ccbxg1;
//                 for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
//                 {
//                     if (m_lobc[0][idim] == LinOpBCType::inflow)
//                     {
//                         if (ovlp.smallEnd(idim) == ccdom.smallEnd(idim)) {
//                             ovlp.growLo(idim, 1);
//                         }
//                     }
//                     if (m_hibc[0][idim] == LinOpBCType::inflow)
//                     {
//                         if (ovlp.bigEnd(idim) == ccdom.bigEnd(idim)) {
//                             ovlp.growHi(idim, 1);
//                         }
//                     }
//                 }
//
//                 Array4<Real const> const& voarr = vold.const_array(mfi);
//                 AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
//                 {
//                     if (ovlp.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
//                         AMREX_D_TERM(uarr(i,j,k,0) = voarr(i,j,k,0);,
//                                      uarr(i,j,k,1) = voarr(i,j,k,1);,
//                                      uarr(i,j,k,2) = voarr(i,j,k,2););
//                     } else {
//                         AMREX_D_TERM(uarr(i,j,k,0) = 0.0;,
//                                      uarr(i,j,k,1) = 0.0;,
//                                      uarr(i,j,k,2) = 0.0;);
//                     }
//                 });
//
//                 tmpmask.resize(bx);
//                 Elixir tmeli = tmpmask.elixir();
//                 Array4<int> const& tmpmaskarr = tmpmask.array();
//                 Array4<int const> const& dmskarr = dmsk.const_array(mfi);
//                 AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
//                 {
//                     tmpmaskarr(i,j,k) = 1-dmskarr(i,j,k);
//                 });
//
//                 rhs.resize(bx);
//                 Elixir rhseli = rhs.elixir();
//                 Array4<Real> const& rhsarr = rhs.array();
//
//                {
//                     AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
//                     {
//                         mlndhelm_divu(i,j,k,rhsarr,uarr,tmpmaskarr,dxinv);
//                     });
//                 }
//
//                 if (rhcc)
//                 {
//                     Array4<Real> rhccarr = uarr;
//                     Array4<Real const> const& rhccarr_orig = rhcc->const_array(mfi);
//                     const Box& ovlp3 = ccvbx & ccbxg1;
//                     AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
//                     {
//                         if (ovlp3.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
//                             rhccarr(i,j,k) = rhccarr_orig(i,j,k);
//                         } else {
//                             rhccarr(i,j,k) = 0.0;
//                         }
//                     });
//                     {
//                         AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
//                         {
//                             Real rhs2 = mlndhelm_rhcc(i,j,k,rhccarr,tmpmaskarr);
//                             rhsarr(i,j,k) += rhs2;
//                         });
//                     }
//                 }
//
//                 Array4<Real> const& sync_resid_a = sync_resid.array(mfi);
//                 Array4<Real const> const& phiarr = phi.const_array(mfi);
//                 Array4<Real const> const& sigmaarr_orig = sigma_orig.const_array(mfi);
//                 {
//                     Array4<Real> sigmaarr = uarr;
//                     const Box& ovlp2 = ccvbx & ccbxg1;
//                     AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
//                     {
//                         if (ovlp2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
//                             sigmaarr(i,j,k) = sigmaarr_orig(i,j,k);
//                         } else {
//                             sigmaarr(i,j,k) = 0.0;
//                         }
//                     });
//
//                     AMREX_HOST_DEVICE_FOR_3D(gbx, i, j, k,
//                     {
//                         if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
//                             mlndhelm_adotx_aa(i,j,k, sync_resid_a, phiarr, sigmaarr, tmpmaskarr,
//                                              dxinv);
//                             sync_resid_a(i,j,k) = rhsarr(i,j,k) - sync_resid_a(i,j,k);
//                         } else {
//                             sync_resid_a(i,j,k) = 0.0;
//                         }
//                     });
//                 }
//             }
//
//             // Do not impose neumann bc here because how SyncRegister works.
//         }
//     }
// }

void
MLNodeHelmDualLinVel::reflux (int crse_amrlev,
                         MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                         MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const
{
    BL_PROFILE("MLNodeHelmDualLinVel::reflux()");

    const Geometry& cgeom = m_geom[crse_amrlev  ][0];
    const Geometry& fgeom = m_geom[crse_amrlev+1][0];
    const auto cdxinv = cgeom.InvCellSizeArray();
    const auto fdxinv = fgeom.InvCellSizeArray();
    const Box& c_cc_domain = cgeom.Domain();
    const Box& c_nd_domain = amrex::surroundingNodes(c_cc_domain);

    const BoxArray& fba = fine_sol.boxArray();
    const DistributionMapping& fdm = fine_sol.DistributionMap();

    const iMultiFab& fdmsk = *m_dirichlet_mask[crse_amrlev+1][0];

    MultiFab fine_res_for_coarse(amrex::coarsen(fba, 2), fdm, 1, 0);

    applyBC(crse_amrlev+1, 0, fine_res, BCMode::Inhomogeneous, StateMode::Solution);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine_res_for_coarse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& cfab = fine_res_for_coarse.array(mfi);
        Array4<Real const> const& ffab = fine_res.const_array(mfi);
        Array4<int const> const& mfab = fdmsk.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
            mlndhelm_restriction(i,j,k,cfab,ffab,mfab);
        });
    }
    res.ParallelCopy(fine_res_for_coarse, cgeom.periodicity());

    MultiFab fine_contrib(amrex::coarsen(fba, 2), fdm, 1, 0);
    fine_contrib.setVal(0.0);

    const auto& fsigma = *m_sigma[crse_amrlev+1][0][0];

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox sigfab;
        FArrayBox Axfab;
        for (MFIter mfi(fine_contrib,mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& cvbx = mfi.validbox();
            const Box& fvbx = amrex::refine(cvbx,2);
            const Box& cbx = mfi.tilebox();
            const Box& fbx = amrex::refine(cbx,2);

            const Box& cc_fbx = amrex::enclosedCells(fbx);
            const Box& cc_fvbx = amrex::enclosedCells(fvbx);
            const Box& bx_sig = amrex::grow(cc_fbx,2) & amrex::grow(cc_fvbx,1);
            const Box& b = bx_sig & cc_fvbx;

            sigfab.resize(bx_sig, 1);
            Elixir sigeli = sigfab.elixir();
            Array4<Real> const& sigarr = sigfab.array();
            Array4<Real const> const& sigarr_orig = fsigma.const_array(mfi);
            AMREX_HOST_DEVICE_FOR_3D(bx_sig, i, j, k,
            {
                if (b.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    sigarr(i,j,k) = sigarr_orig(i,j,k);
                } else {
                    sigarr(i,j,k) = 0.0;
                }
            });

            const Box& bx_Ax = amrex::grow(fbx,1);
            const Box& b2 = bx_Ax & amrex::grow(fvbx,-1);
            Axfab.resize(bx_Ax);
            Elixir Axeli = Axfab.elixir();
            Array4<Real> const& Axarr = Axfab.array();
            Array4<Real const> const& rhsarr = fine_rhs.const_array(mfi);
            Array4<Real const> const& resarr = fine_res.const_array(mfi);
            Array4<Real const> const& solarr = fine_sol.const_array(mfi);
            AMREX_HOST_DEVICE_FOR_3D(bx_Ax, i, j, k,
            {
                if (b2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    Axarr(i,j,k) = rhsarr(i,j,k) - resarr(i,j,k);
                } else {
                    Axarr(i,j,k) = 0.0;
                }
                mlndhelm_res_fine_Ax(i,j,k, fvbx, Axarr, solarr, sigarr,
                                    fdxinv);
            });

            Array4<Real> const& farr = fine_contrib.array(mfi);
            Array4<int const> const& marr = fdmsk.const_array(mfi);
            AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
            {
                mlndhelm_res_fine_contrib(i,j,k,farr,Axarr,marr);
            });
        }
    }

    MultiFab fine_contrib_on_crse(crse_sol.boxArray(), crse_sol.DistributionMap(), 1, 0);
    fine_contrib_on_crse.setVal(0.0);
    fine_contrib_on_crse.ParallelAdd(fine_contrib, cgeom.periodicity());

    const iMultiFab& cdmsk = *m_dirichlet_mask[crse_amrlev][0];
    const auto& nd_mask     = m_nd_fine_mask[crse_amrlev];
    const auto& cc_mask     = m_cc_fine_mask[crse_amrlev];
    const auto& has_fine_bndry = m_has_fine_bndry[crse_amrlev];

    const auto lobc = LoBC();
    const auto hibc = HiBC();

    const auto& csigma = *m_sigma[crse_amrlev][0][0];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(res,mfi_info); mfi.isValid(); ++mfi)
    {
        if ((*has_fine_bndry)[mfi])
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& resarr = res.array(mfi);
            Array4<Real const> const& csolarr = crse_sol.const_array(mfi);
            Array4<Real const> const& crhsarr = crse_rhs.const_array(mfi);
            Array4<Real const> const& csigarr = csigma.const_array(mfi);
            Array4<int const> const& cdmskarr = cdmsk.const_array(mfi);
            Array4<int const> const& ndmskarr = nd_mask->const_array(mfi);
            Array4<int const> const& ccmskarr = cc_mask->const_array(mfi);
            Array4<Real const> const& fcocarr = fine_contrib_on_crse.const_array(mfi);

            AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
            {
                mlndhelm_res_cf_contrib(i,j,k,resarr,csolarr,crhsarr,csigarr,
                                       cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                       cdxinv,c_nd_domain,
                                       lobc,hibc);
            });
        }
    }
}

void
MLNodeHelmDualLinVel::checkPoint (std::string const& file_name) const
{
    if (ParallelContext::IOProcessorSub())
    {
        UtilCreateCleanDirectory(file_name, false);
        {
            std::string HeaderFileName(file_name+"/Header");
            std::ofstream HeaderFile;
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) {
                FileOpenFailed(HeaderFileName);
            }

            HeaderFile.precision(17);

            // MLLinop stuff
            HeaderFile << "verbose = " << verbose << "\n"
                       << "nlevs = " << NAMRLevels() << "\n"
                       << "do_agglomeration = " << info.do_agglomeration << "\n"
                       << "do_consolidation = " << info.do_consolidation << "\n"
                       << "agg_grid_size = " << info.agg_grid_size << "\n"
                       << "con_grid_size = " << info.con_grid_size << "\n"
                       << "has_metric_term = " << info.has_metric_term << "\n"
                       << "max_coarsening_level = " << info.max_coarsening_level << "\n";
#if (AMREX_SPACEDIM == 1)
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0]) << "\n";
#elif (AMREX_SPACEDIM == 2)
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0])
                       << " "       << static_cast<int>(m_lobc[0][1]) << "\n";
#else
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0])
                       << " "       << static_cast<int>(m_lobc[0][1])
                       << " "       << static_cast<int>(m_lobc[0][2]) << "\n";
#endif
#if (AMREX_SPACEDIM == 1)
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0]) << "\n";
#elif (AMREX_SPACEDIM == 2)
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0])
                       << " "       << static_cast<int>(m_hibc[0][1]) << "\n";
#else
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0])
                       << " "       << static_cast<int>(m_hibc[0][1])
                       << " "       << static_cast<int>(m_hibc[0][2]) << "\n";
#endif
            // m_coarse_data_for_bc: not used
            HeaderFile << "maxorder = " << getMaxOrder() << "\n";

            // MLNodeHelmDualLinVel stuff
            HeaderFile << "use_gauss_seidel = " << m_use_gauss_seidel << "\n";
            HeaderFile << "use_harmonic_average = " << m_use_harmonic_average << "\n";
            HeaderFile << "coarsen_strategy = " << static_cast<int>(m_coarsening_strategy) << "\n";
            // No level bc multifab
        }

        for (int ilev = 0; ilev < NAMRLevels(); ++ilev)
        {
            UtilCreateCleanDirectory(file_name+"/Level_"+std::to_string(ilev), false);
            std::string HeaderFileName(file_name+"/Level_"+std::to_string(ilev)+"/Header");
            std::ofstream HeaderFile;
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) {
                FileOpenFailed(HeaderFileName);
            }

            HeaderFile.precision(17);

            HeaderFile << Geom(ilev) << "\n";
            m_grids[ilev][0].writeOn(HeaderFile);  HeaderFile << "\n";
        }
    }

    ParallelContext::BarrierSub();

    for (int ilev = 0; ilev < NAMRLevels(); ++ilev)
    {
        VisMF::Write(*m_sigma[ilev][0][0], file_name+"/Level_"+std::to_string(ilev)+"/sigma");
    }
}

}
// clang-format on

#pragma GCC diagnostic pop
