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
#ifndef AMREX_ML_NODEHELM_DUAL_CSTVEL_HPP
#define AMREX_ML_NODEHELM_DUAL_CSTVEL_HPP

#include <AMReX.H>
#include "fub/AMReX/MLMG/MLNodeHelmholtz.hpp"

namespace amrex {

class MLNodeHelmDualCstVel
    : public MLNodeHelmholtz
{
public:

    MLNodeHelmDualCstVel () noexcept {}
    MLNodeHelmDualCstVel (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info = LPInfo(),
                     const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    virtual ~MLNodeHelmDualCstVel ();

    MLNodeHelmDualCstVel (const MLNodeHelmDualCstVel&) = delete;
    MLNodeHelmDualCstVel (MLNodeHelmDualCstVel&&) = delete;
    MLNodeHelmDualCstVel& operator= (const MLNodeHelmDualCstVel&) = delete;
    MLNodeHelmDualCstVel& operator= (MLNodeHelmDualCstVel&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    virtual std::string name () const override { return std::string("MLNodeHelmDualCstVel"); }

    void setNormalizationThreshold (Real t) noexcept { m_normalization_threshold = t; }

    void setSigma (int amrlev, const MultiFab& a_sigma);

    void setSigmaCross (int amrlev, const MultiFab& a_sigmacross);

    void setAlpha (int amrlev, const MultiFab& a_alpha);

    void compDivergence (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel);

    void compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
                  const Vector<const MultiFab*>& rhnd,
                  const Vector<MultiFab*>& rhcc);

    void updateVelocity (const Vector<MultiFab*>& vel, const Vector<MultiFab const*>& sol) const;

//     void compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& phi,
//                                  const MultiFab& vold, const MultiFab* rhcc,
//                                  const BoxArray& fine_grids, const IntVect& ref_ratio);
//
//     void compSyncResidualFine (MultiFab& sync_resid, const MultiFab& phi, const MultiFab& vold,
//                                const MultiFab* rhcc);

    void setGaussSeidel (bool flag) noexcept { m_use_gauss_seidel = flag; }
    void setHarmonicAverage (bool flag) noexcept { m_use_harmonic_average = flag; }

    void setCoarseningStrategy (CoarseningStrategy cs) noexcept { m_coarsening_strategy = cs; }

    virtual BottomSolver getDefaultBottomSolver () const final override {
        return BottomSolver::bicgstab;
    }

    virtual void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final override;
    virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final override;
    virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final override;

    virtual void reflux (int crse_amrlev,
                         MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                         MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const final override;

    virtual void prepareForSolve () final override;
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final override;
    virtual void normalize (int amrlev, int mglev, MultiFab& mf) const final override;

    virtual void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final override;

    virtual void getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& ,
                            const Vector<MultiFab*>& ,
                            Location ) const final override {
        amrex::Abort("MLNodeHelmDualCstVel::getFluxes: How did we get here?");
    }
    virtual void getFluxes (const Vector<MultiFab*>& a_flux,
                            const Vector<MultiFab*>& a_sol) const final override;

    // virtual void unimposeNeumannBC (int , MultiFab& ) const final override {};

    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);
    void averageDownCoeffsSameAmrLevel (int amrlev);

    void restrictInteriorNodes (int camrlev, MultiFab& crhs, MultiFab& frhs) const;

    void FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom);

private:

    Vector<Vector<Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> > > m_sigma;
    Vector<Vector<Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> > > m_sigmacross;
    Vector<Vector<MultiFab> > m_alpha;

    Real m_normalization_threshold = 1.e-10;

    bool m_use_gauss_seidel = true;
    bool m_use_harmonic_average = false;

    virtual void checkPoint (std::string const& file_name) const final;
};

} // namespace amrex

#endif
// clang-format on
