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

#ifndef AMREX_ML_NODEHELM_DUAL_CSTVEL_2D_K_H_
#define AMREX_ML_NODEHELM_DUAL_CSTVEL_2D_K_H_

namespace amrex {

//
// masks
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_set_nodal_mask (int i, int j, int k, Array4<int> const& nmsk,
                             Array4<int const> const& cmsk) noexcept
{
    int s = cmsk(i-1,j-1,k) + cmsk(i  ,j-1,k)
        +   cmsk(i-1,j  ,k) + cmsk(i  ,j  ,k);
    if (s == 4*crse_cell) {
        nmsk(i,j,k) = crse_node;
    }
    else if (s == 4*fine_cell) {
        nmsk(i,j,k) = fine_node;
    } else {
        nmsk(i,j,k) = crse_fine_node;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_set_dirichlet_mask (Box const& bx, Array4<int> const& dmsk,
                                 Array4<int const> const& omsk, Box const& dom,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        dmsk(i,j,0) = (omsk(i-1,j-1,0) == 1 or omsk(i,j-1,0) == 1 or
                       omsk(i-1,j  ,0) == 1 or omsk(i,j  ,0) == 1);
    }}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if (bclo[0] == LinOpBCType::Dirichlet and lo.x == domlo.x) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,0) = 1;
        }
    }

    if (bchi[0] == LinOpBCType::Dirichlet and hi.x == domhi.x) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,0) = 1;
        }
    }

    if (bclo[1] == LinOpBCType::Dirichlet and lo.y == domlo.y) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,0) = 1;
        }
    }

    if (bchi[1] == LinOpBCType::Dirichlet and hi.y == domhi.y) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,0) = 1;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_set_dot_mask (Box const& bx, Array4<Real> const& dmsk,
                           Array4<int const> const& omsk, Box const& dom,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        dmsk(i,j,0) = static_cast<Real>(omsk(i,j,0));
    }}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if ((bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow)
        and lo.x == domlo.x)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,0) *= 0.5;
        }
    }

    if ((bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow)
        and hi.x == domhi.x)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,0) *= 0.5;
        }
    }

    if ((bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow)
        and lo.y == domlo.y)
    {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,0) *= 0.5;
        }
    }

    if ((bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow)
        and hi.y == domhi.y)
    {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,0) *= 0.5;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_zero_fine (int i, int j, int, Array4<Real> const& phi,
                        Array4<int const> const& msk, int fine_flag) noexcept
{
    // Testing if the node is covered by a fine level in computing
    // coarse sync residual
    if (msk(i-1,j-1,0) == fine_flag and
        msk(i  ,j-1,0) == fine_flag and
        msk(i-1,j  ,0) == fine_flag and
        msk(i  ,j  ,0) == fine_flag)
    {
        phi(i,j,0) = 0.0;
    }
}

//
// coeffs
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_avgdown_coeff_x (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    Real a = fine(2*i  ,2*j,k) + fine(2*i  ,2*j+1,k);
    Real b = fine(2*i+1,2*j,k) + fine(2*i+1,2*j+1,k);
    crse(i,j,k) = a*b/(a+b);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_avgdown_coeff_y (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    Real a = fine(2*i,2*j  ,k) + fine(2*i+1,2*j  ,k);
    Real b = fine(2*i,2*j+1,k) + fine(2*i+1,2*j+1,k);
    crse(i,j,k) = a*b/(a+b);
}

//
// bc
//

template <typename T>
void mlndhelm_bc_doit (Box const& vbx, Array4<T> const& a, Box const& domain,
                      GpuArray<bool,AMREX_SPACEDIM> const& bflo,
                      GpuArray<bool,AMREX_SPACEDIM> const& bfhi) noexcept
{
    if (domain.strictly_contains(vbx)) return;

    int offset = domain.cellCentered() ? 0 : 1;

    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    Box const& sbox = amrex::grow(vbx,1);

    Box xlo_face = domain;
    xlo_face.setSmall(0,dlo.x-1);
    xlo_face.setBig  (0,dlo.x-1);
    xlo_face &= sbox;
    Box ylo_face = domain;
    ylo_face.setSmall(1,dlo.y-1);
    ylo_face.setBig  (1,dlo.y-1);
    ylo_face &= sbox;
    int xoffset = vbx.length(0)+1;
    int yoffset = vbx.length(1)+1;

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
    xlo_face, txbxlo,
    {
        auto lo = amrex::lbound(txbxlo);
        auto hi = amrex::ubound(txbxlo);
        if (lo.x == dlo.x-1 and bflo[0]) {
            for (int j = lo.y; j <= hi.y; ++j) {
                a(dlo.x-1,j,0) = a(dlo.x+offset,j,0);
            }
        }
        if (lo.x+xoffset == dhi.x+1 and bfhi[0]) {
            for (int j = lo.y; j <= hi.y; ++j) {
                a(dhi.x+1,j,0) = a(dhi.x-offset,j,0);
            }
        }
    },
    ylo_face, tybxlo,
    {
        auto lo = amrex::lbound(tybxlo);
        auto hi = amrex::ubound(tybxlo);
        if (lo.y == dlo.y-1 and bflo[1]) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,dlo.y-1,0) = a(i,dlo.y+offset,0);
            }
        }
        if (lo.y+yoffset == dhi.y+1 and bfhi[1]) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,dhi.y+1,0) = a(i,dhi.y-offset,0);
            }
        }
    });

    const auto lo = amrex::lbound(sbox);
    const auto hi = amrex::ubound(sbox);

    AMREX_HOST_DEVICE_FOR_1D ( 4, icorner,
    {
        switch (icorner) {
        case 0: {
            // xlo & ylo
            if (lo.x == dlo.x-1 and lo.y == dlo.y-1) {
                if (bflo[0]) {
                    a(dlo.x-1,dlo.y-1,0) = a(dlo.x+offset,dlo.y-1,0);
                } else if (bflo[1]) {
                    a(dlo.x-1,dlo.y-1,0) = a(dlo.x-1,dlo.y+offset,0);
                }
            }
            break;
        }
        case 1: {
            // xhi & ylo
            if (hi.x == dhi.x+1 and lo.y == dlo.y-1) {
                if (bfhi[0]) {
                    a(dhi.x+1,dlo.y-1,0) = a(dhi.x-offset,dlo.y-1,0);
                } else if (bflo[1]) {
                    a(dhi.x+1,dlo.y-1,0) = a(dhi.x+1,dlo.y+offset,0);
                }
            }
            break;
        }
        case 2: {
            // xlo & yhi
            if (lo.x == dlo.x-1 and hi.y == dhi.y+1) {
                if (bflo[0]) {
                    a(dlo.x-1,dhi.y+1,0) = a(dlo.x+offset,dhi.y+1,0);
                } else if (bfhi[1]) {
                    a(dlo.x-1,dhi.y+1,0) = a(dlo.x-1,dhi.y-offset,0);
                }
            }
            break;
        }
        case 3: {
            // xhi & yhi
            if (hi.x == dhi.x+1 and hi.y == dhi.y+1) {
                if (bfhi[0]) {
                    a(dhi.x+1,dhi.y+1,0) = a(dhi.x-offset,dhi.y+1,0);
                } else if (bfhi[1]) {
                    a(dhi.x+1,dhi.y+1,0) = a(dhi.x+1,dhi.y-offset,0);
                }
            }
            break;
        }
        default: {}
        }
    });
}

//
// operator
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_adotx_ha (int i, int j, int k, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sx, Array4<Real const> const& sy,
                       Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        y(i,j,k) = 0.0;
    } else {
        Real facx = (1./4.)*dxinv[0]*dxinv[0];
        Real facy = (1./4.)*dxinv[1]*dxinv[1];
        y(i,j,k) = x(i-1,j-1,k)*(facx*sx(i-1,j-1,k)+facy*sy(i-1,j-1,k))
               +   x(i+1,j-1,k)*(facx*sx(i  ,j-1,k)+facy*sy(i  ,j-1,k))
               +   x(i-1,j+1,k)*(facx*sx(i-1,j  ,k)+facy*sy(i-1,j  ,k))
               +   x(i+1,j+1,k)*(facx*sx(i  ,j  ,k)+facy*sy(i  ,j  ,k))
               +   x(i-1,j,k)*(  facx*(sx(i-1,j-1,k)+sx(i-1,j,k))
                               - facy*(sy(i-1,j-1,k)+sx(i-1,j,k)))
               +   x(i+1,j,k)*(  facx*(sx(i  ,j-1,k)+sx(i  ,j,k))
                               - facy*(sy(i  ,j-1,k)+sx(i  ,j,k)))
               +   x(i,j-1,k)*(- facx*(sx(i-1,j-1,k)+sx(i,j-1,k))
                               + facy*(sy(i-1,j-1,k)+sy(i,j-1,k)))
               +   x(i,j+1,k)*(- facx*(sx(i-1,j  ,k)+sx(i,j  ,k))
                               + facy*(sy(i-1,j  ,k)+sy(i,j  ,k)))
               +   x(i,j,k)*(-1.0)*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                                   +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_adotx_aa (int i, int j, int k, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sig, Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        y(i,j,k) = 0.0;
    } else {
        Real facx = (1.0/4.0)*dxinv[0]*dxinv[0];
        Real facy = (1.0/4.0)*dxinv[1]*dxinv[1];
        Real fxy = facx + facy;
        Real fxmy = facx - facy;
        Real fmxy = facy - facx;
        y(i,j,k) = x(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
               +   x(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
               +   x(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
               +   x(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
               +   x(i-1,j,k)*fxmy*(sig(i-1,j-1,k)+sig(i-1,j,k))
               +   x(i+1,j,k)*fxmy*(sig(i  ,j-1,k)+sig(i  ,j,k))
               +   x(i,j-1,k)*fmxy*(sig(i-1,j-1,k)+sig(i,j-1,k))
               +   x(i,j+1,k)*fmxy*(sig(i-1,j  ,k)+sig(i,j  ,k))
               +   x(i,j,k)*(-1.0)*fxy*
                      (sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_normalize_ha (Box const& bx, Array4<Real> const& x, Array4<Real const> const& sx,
                           Array4<Real const> const& sy, Array4<int const> const& msk,
                           GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/4.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/4.0)*dxinv[1]*dxinv[1];

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k)) {
            x(i,j,k) /= (-1.0)*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                               +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_normalize_aa (Box const& bx, Array4<Real> const& x, Array4<Real const> const& sig,
                           Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/4.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/4.0)*dxinv[1]*dxinv[1];
    Real fxy = facx + facy;

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k)) {
            x(i,j,k) /= (-1.0)*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_jacobi_ha (Box const& bx, Array4<Real> const& sol, Array4<Real const> const& Ax,
                        Array4<Real const> const& rhs, Array4<Real const> const& sx,
                        Array4<Real const> const& sy, Array4<int const> const& msk,
                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = -1.0 * (1.0/4.0)*dxinv[0]*dxinv[0];
    Real facy = -1.0 * (1.0/4.0)*dxinv[1]*dxinv[1];

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            sol(i,j,k) += (2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                +  facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_jacobi_aa (Box const& bx, Array4<Real> const& sol, Array4<Real const> const& Ax,
                        Array4<Real const> const& rhs, Array4<Real const> const& sig,
                        Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real fac = -1.0 * (1.0/4.0)*(dxinv[0]*dxinv[0] + dxinv[1]*dxinv[1]);

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            sol(i,j,k) += (2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (fac*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k)));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_gauss_seidel_ha (Box const& bx, Array4<Real> const& sol,
                              Array4<Real const> const& rhs, Array4<Real const> const& sx,
                              Array4<Real const> const& sy, Array4<int const> const& msk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/4.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/4.0)*dxinv[1]*dxinv[1];

    amrex::Loop(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            Real s0 = (-1.0)*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                             +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));

            Real Ax = sol(i-1,j-1,k)*(facx*sx(i-1,j-1,k)+facy*sy(i-1,j-1,k))
                    + sol(i+1,j-1,k)*(facx*sx(i  ,j-1,k)+facy*sy(i  ,j-1,k))
                    + sol(i-1,j+1,k)*(facx*sx(i-1,j  ,k)+facy*sy(i-1,j  ,k))
                    + sol(i+1,j+1,k)*(facx*sx(i  ,j  ,k)+facy*sy(i  ,j  ,k))
                    + sol(i-1,j,k)*(  facx*(sx(i-1,j-1,k)+sx(i-1,j,k))
                                    - facy*(sy(i-1,j-1,k)+sx(i-1,j,k)))
                    + sol(i+1,j,k)*(  facx*(sx(i  ,j-1,k)+sx(i  ,j,k))
                                    - facy*(sy(i  ,j-1,k)+sx(i  ,j,k)))
                    + sol(i,j-1,k)*(- facx*(sx(i-1,j-1,k)+sx(i,j-1,k))
                                    + facy*(sy(i-1,j-1,k)+sy(i,j-1,k)))
                    + sol(i,j+1,k)*(- facx*(sx(i-1,j  ,k)+sx(i,j  ,k))
                                    + facy*(sy(i-1,j  ,k)+sy(i,j  ,k)))
                    + sol(i,j,k)*s0;

            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_gauss_seidel_aa (Box const& bx, Array4<Real> const& sol,
                              Array4<Real const> const& rhs, Array4<Real const> const& sig,
                              Array4<int const> const& msk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/4.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/4.0)*dxinv[1]*dxinv[1];
    Real fxy = facx + facy;
    Real fxmy = facx - facy;
    Real fmxy = facy - facx;

    amrex::Loop(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            Real s0 = (-1.0)*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
            Real Ax =   sol(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
                      + sol(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
                      + sol(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
                      + sol(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
                      + sol(i-1,j,k)*fxmy*(sig(i-1,j-1,k)+sig(i-1,j,k))
                      + sol(i+1,j,k)*fxmy*(sig(i  ,j-1,k)+sig(i  ,j,k))
                      + sol(i,j-1,k)*fmxy*(sig(i-1,j-1,k)+sig(i,j-1,k))
                      + sol(i,j+1,k)*fmxy*(sig(i-1,j  ,k)+sig(i,j  ,k))
                      + sol(i,j,k)*s0;

            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

//
// restriction
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_restriction (int i, int j, int k, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk) noexcept
{
    int ii = i*2;
    int jj = j*2;
    int kk = 0;
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = 0.0;
    } else {
        crse(i,j,k) = (1./16.)*(fine(ii-1,jj-1,kk) + 2.*fine(ii  ,jj-1,kk) +    fine(ii+1,jj-1,kk)
                           + 2.*fine(ii-1,jj  ,kk) + 4.*fine(ii  ,jj  ,kk) + 2.*fine(ii+1,jj  ,kk)
                              + fine(ii-1,jj+1,kk) + 2.*fine(ii  ,jj+1,kk) +    fine(ii+1,jj+1,kk));
    }
}

//
// interpolation
//

namespace {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_x (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int ic, int jc) noexcept
    {
        Real w1 = sig(i-1,j-1,0) + sig(i-1,j,0);
        Real w2 = sig(i  ,j-1,0) + sig(i  ,j,0);
        return (w1*crse(ic,jc,0)+w2*crse(ic+1,jc,0))/(w1+w2);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_y (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int ic, int jc) noexcept
    {
        Real w1 = sig(i-1,j-1,0) + sig(i,j-1,0);
        Real w2 = sig(i-1,j  ,0) + sig(i,j  ,0);
        return (w1*crse(ic,jc,0)+w2*crse(ic,jc+1,0))/(w1+w2);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_xy (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int ic, int jc) noexcept
    {
        Real w1 = sig(i-1,j-1,0) + sig(i-1,j,0);
        Real w2 = sig(i  ,j-1,0) + sig(i  ,j,0);
        Real w3 = sig(i-1,j-1,0) + sig(i,j-1,0);
        Real w4 = sig(i-1,j  ,0) + sig(i,j  ,0);
        return (w1 * aa_interp_line_y(crse,sig,i-1,j  ,ic  ,jc  ) +
                w2 * aa_interp_line_y(crse,sig,i+1,j  ,ic+1,jc  ) +
                w3 * aa_interp_line_x(crse,sig,i  ,j-1,ic  ,jc  ) +
                w4 * aa_interp_line_x(crse,sig,i  ,j+1,ic  ,jc+1)) / (w1+w2+w3+w4);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_interpadd_aa (int i, int j, int, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sig,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,0)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        if (i_is_odd and j_is_odd) {
            // Node on a X-Y face
            fine(i,j,0) += aa_interp_face_xy(crse,sig,i,j,ic,jc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,0) += aa_interp_line_x(crse,sig,i,j,ic,jc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,0) += aa_interp_line_y(crse,sig,i,j,ic,jc);
        } else {
            // Node coincident with coarse node
            fine(i,j,0) += crse(ic,jc,0);
        }
    }
}

namespace {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ha_interp_face_xy (Array4<Real const> const& crse,
                            Array4<Real const> const& sigx, Array4<Real const> const& sigy,
                            int i, int j, int ic, int jc) noexcept
    {
        Real w1 = sigx(i-1,j-1,0) + sigx(i-1,j,0);
        Real w2 = sigx(i  ,j-1,0) + sigx(i  ,j,0);
        Real w3 = sigy(i-1,j-1,0) + sigy(i,j-1,0);
        Real w4 = sigy(i-1,j  ,0) + sigy(i,j  ,0);
        return (w1 * aa_interp_line_y(crse,sigy,i-1,j  ,ic  ,jc  ) +
                w2 * aa_interp_line_y(crse,sigy,i+1,j  ,ic+1,jc  ) +
                w3 * aa_interp_line_x(crse,sigx,i  ,j-1,ic  ,jc  ) +
                w4 * aa_interp_line_x(crse,sigx,i  ,j+1,ic  ,jc+1)) / (w1+w2+w3+w4);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_interpadd_ha (int i, int j, int,
                           Array4<Real> const& fine, Array4<Real const> const& crse,
                           Array4<Real const> const& sigx, Array4<Real const> const& sigy,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,0)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        if (i_is_odd and j_is_odd) {
            // Node on a X-Y face
            fine(i,j,0) += ha_interp_face_xy(crse,sigx,sigy,i,j,ic,jc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,0) += aa_interp_line_x(crse,sigx,i,j,ic,jc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,0) += aa_interp_line_y(crse,sigy,i,j,ic,jc);
        } else {
            // Node coincident with coarse node
            fine(i,j,0) += crse(ic,jc,0);
        }
    }
}

//
// rhs & u
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_divu (int i, int j, int k, Array4<Real> const& rhs, Array4<Real const> const& vel,
                   Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = 0.5*dxinv[0];
    Real facy = 0.5*dxinv[1];

    if (msk(i,j,k)) {
        rhs(i,j,k) = 0.0;
    } else {
        rhs(i,j,k) = facx*(-vel(i-1,j-1,k,0) + vel(i,j-1,k,0)
                           -vel(i-1,j  ,k,0) + vel(i,j  ,k,0))
                   + facy*(-vel(i-1,j-1,k,1) - vel(i,j-1,k,1)
                           +vel(i-1,j  ,k,1) + vel(i,j  ,k,1));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndhelm_rhcc (int i, int j, int k, Array4<Real const> const& rhcc,
                   Array4<int const> const& msk) noexcept
{
    Real r;
    if (msk(i,j,k)) {
        r = 0.0;
    } else {
        r = 0.25 * (rhcc(i-1,j-1,k)+rhcc(i,j-1,k)+rhcc(i-1,j,k)+rhcc(i,j,k));
    }
    return r;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_mknewu (int i, int j, int k, Array4<Real> const& u, Array4<Real const> const& p,
                     Array4<Real const> const& sig, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = 0.5*dxinv[0];
    Real facy = 0.5*dxinv[1];
    u(i,j,k,0) -= sig(i,j,k)*facx*(-p(i,j,k)+p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    u(i,j,k,1) -= sig(i,j,k)*facy*(-p(i,j,k)-p(i+1,j,k)+p(i,j+1,k)+p(i+1,j+1,k));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_divu_compute_fine_contrib (int i, int j, int, Box const& fvbx,
                                        Array4<Real> const& frh, Array4<Real const> const& vel,
                                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    IntVect iv(i,j);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv))
    {
        Real facx = 0.5_rt*dxinv[0];
        Real facy = 0.5_rt*dxinv[1];
        frh(i,j,0) = facx*(-vel(i-1,j-1,0,0)+vel(i,j-1,0,0)-vel(i-1,j,0,0)+vel(i,j,0,0))
            +        facy*(-vel(i-1,j-1,0,1)-vel(i,j-1,0,1)+vel(i-1,j,0,1)+vel(i,j,0,1));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_divu_add_fine_contrib (int i, int j, int /*k*/, Box const& fvbx,
                                    Array4<Real> const& rhs, Array4<Real const> const& frh,
                                    Array4<int const> const& msk) noexcept
{
    constexpr Real rfd = 0.25_rt;
    constexpr Real chip = 0.5_rt;
    constexpr Real chip2 = 0.25_rt;

    int ii = 2*i;
    int jj = 2*j;
    IntVect iv(ii,jj);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv) and msk(ii,jj,0))
    {
        rhs(i,j,0) +=
            rfd*(frh(ii,jj,0)
                 + chip*(frh(ii-1,jj,0)+frh(ii+1,jj,0)+frh(ii,jj-1,0)+frh(ii,jj+1,0))
                 + chip2*(frh(ii-1,jj-1,0)+frh(ii+1,jj-1,0)+frh(ii-1,jj+1,0)+frh(ii+1,jj+1,0)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_rhcc_fine_contrib (int i, int j, int, Box const& fvbx,
                                Array4<Real> const& rhs, Array4<Real const> const& cc,
                                Array4<int const> const& msk) noexcept
{
    constexpr Real w1 = 9._rt/64._rt;
    constexpr Real w2 = 3._rt/64._rt;
    constexpr Real w3 = 1._rt/64._rt;

    int ii = 2*i;
    int jj = 2*j;
    IntVect iv(ii,jj);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv) and msk(ii,jj,0))
    {
        rhs(i,j,0) += w1*(cc(ii-1,jj-1,0)+cc(ii  ,jj-1,0)+cc(ii-1,jj  ,0)+cc(ii  ,jj  ,0))
                    + w2*(cc(ii-2,jj-1,0)+cc(ii+1,jj-1,0)+cc(ii-2,jj  ,0)+cc(ii+1,jj  ,0)
                         +cc(ii-1,jj-2,0)+cc(ii  ,jj-2,0)+cc(ii-1,jj+1,0)+cc(ii  ,jj+1,0))
                    + w3*(cc(ii-2,jj-2,0)+cc(ii+1,jj-2,0)+cc(ii-2,jj+1,0)+cc(ii+1,jj+1,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_divu_cf_contrib (int i, int j, int, Array4<Real> const& rhs,
                              Array4<Real const> const& vel, Array4<Real const> const& fc,
                              Array4<Real const> const& rhcc, Array4<int const> const& dmsk,
                              Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi) noexcept
{
    if (!dmsk(i,j,0) and ndmsk(i,j,0) == crse_fine_node) {
        Real facx = 0.5_rt * dxinv[0];
        Real facy = 0.5_rt * dxinv[1];
        const auto ndlo = amrex::lbound(nddom);
        const auto ndhi = amrex::ubound(nddom);
        Real r = fc(i,j,0);
        if (rhcc) {
            r += 0.25_rt*( (1._rt-ccmsk(i-1,j-1,0)) * rhcc(i-1,j-1,0)
                         + (1._rt-ccmsk(i  ,j-1,0)) * rhcc(i  ,j-1,0)
                         + (1._rt-ccmsk(i-1,j  ,0)) * rhcc(i-1,j  ,0)
                         + (1._rt-ccmsk(i  ,j  ,0)) * rhcc(i  ,j  ,0));
        }
        r += (1._rt-ccmsk(i-1,j-1,0)) * (-facx*vel(i-1,j-1,0,0) - facy*vel(i-1,j-1,0,1))
           + (1._rt-ccmsk(i  ,j-1,0)) * ( facx*vel(i  ,j-1,0,0) - facy*vel(i  ,j-1,0,1))
           + (1._rt-ccmsk(i-1,j  ,0)) * (-facx*vel(i-1,j  ,0,0) + facy*vel(i-1,j  ,0,1))
           + (1._rt-ccmsk(i  ,j  ,0)) * ( facx*vel(i  ,j  ,0,0) + facy*vel(i  ,j  ,0,1));

        if (i == ndlo.x and ( bclo[0] == LinOpBCType::Neumann or
                              bclo[0] == LinOpBCType::inflow)) {
            r *= 2._rt;
        } else if (i== ndhi.x and ( bchi[0] == LinOpBCType::Neumann or
                                    bchi[0] == LinOpBCType::inflow)) {
            r *= 2._rt;
        }

        if (j == ndlo.y and ( bclo[1] == LinOpBCType::Neumann or
                              bclo[1] == LinOpBCType::inflow)) {
            r *= 2._rt;
        } else if (j == ndhi.y and ( bchi[1] == LinOpBCType::Neumann or
                                     bchi[1] == LinOpBCType::inflow)) {
            r *= 2._rt;
        }

        rhs(i,j,0) = r;
    }
}

//
// residual
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_crse_resid (int i, int j, int k, Array4<Real> const& resid,
                         Array4<Real const> const& rhs, Array4<int const> const& msk,
                         Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                         GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi) noexcept
{
    if ((msk(i-1,j-1,k  ) == 0 or
         msk(i  ,j-1,k  ) == 0 or
         msk(i-1,j  ,k  ) == 0 or
         msk(i  ,j  ,k  ) == 0) and
        (msk(i-1,j-1,k  ) == 0 or
         msk(i  ,j-1,k  ) == 0 or
         msk(i-1,j  ,k  ) == 0 or
         msk(i  ,j  ,k  ) == 0))
    {
        const auto ndlo = amrex::lbound(nddom);
        const auto ndhi = amrex::ubound(nddom);
        Real fac = 1.0_rt;
        if (i == ndlo.x and ( bclo[0] == LinOpBCType::Neumann or
                              bclo[0] == LinOpBCType::inflow)) {
            fac *= 2._rt;
        } else if (i== ndhi.x and ( bchi[0] == LinOpBCType::Neumann or
                                    bchi[0] == LinOpBCType::inflow)) {
            fac *= 2._rt;
        }

        if (j == ndlo.y and ( bclo[1] == LinOpBCType::Neumann or
                              bclo[1] == LinOpBCType::inflow)) {
            fac *= 2._rt;
        } else if (j == ndhi.y and ( bchi[1] == LinOpBCType::Neumann or
                                     bchi[1] == LinOpBCType::inflow)) {
            fac *= 2._rt;
        }

        resid(i,j,k) = (rhs(i,j,k) - resid(i,j,k)) * fac;
    } else {
        resid(i,j,k) = 0._rt;
    }
}

//
// sync residual
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_res_fine_Ax (int i, int j, int, Box const& fvbx, Array4<Real> const& Ax,
                          Array4<Real const> const& x, Array4<Real const> const& sig,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    IntVect iv(i,j);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv)) {
        Real facx = (1._rt/4._rt)*dxinv[0]*dxinv[0];
        Real facy = (1._rt/4._rt)*dxinv[1]*dxinv[1];
        Real fxy = facx + facy;
        Real fxmy = facx - facy;
        Real fmxy = facy - facx;
        Ax(i,j,0) = x(i-1,j-1,0)*fxy*sig(i-1,j-1,0)
            +       x(i+1,j-1,0)*fxy*sig(i  ,j-1,0)
            +       x(i-1,j+1,0)*fxy*sig(i-1,j  ,0)
            +       x(i+1,j+1,0)*fxy*sig(i  ,j  ,0)
            +       x(i-1,j,0)*fxmy*(sig(i-1,j-1,0)+sig(i-1,j  ,0))
            +       x(i+1,j,0)*fxmy*(sig(i  ,j-1,0)+sig(i  ,j  ,0))
            +       x(i,j-1,0)*fmxy*(sig(i-1,j-1,0)+sig(i  ,j-1,0))
            +       x(i,j+1,0)*fmxy*(sig(i-1,j  ,0)+sig(i  ,j  ,0))
            +       x(i,j,0)*(-1._rt)*fxy*(sig(i-1,j-1,0)+sig(i,j-1,0)+sig(i-1,j,0)+sig(i,j,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_res_fine_contrib (int i, int j, int k, Array4<Real> const& f,
                               Array4<Real const> const& Ax,
                               Array4<int const> const& msk) noexcept
{
    constexpr Real rfd = 0.25_rt;
    constexpr Real chip = 0.5_rt;
    constexpr Real chip2 = 0.25_rt;

    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    if (msk(ii,jj,kk)) {
        f(i,j,0) += rfd*(Ax(ii,jj,0)
                 + chip*(Ax(ii-1,jj,0)+Ax(ii+1,jj,0)+Ax(ii,jj-1,0)+Ax(ii,jj+1,0))
                 + chip2*(Ax(ii-1,jj-1,0)+Ax(ii+1,jj-1,0)+Ax(ii-1,jj+1,0)+Ax(ii+1,jj+1,0)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_res_cf_contrib (int i, int j, int, Array4<Real> const& res,
                             Array4<Real const> const& phi, Array4<Real const> const& rhs,
                             Array4<Real const> const& sig, Array4<int const> const& dmsk,
                             Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                             Array4<Real const> const& fc,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv, Box const& nddom,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    if (!dmsk(i,j,0) and ndmsk(i,j,0) == crse_fine_node) {
        Real facx = (1._rt/6._rt)*dxinv[0]*dxinv[0];
        Real facy = (1._rt/6._rt)*dxinv[1]*dxinv[1];

        Real Ax = 0._rt;
        if (ccmsk(i-1,j-1,0) == crse_cell) {
            Ax += sig(i-1,j-1,0)*(facx*(2._rt*(phi(i-1,j  ,0)-phi(i  ,j  ,0))
                                        +     (phi(i-1,j-1,0)-phi(i  ,j-1,0)))
                                + facy*(2._rt*(phi(i  ,j-1,0)-phi(i  ,j  ,0))
                                        +     (phi(i-1,j-1,0)-phi(i-1,j  ,0))));
        }
        if (ccmsk(i,j-1,0) == crse_cell) {
            Ax += sig(i,j-1,0)*(facx*(2._rt*(phi(i+1,j  ,0)-phi(i  ,j  ,0))
                                      +     (phi(i+1,j-1,0)-phi(i  ,j-1,0)))
                              + facy*(2._rt*(phi(i  ,j-1,0)-phi(i  ,j  ,0))
                                      +     (phi(i+1,j-1,0)-phi(i+1,j  ,0))));
        }
        if (ccmsk(i-1,j,0) == crse_cell) {
            Ax += sig(i-1,j,0)*(facx*(2._rt*(phi(i-1,j  ,0)-phi(i  ,j  ,0))
                                      +     (phi(i-1,j+1,0)-phi(i  ,j+1,0)))
                              + facy*(2._rt*(phi(i  ,j+1,0)-phi(i  ,j  ,0))
                                      +     (phi(i-1,j+1,0)-phi(i-1,j  ,0))));
        }
        if (ccmsk(i,j,0) == crse_cell) {
            Ax += sig(i,j,0)*(facx*(2._rt*(phi(i+1,j  ,0)-phi(i  ,j  ,0))
                                   +      (phi(i+1,j+1,0)-phi(i  ,j+1,0)))
                            + facy*(2._rt*(phi(i  ,j+1,0)-phi(i  ,j  ,0))
                                   +      (phi(i+1,j+1,0)-phi(i+1,j  ,0))));
        }

        Real Axf = fc(i,j,0);
        const auto ndlo = amrex::lbound(nddom);
        const auto ndhi = amrex::ubound(nddom);

        if (i == ndlo.x and (bclo[0] == LinOpBCType::Neumann or
                             bclo[0] == LinOpBCType::inflow)) {
            Axf *= 2._rt;
        } else if (i== ndhi.x and (bchi[0] == LinOpBCType::Neumann or
                                   bchi[0] == LinOpBCType::inflow)) {
            Axf *= 2._rt;
        }

        if (j == ndlo.y and (bclo[1] == LinOpBCType::Neumann or
                             bclo[1] == LinOpBCType::inflow)) {
            Axf *= 2._rt;
        } else if (j == ndhi.y and (bchi[1] == LinOpBCType::Neumann or
                                    bchi[1] == LinOpBCType::inflow)) {
            Axf *= 2._rt;
        }

        res(i,j,0) = rhs(i,j,0) - (Ax + Axf);
    }
}

}
#endif
