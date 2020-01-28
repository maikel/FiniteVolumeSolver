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
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Notice: This file is copied and modified from AMReX-Codes
// (https://amrex-codes.github.io/).

#ifndef AMREX_ML_NODEHELM_DUAL_CSTVEL_1D_K_H_
#define AMREX_ML_NODEHELM_DUAL_CSTVEL_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_set_nodal_mask(int i, int j, int k, Array4<int> const& nmsk,
                        Array4<int const> const& cmsk) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void mlndhelm_set_dirichlet_mask(
    Box const& bx, Array4<int> const& dmsk, Array4<int const> const& omsk,
    Box const& dom, GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void mlndhelm_set_dot_mask(
    Box const& bx, Array4<Real> const& dmsk, Array4<int const> const& omsk,
    Box const& dom, GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_zero_fine(int i, int j, int k, Array4<Real> const& phi,
                   Array4<int const> const& msk, int fine_flag) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_avgdown_coeff_x(int i, int j, int k, Array4<Real> const& crse,
                         Array4<Real const> const& fine) noexcept {}

template <typename T>
inline void
mlndhelm_bc_doit(Box const& vbx, Array4<T> const& a, Box const& domain,
                 GpuArray<bool, AMREX_SPACEDIM> const& bflo,
                 GpuArray<bool, AMREX_SPACEDIM> const& bfhi) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_adotx_ha(int i, int j, int k, Array4<Real> const& y,
                  Array4<Real const> const& x, Array4<Real const> const& sx,
                  Array4<int const> const& msk,
                  GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_adotx_aa(int i, int j, int k, Array4<Real> const& y,
                  Array4<Real const> const& x, Array4<Real const> const& sig,
                  Array4<int const> const& msk,
                  GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_normalize_ha(Box const& bx, Array4<Real> const& x,
                      Array4<Real const> const& sx,
                      Array4<int const> const& msk,
                      GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_normalize_aa(Box const& bx, Array4<Real> const& x,
                      Array4<Real const> const& sig,
                      Array4<int const> const& msk,
                      GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_jacobi_ha(Box const& bx, Array4<Real> const& sol,
                   Array4<Real const> const& Ax, Array4<Real const> const& rhs,
                   Array4<Real const> const& sx, Array4<int const> const& msk,
                   GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_jacobi_aa(Box const& bx, Array4<Real> const& sol,
                   Array4<Real const> const& Ax, Array4<Real const> const& rhs,
                   Array4<Real const> const& sig, Array4<int const> const& msk,
                   GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void mlndhelm_gauss_seidel_ha(
    Box const& bx, Array4<Real> const& sol, Array4<Real const> const& rhs,
    Array4<Real const> const& sx, Array4<int const> const& msk,
    GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void mlndhelm_gauss_seidel_aa(
    Box const& bx, Array4<Real> const& sol, Array4<Real const> const& rhs,
    Array4<Real const> const& sig, Array4<int const> const& msk,
    GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_restriction(int i, int j, int k, Array4<Real> const& crse,
                     Array4<Real const> const& fine,
                     Array4<int const> const& msk) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_interpadd_aa(int i, int j, int k, Array4<Real> const& fine,
                      Array4<Real const> const& crse,
                      Array4<Real const> const& sig,
                      Array4<int const> const& msk) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_interpadd_ha(int i, int j, int k, Array4<Real> const& fine,
                      Array4<Real const> const& crse,
                      Array4<Real const> const& sigx,
                      Array4<int const> const& msk) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_divu(int i, int j, int k, Array4<Real> const& rhs,
              Array4<Real const> const& vel, Array4<int const> const& msk,
              GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
mlndhelm_rhcc(int i, int j, int k, Array4<Real const> const& rhcc,
              Array4<int const> const& msk) noexcept {
  return 0._rt;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_mknewu(int i, int j, int k, Array4<Real> const& u,
                Array4<Real const> const& p, Array4<Real const> const& sig,
                GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_divu_compute_fine_contrib(
    int i, int j, int k, Box const& fvbx, Array4<Real> const& frh,
    Array4<Real const> const& vel,
    GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void mlndhelm_divu_add_fine_contrib(
    int i, int j, int k, Box const& fvbx, Array4<Real> const& rhs,
    Array4<Real const> const& frh, Array4<int const> const& msk) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void mlndhelm_rhcc_fine_contrib(
    int i, int j, int k, Box const& fvbx, Array4<Real> const& rhs,
    Array4<Real const> const& cc, Array4<int const> const& msk) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void mlndhelm_divu_cf_contrib(
    int i, int j, int, Array4<Real> const& rhs, Array4<Real const> const& vel,
    Array4<Real const> const& fc, Array4<Real const> const& rhcc,
    Array4<int const> const& dmsk, Array4<int const> const& ndmsk,
    Array4<int const> const& ccmsk, GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
    Box const& nddom, GpuArray<LinOpBCType, AMREX_SPACEDIM> const& lobc,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& hibc) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void mlndhelm_crse_resid(
    int i, int j, int k, Array4<Real> const& resid,
    Array4<Real const> const& rhs, Array4<int const> const& msk,
    Box const& nddom, GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_res_fine_Ax(int i, int j, int k, Box const& fvbx,
                     Array4<Real> const& Ax, Array4<Real const> const& x,
                     Array4<Real const> const& sig,
                     GpuArray<Real, AMREX_SPACEDIM> const& dxinv) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlndhelm_res_fine_contrib(int i, int j, int k, Array4<Real> const& f,
                          Array4<Real const> const& Ax,
                          Array4<int const> const& msk) noexcept {}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void mlndhelm_res_cf_contrib(
    int i, int j, int k, Array4<Real> const& res, Array4<Real const> const& phi,
    Array4<Real const> const& rhs, Array4<Real const> const& sig,
    Array4<int const> const& dmsk, Array4<int const> const& ndmsk,
    Array4<int const> const& ccmsk, Array4<Real const> const& fc,
    GpuArray<Real, AMREX_SPACEDIM> const& dxinv, Box const& nddom,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& lobc,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& hibc) noexcept {}

} // namespace amrex

#endif
