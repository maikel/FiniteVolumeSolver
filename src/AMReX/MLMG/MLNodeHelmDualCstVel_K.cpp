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

#ifndef AMREX_ML_NODEHELM_DUAL_CSTVEL_K_H_
#define AMREX_ML_NODEHELM_DUAL_CSTVEL_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_LO_BCTYPES.H>

namespace amrex {
namespace {

[[maybe_unused]] constexpr int crse_cell = 0;
[[maybe_unused]] constexpr int fine_cell = 1;
[[maybe_unused]] constexpr int crse_node = 0;
[[maybe_unused]] constexpr int crse_fine_node = 1;
[[maybe_unused]] constexpr int fine_node = 2;
//#if (BL_USE_FLOAT)
//constexpr double eps = 1.e-30;
//#else
//constexpr double eps = 1.e-100;
//#endif
//constexpr Real almostone =
//    1._rt - 100._rt * std::numeric_limits<Real>::epsilon();

} // namespace

namespace {
inline void mlndhelm_scale_neumann_bc(
    Real s, Box const& bx, Array4<Real> const& rhs, Box const& nddom,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& lobc,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& hibc) noexcept {
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (lobc[idim] == LinOpBCType::Neumann or
        lobc[idim] == LinOpBCType::inflow) {
      Box const& blo = amrex::bdryLo(bx, idim);
      if (blo.smallEnd(idim) == nddom.smallEnd(idim)) {
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(blo, i, j, k, { rhs(i, j, k) *= s; });
      }
    }
    if (hibc[idim] == LinOpBCType::Neumann or
        hibc[idim] == LinOpBCType::inflow) {
      Box const& bhi = amrex::bdryHi(bx, idim);
      if (bhi.bigEnd(idim) == nddom.bigEnd(idim)) {
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bhi, i, j, k, { rhs(i, j, k) *= s; });
      }
    }
  }
}
} // namespace

inline void mlndhelm_impose_neumann_bc(
    Box const& bx, Array4<Real> const& rhs, Box const& nddom,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& lobc,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& hibc) noexcept {
  mlndhelm_scale_neumann_bc(2.0, bx, rhs, nddom, lobc, hibc);
}

inline void mlndhelm_unimpose_neumann_bc(
    Box const& bx, Array4<Real> const& rhs, Box const& nddom,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& lobc,
    GpuArray<LinOpBCType, AMREX_SPACEDIM> const& hibc) noexcept {
  mlndhelm_scale_neumann_bc(0.5, bx, rhs, nddom, lobc, hibc);
}

} // namespace amrex

#if (AMREX_SPACEDIM == 1)
#include <src/AMReX/MLMG/MLNodeHelmDualCstVel_1D_K.cpp>
#elif (AMREX_SPACEDIM == 2)
#include <src/AMReX/MLMG/MLNodeHelmDualCstVel_2D_K.cpp>
#else
#include <src/AMReX/MLMG/MLNodeHelmDualCstVel_3D_K.cpp>
#endif

namespace amrex {

template <typename T>
void mlndhelm_fillbc_cc(Box const& vbx, Array4<T> const& sigma,
                        Box const& domain,
                        GpuArray<LinOpBCType, AMREX_SPACEDIM> bclo,
                        GpuArray<LinOpBCType, AMREX_SPACEDIM> bchi) noexcept {
  GpuArray<bool, AMREX_SPACEDIM> bflo{AMREX_D_DECL(
      bclo[0] != LinOpBCType::Periodic, bclo[1] != LinOpBCType::Periodic,
      bclo[2] != LinOpBCType::Periodic)};
  GpuArray<bool, AMREX_SPACEDIM> bfhi{AMREX_D_DECL(
      bchi[0] != LinOpBCType::Periodic, bchi[1] != LinOpBCType::Periodic,
      bchi[2] != LinOpBCType::Periodic)};
  mlndhelm_bc_doit(vbx, sigma, domain, bflo, bfhi);
}

template <typename T>
void mlndhelm_applybc(Box const& vbx, Array4<T> const& phi, Box const& domain,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> bclo,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> bchi) noexcept {
  GpuArray<bool, AMREX_SPACEDIM> bflo{AMREX_D_DECL(
      bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow,
      bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow,
      bclo[2] == LinOpBCType::Neumann or bclo[2] == LinOpBCType::inflow)};
  GpuArray<bool, AMREX_SPACEDIM> bfhi{AMREX_D_DECL(
      bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow,
      bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow,
      bchi[2] == LinOpBCType::Neumann or bchi[2] == LinOpBCType::inflow)};
  mlndhelm_bc_doit(vbx, phi, domain, bflo, bfhi);
}

AMREX_FORCE_INLINE
bool mlndhelm_any_fine_sync_cells(Box const& bx, Array4<int const> const& msk,
                                  int fine_flag) noexcept {
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);
  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        if (msk(i, j, k) == fine_flag)
          return true;
      }
    }
  }
  return false;
}

} // namespace amrex

#endif
