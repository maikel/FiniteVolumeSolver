#ifndef AMREX_ML_NODEHELM_DUAL_CSTVEL_3D_K_H_
#define AMREX_ML_NODEHELM_DUAL_CSTVEL_3D_K_H_

namespace amrex {

//
// masks
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_set_nodal_mask (int i, int j, int k, Array4<int> const& nmsk,
                             Array4<int const> const& cmsk) noexcept
{
    int s = cmsk(i-1,j-1,k-1) + cmsk(i  ,j-1,k-1)
        +   cmsk(i-1,j  ,k-1) + cmsk(i  ,j  ,k-1)
        +   cmsk(i-1,j-1,k  ) + cmsk(i  ,j-1,k  )
        +   cmsk(i-1,j  ,k  ) + cmsk(i  ,j  ,k  );
    if (s == 8*crse_cell) {
        nmsk(i,j,k) = crse_node;
    }
    else if (s == 8*fine_cell) {
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
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        dmsk(i,j,k) = (omsk(i-1,j-1,k-1) == 1 or omsk(i,j-1,k-1) == 1 or
                       omsk(i-1,j  ,k-1) == 1 or omsk(i,j  ,k-1) == 1 or
                       omsk(i-1,j-1,k  ) == 1 or omsk(i,j-1,k  ) == 1 or
                       omsk(i-1,j  ,k  ) == 1 or omsk(i,j  ,k  ) == 1);
    }}}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if (bclo[0] == LinOpBCType::Dirichlet and lo.x == domlo.x) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,k) = 1;
        }}
    }

    if (bchi[0] == LinOpBCType::Dirichlet and hi.x == domhi.x) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,k) = 1;
        }}
    }

    if (bclo[1] == LinOpBCType::Dirichlet and lo.y == domlo.y) {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,k) = 1;
        }}
    }

    if (bchi[1] == LinOpBCType::Dirichlet and hi.y == domhi.y) {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,k) = 1;
        }}
    }

    if (bclo[2] == LinOpBCType::Dirichlet and lo.z == domlo.z) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,lo.z) = 1;
        }}
    }

    if (bchi[2] == LinOpBCType::Dirichlet and hi.z == domhi.z) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,hi.z) = 1;
        }}
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
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        dmsk(i,j,k) = static_cast<Real>(omsk(i,j,k));
    }}}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if ((bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow)
        and lo.x == domlo.x)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,k) *= 0.5;
        }}
    }

    if ((bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow)
        and hi.x == domhi.x)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,k) *= 0.5;
        }}
    }

    if ((bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow)
        and lo.y == domlo.y)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,k) *= 0.5;
        }}
    }

    if ((bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow)
        and hi.y == domhi.y)
    {
        for (int k = lo.z; k <= hi.z; ++k) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,k) *= 0.5;
        }}
    }

    if ((bclo[2] == LinOpBCType::Neumann or bclo[2] == LinOpBCType::inflow)
        and lo.z == domlo.z)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,lo.z) *= 0.5;
        }}
    }

    if ((bchi[2] == LinOpBCType::Neumann or bchi[2] == LinOpBCType::inflow)
        and hi.z == domhi.z)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,j,hi.z) *= 0.5;
        }}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_zero_fine (int i, int j, int k, Array4<Real> const& phi,
                        Array4<int const> const& msk, int fine_flag) noexcept
{
    // Testing if the node is covered by a fine level in computing
    // coarse sync residual
    if (msk(i-1,j-1,k-1) == fine_flag and
        msk(i  ,j-1,k-1) == fine_flag and
        msk(i-1,j  ,k-1) == fine_flag and
        msk(i  ,j  ,k-1) == fine_flag and
        msk(i-1,j-1,k  ) == fine_flag and
        msk(i  ,j-1,k  ) == fine_flag and
        msk(i-1,j  ,k  ) == fine_flag and
        msk(i  ,j  ,k  ) == fine_flag)
    {
        phi(i,j,k) = 0.0;
    }
}

//
// coeffs
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_avgdown_coeff_x (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    Real cl = fine(2*i  ,2*j,2*k  )+fine(2*i  ,2*j+1,2*k  )+
              fine(2*i  ,2*j,2*k+1)+fine(2*i  ,2*j+1,2*k+1);
    Real cr = fine(2*i+1,2*j,2*k  )+fine(2*i+1,2*j+1,2*k  )+
              fine(2*i+1,2*j,2*k+1)+fine(2*i+1,2*j+1,2*k+1);
    crse(i,j,k) = 0.5*cl*cr/(cl+cr);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_avgdown_coeff_y (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    Real cl = fine(2*i,2*j  ,2*k  )+fine(2*i+1,2*j  ,2*k  )+
              fine(2*i,2*j  ,2*k+1)+fine(2*i+1,2*j  ,2*k+1);
    Real cr = fine(2*i,2*j+1,2*k  )+fine(2*i+1,2*j+1,2*k  )+
              fine(2*i,2*j+1,2*k+1)+fine(2*i+1,2*j+1,2*k+1);
    crse(i,j,k) = 0.5*cl*cr/(cl+cr);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_avgdown_coeff_z (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    Real cl = fine(2*i,2*j  ,2*k  )+fine(2*i+1,2*j  ,2*k  )+
              fine(2*i,2*j+1,2*k  )+fine(2*i+1,2*j+1,2*k  );
    Real cr = fine(2*i,2*j  ,2*k+1)+fine(2*i+1,2*j  ,2*k+1)+
              fine(2*i,2*j+1,2*k+1)+fine(2*i+1,2*j+1,2*k+1);
    crse(i,j,k) = 0.5*cl*cr/(cl+cr);
}

//
// bc
//

template <typename T>
inline void mlndhelm_bc_doit (Box const& vbx, Array4<T> const& a, Box const& domain,
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
    Box zlo_face = domain;
    zlo_face.setSmall(2,dlo.z-1);
    zlo_face.setBig  (2,dlo.z-1);
    zlo_face &= sbox;
    int xoffset = vbx.length(0)+1;
    int yoffset = vbx.length(1)+1;
    int zoffset = vbx.length(2)+1;

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
    xlo_face, txbxlo,
    {
        auto lo = amrex::lbound(txbxlo);
        auto hi = amrex::ubound(txbxlo);
        if (lo.x == dlo.x-1 and bflo[0]) {
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                a(dlo.x-1,j,k) = a(dlo.x+offset,j,k);
            }}
        }
        if (lo.x+xoffset == dhi.x+1 and bfhi[0]) {
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                a(dhi.x+1,j,k) = a(dhi.x-offset,j,k);
            }}
        }
    },
    ylo_face, tybxlo,
    {
        auto lo = amrex::lbound(tybxlo);
        auto hi = amrex::ubound(tybxlo);
        if (lo.y == dlo.y-1 and bflo[1]) {
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,dlo.y-1,k) = a(i,dlo.y+offset,k);
            }}
        }
        if (lo.y+yoffset == dhi.y+1 and bfhi[1]) {
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,dhi.y+1,k) = a(i,dhi.y-offset,k);
            }}
        }
    },
    zlo_face, tzbxlo,
    {
        auto lo = amrex::lbound(tzbxlo);
        auto hi = amrex::ubound(tzbxlo);
        if (lo.z == dlo.z-1 and bflo[2]) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,j,dlo.z-1) = a(i,j,dlo.z+offset);
            }}
        }
        if (lo.z+zoffset == dhi.z+1 and bfhi[2]) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,j,dhi.z+1) = a(i,j,dhi.z-offset);
            }}
        }
    });

    const auto lo = amrex::lbound(sbox);
    const auto hi = amrex::ubound(sbox);

    AMREX_HOST_DEVICE_FOR_1D ( 12, iedge,
    {
        switch (iedge) {
        case 0: {
            // xlo & ylo
            if (lo.x == dlo.x-1 and lo.y == dlo.y-1) {
                if (bflo[0]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dlo.x-1,dlo.y-1,k) = a(dlo.x+offset,dlo.y-1,k);
                    }
                } else if (bflo[1]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dlo.x-1,dlo.y-1,k) = a(dlo.x-1,dlo.y+offset,k);
                    }
                }
            }
            break;
        }
        case 1: {
            // xhi & ylo
            if (hi.x == dhi.x+1 and lo.y == dlo.y-1) {
                if (bfhi[0]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dhi.x+1,dlo.y-1,k) = a(dhi.x-offset,dlo.y-1,k);
                    }
                } else if (bflo[1]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dhi.x+1,dlo.y-1,k) = a(dhi.x+1,dlo.y+offset,k);
                    }
                }
            }
            break;
        }
        case 2: {
            // xlo & yhi
            if (lo.x == dlo.x-1 and hi.y == dhi.y+1) {
                if (bflo[0]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dlo.x-1,dhi.y+1,k) = a(dlo.x+offset,dhi.y+1,k);
                    }
                } else if (bfhi[1]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dlo.x-1,dhi.y+1,k) = a(dlo.x-1,dhi.y-offset,k);
                    }
                }
            }
            break;
        }
        case 3: {
            // xhi & yhi
            if (hi.x == dhi.x+1 and hi.y == dhi.y+1) {
                if (bfhi[0]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dhi.x+1,dhi.y+1,k) = a(dhi.x-offset,dhi.y+1,k);
                    }
                } else if (bfhi[1]) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        a(dhi.x+1,dhi.y+1,k) = a(dhi.x+1,dhi.y-offset,k);
                    }
                }
            }
            break;
        }
        case 4: {
            // xlo & zlo
            if (lo.x == dlo.x-1 and lo.z == dlo.z-1) {
                if (bflo[0]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dlo.x-1,j,dlo.z-1) = a(dlo.x+offset,j,dlo.z-1);
                    }
                } else if (bflo[2]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dlo.x-1,j,dlo.z-1) = a(dlo.x-1,j,dlo.z+offset);
                    }
                }
            }
            break;
        }
        case 5: {
            // xhi & zlo
            if (hi.x == dhi.x+1 and lo.z == dlo.z-1) {
                if (bfhi[0]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dhi.x+1,j,dlo.z-1) = a(dhi.x-offset,j,dlo.z-1);
                    }
                } else if (bflo[2]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dhi.x+1,j,dlo.z-1) = a(dhi.x+1,j,dlo.z+offset);
                    }
                }
            }
            break;
        }
        case 6: {
            // xlo & zhi
            if (lo.x == dlo.x-1 and hi.z == dhi.z+1) {
                if (bflo[0]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dlo.x-1,j,dhi.z+1) = a(dlo.x+offset,j,dhi.z+1);
                    }
                } else if (bfhi[2]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dlo.x-1,j,dhi.z+1) = a(dlo.x-1,j,dhi.z-offset);
                    }
                }
            }
            break;
        }
        case 7: {
            // xhi & zhi
            if (hi.x == dhi.x+1 and hi.z == dhi.z+1) {
                if (bfhi[0]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dhi.x+1,j,dhi.z+1) = a(dhi.x-offset,j,dhi.z+1);
                    }
                } else if (bfhi[2]) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        a(dhi.x+1,j,dhi.z+1) = a(dhi.x+1,j,dhi.z-offset);
                    }
                }
            }
            break;
        }
        case 8: {
            // ylo & zlo
            if (lo.y == dlo.y-1 and lo.z == dlo.z-1) {
                if (bflo[1]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dlo.y-1,dlo.z-1) = a(i,dlo.y+offset,dlo.z-1);
                    }
                } else if (bflo[2]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dlo.y-1,dlo.z-1) = a(i,dlo.y-1,dlo.z+offset);
                    }
                }
            }
            break;
        }
        case 9: {
            // yhi & zlo
            if (hi.y == dhi.y+1 and lo.z == dlo.z-1) {
                if (bfhi[1]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dhi.y+1,dlo.z-1) = a(i,dhi.y-offset,dlo.z-1);
                    }
                } else if (bflo[2]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dhi.y+1,dlo.z-1) = a(i,dhi.y+1,dlo.z+offset);
                    }
                }
            }
            break;
        }
        case 10: {
            // ylo & zhi
            if (lo.y == dlo.y-1 and hi.z == dhi.z+1) {
                if (bflo[1]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dlo.y-1,dhi.z+1) = a(i,dlo.y+offset,dhi.z+1);
                    }
                } else if (bfhi[2]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dlo.y-1,dhi.z+1) = a(i,dlo.y-1,dhi.z-offset);
                    }
                }
            }
            break;
        }
        case 11: {
            // yhi & zhi
            if (hi.y == dhi.y+1 and hi.z == dhi.z+1) {
                if (bfhi[1]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dhi.y+1,dhi.z+1) = a(i,dhi.y-offset,dhi.z+1);
                    }
                } else if (bfhi[2]) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        a(i,dhi.y+1,dhi.z+1) = a(i,dhi.y+1,dhi.z-offset);
                    }
                }
            }
            break;
        }
        default: {}
        }
    });

    AMREX_HOST_DEVICE_FOR_1D ( 8, icorner,
    {
        switch (icorner) {
        case 0: {
            // xlo & ylo & zlo
            if (lo.x == dlo.x-1 and lo.y == dlo.y-1 and lo.z == dlo.z-1) {
                if (bflo[0]) {
                    a(dlo.x-1,dlo.y-1,dlo.z-1) = a(dlo.x+offset,dlo.y-1,dlo.z-1);
                } else if (bflo[1]) {
                    a(dlo.x-1,dlo.y-1,dlo.z-1) = a(dlo.x-1,dlo.y+offset,dlo.z-1);
                } else if (bflo[2]) {
                    a(dlo.x-1,dlo.y-1,dlo.z-1) = a(dlo.x-1,dlo.y-1,dlo.z+offset);
                }
            }
            break;
        }
        case 1: {
            // xhi & ylo & zlo
            if (hi.x == dhi.x+1 and lo.y == dlo.y-1 and lo.z == dlo.z-1) {
                if (bfhi[0]) {
                    a(dhi.x+1,dlo.y-1,dlo.z-1) = a(dhi.x-offset,dlo.y-1,dlo.z-1);
                } else if (bflo[1]) {
                    a(dhi.x+1,dlo.y-1,dlo.z-1) = a(dhi.x+1,dlo.y+offset,dlo.z-1);
                } else if (bflo[2]) {
                    a(dhi.x+1,dlo.y-1,dlo.z-1) = a(dhi.x+1,dlo.y-1,dlo.z+offset);
                }
            }
            break;
        }
        case 2: {
            // xlo & yhi & zlo
            if (lo.x == dlo.x-1 and hi.y == dhi.y+1 and lo.z == dlo.z-1) {
                if (bflo[0]) {
                    a(dlo.x-1,dhi.y+1,dlo.z-1) = a(dlo.x+offset,dhi.y+1,dlo.z-1);
                } else if (bfhi[1]) {
                    a(dlo.x-1,dhi.y+1,dlo.z-1) = a(dlo.x-1,dhi.y-offset,dlo.z-1);
                } else if (bflo[2]) {
                    a(dlo.x-1,dhi.y+1,dlo.z-1) = a(dlo.x-1,dhi.y+1,dlo.z+offset);
                }
            }
            break;
        }
        case 3: {
            // xhi & yhi & zlo
            if (hi.x == dhi.x+1 and hi.y == dhi.y+1 and lo.z == dlo.z-1) {
                if (bfhi[0]) {
                    a(dhi.x+1,dhi.y+1,dlo.z-1) = a(dhi.x-offset,dhi.y+1,dlo.z-1);
                } else if (bfhi[1]) {
                    a(dhi.x+1,dhi.y+1,dlo.z-1) = a(dhi.x+1,dhi.y-offset,dlo.z-1);
                } else if (bflo[2]) {
                    a(dhi.x+1,dhi.y+1,dlo.z-1) = a(dhi.x+1,dhi.y+1,dlo.z+offset);
                }
            }
            break;
        }
        case 4: {
            // xlo & ylo & zhi
            if (lo.x == dlo.x-1 and lo.y == dlo.y-1 and hi.z == dhi.z+1) {
                if (bflo[0]) {
                    a(dlo.x-1,dlo.y-1,dhi.z+1) = a(dlo.x+offset,dlo.y-1,dhi.z+1);
                } else if (bflo[1]) {
                    a(dlo.x-1,dlo.y-1,dhi.z+1) = a(dlo.x-1,dlo.y+offset,dhi.z+1);
                } else if (bfhi[2]) {
                    a(dlo.x-1,dlo.y-1,dhi.z+1) = a(dlo.x-1,dlo.y-1,dhi.z-offset);
                }
            }
            break;
        }
        case 5: {
            // xhi & ylo & zhi
            if (hi.x == dhi.x+1 and lo.y == dlo.y-1 and hi.z == dhi.z+1) {
                if (bfhi[0]) {
                    a(dhi.x+1,dlo.y-1,dhi.z+1) = a(dhi.x-offset,dlo.y-1,dhi.z+1);
                } else if (bflo[1]) {
                    a(dhi.x+1,dlo.y-1,dhi.z+1) = a(dhi.x+1,dlo.y+offset,dhi.z+1);
                } else if (bfhi[2]) {
                    a(dhi.x+1,dlo.y-1,dhi.z+1) = a(dhi.x+1,dlo.y-1,dhi.z-offset);
                }
            }
            break;
        }
        case 6: {
            // xlo & yhi & zhi
            if (lo.x == dlo.x-1 and hi.y == dhi.y+1 and hi.z == dhi.z+1) {
                if (bflo[0]) {
                    a(dlo.x-1,dhi.y+1,dhi.z+1) = a(dlo.x+offset,dhi.y+1,dhi.z+1);
                } else if (bfhi[1]) {
                    a(dlo.x-1,dhi.y+1,dhi.z+1) = a(dlo.x-1,dhi.y-offset,dhi.z+1);
                } else if (bfhi[2]) {
                    a(dlo.x-1,dhi.y+1,dhi.z+1) = a(dlo.x-1,dhi.y+1,dhi.z-offset);
                }
            }
            break;
        }
        case 7: {
            // xhi & yhi & zhi
            if (hi.x == dhi.x+1 and hi.y == dhi.y+1 and hi.z == dhi.z+1) {
                if (bfhi[0] ) {
                    a(dhi.x+1,dhi.y+1,dhi.z+1) = a(dhi.x-offset,dhi.y+1,dhi.z+1);
                } else if (bfhi[1]) {
                    a(dhi.x+1,dhi.y+1,dhi.z+1) = a(dhi.x+1,dhi.y-offset,dhi.z+1);
                } else if (bfhi[2]) {
                    a(dhi.x+1,dhi.y+1,dhi.z+1) = a(dhi.x+1,dhi.y+1,dhi.z-offset);
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
                       Array4<Real const> const& sz, Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        y(i,j,k) = 0.0;
    } else {
        Real facx = (1./36.)*dxinv[0]*dxinv[0];
        Real facy = (1./36.)*dxinv[1]*dxinv[1];
        Real facz = (1./36.)*dxinv[2]*dxinv[2];
        y(i,j,k) = x(i,j,k)*(-4.0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                                         +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                   +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                                         +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                   +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                                         +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
        y(i,j,k) += x(i-1,j-1,k-1)*(facx*sx(i-1,j-1,k-1)
                                   +facy*sy(i-1,j-1,k-1)
                                   +facz*sz(i-1,j-1,k-1))
                  + x(i+1,j-1,k-1)*(facx*sx(i  ,j-1,k-1)
                                   +facy*sy(i  ,j-1,k-1)
                                   +facz*sz(i  ,j-1,k-1))
                  + x(i-1,j+1,k-1)*(facx*sx(i-1,j  ,k-1)
                                   +facy*sy(i-1,j  ,k-1)
                                   +facz*sz(i-1,j  ,k-1))
                  + x(i+1,j+1,k-1)*(facx*sx(i  ,j  ,k-1)
                                   +facy*sy(i  ,j  ,k-1)
                                   +facz*sz(i  ,j  ,k-1))
                  + x(i-1,j-1,k+1)*(facx*sx(i-1,j-1,k  )
                                   +facy*sy(i-1,j-1,k  )
                                   +facz*sz(i-1,j-1,k  ))
                  + x(i+1,j-1,k+1)*(facx*sx(i  ,j-1,k  )
                                   +facy*sy(i  ,j-1,k  )
                                   +facz*sz(i  ,j-1,k  ))
                  + x(i-1,j+1,k+1)*(facx*sx(i-1,j  ,k  )
                                   +facy*sy(i-1,j  ,k  )
                                   +facz*sz(i-1,j  ,k  ))
                  + x(i+1,j+1,k+1)*(facx*sx(i  ,j  ,k  )
                                   +facy*sy(i  ,j  ,k  )
                                   +facz*sz(i  ,j  ,k  ));
        y(i,j,k) += x(i  ,j-1,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1))
                                    +2.0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1))
                                    +2.0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)))
                  + x(i  ,j+1,k-1)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1))
                                    +2.0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1))
                                    +2.0*facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)))
                  + x(i  ,j-1,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  ))
                                    +2.0*facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  ))
                                    +2.0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )))
                  + x(i  ,j+1,k+1)*(    -facx*(sx(i-1,j  ,k  )+sx(i,j  ,k  ))
                                    +2.0*facy*(sy(i-1,j  ,k  )+sy(i,j  ,k  ))
                                    +2.0*facz*(sz(i-1,j  ,k  )+sz(i,j  ,k  )))
                  + x(i-1,j  ,k-1)*( 2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1))
                                        -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1))
                                    +2.0*facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)))
                  + x(i+1,j  ,k-1)*( 2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1))
                                        -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1))
                                    +2.0*facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)))
                  + x(i-1,j  ,k+1)*( 2.0*facx*(sx(i-1,j-1,k  )+sx(i-1,j,k  ))
                                        -facy*(sy(i-1,j-1,k  )+sy(i-1,j,k  ))
                                    +2.0*facz*(sz(i-1,j-1,k  )+sz(i-1,j,k  )))
                  + x(i+1,j  ,k+1)*( 2.0*facx*(sx(i  ,j-1,k  )+sx(i  ,j,k  ))
                                        -facy*(sy(i  ,j-1,k  )+sy(i  ,j,k  ))
                                    +2.0*facz*(sz(i  ,j-1,k  )+sz(i  ,j,k  )))
                  + x(i-1,j-1,k  )*( 2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j-1,k))
                                    +2.0*facy*(sy(i-1,j-1,k-1)+sy(i-1,j-1,k))
                                        -facz*(sz(i-1,j-1,k-1)+sz(i-1,j-1,k)))
                  + x(i+1,j-1,k  )*( 2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j-1,k))
                                    +2.0*facy*(sy(i  ,j-1,k-1)+sy(i  ,j-1,k))
                                        -facz*(sz(i  ,j-1,k-1)+sz(i  ,j-1,k)))
                  + x(i-1,j+1,k  )*( 2.0*facx*(sx(i-1,j  ,k-1)+sx(i-1,j  ,k))
                                    +2.0*facy*(sy(i-1,j  ,k-1)+sy(i-1,j  ,k))
                                        -facz*(sz(i-1,j  ,k-1)+sz(i-1,j  ,k)))
                  + x(i+1,j+1,k  )*( 2.0*facx*(sx(i  ,j  ,k-1)+sx(i  ,j  ,k))
                                    +2.0*facy*(sy(i  ,j  ,k-1)+sy(i  ,j  ,k))
                                        -facz*(sz(i  ,j  ,k-1)+sz(i  ,j  ,k)));
            y(i,j,k) += 2.0*x(i-1,j,k)*( 2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)+sx(i-1,j-1,k)+sx(i-1,j,k))
                                            -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)+sy(i-1,j-1,k)+sy(i-1,j,k))
                                            -facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)+sz(i-1,j-1,k)+sz(i-1,j,k)))
                      + 2.0*x(i+1,j,k)*( 2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k))
                                            -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k))
                                            -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k)))
                      + 2.0*x(i,j-1,k)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k))
                                        +2.0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k))
                                            -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k)))
                      + 2.0*x(i,j+1,k)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k))
                                        +2.0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k))
                                            -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k)))
                      + 2.0*x(i,j,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1))
                                            -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1))
                                        +2.0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)))
                      + 2.0*x(i,j,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                            -facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                        +2.0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
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
        Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
        Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
        Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];
        Real fxyz = facx + facy + facz;
        Real fmx2y2z = -facx + 2.0*facy + 2.0*facz;
        Real f2xmy2z = 2.0*facx - facy + 2.0*facz;
        Real f2x2ymz = 2.0*facx + 2.0*facy - facz;
        Real f4xm2ym2z = 4.0*facx - 2.0*facy - 2.0*facz;
        Real fm2x4ym2z = -2.0*facx + 4.0*facy - 2.0*facz;
        Real fm2xm2y4z = -2.0*facx - 2.0*facy + 4.0*facz;
        y(i,j,k) = x(i,j,k)*(-4.0)*fxyz*
            (sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
            +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ))
            + fxyz*(x(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                  + x(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                  + x(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                  + x(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                  + x(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                  + x(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                  + x(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                  + x(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
            + fmx2y2z*(x(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                     + x(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                     + x(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                     + x(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
            + f2xmy2z*(x(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                     + x(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                     + x(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                     + x(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
            + f2x2ymz*(x(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                     + x(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                     + x(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                     + x(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
            + f4xm2ym2z*(x(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                       + x(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
            + fm2x4ym2z*(x(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                       + x(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
            + fm2xm2y4z*(x(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                       + x(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_normalize_ha (Box const& bx, Array4<Real> const& x, Array4<Real const> const& sx,
                           Array4<Real const> const& sy, Array4<Real const> const& sz,
                           Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k)) {
            x(i,j,k) = x(i,j,k)/((-4.0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                                              +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                        +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                                              +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                        +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                                              +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  ))));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_normalize_aa (Box const& bx, Array4<Real> const& x, Array4<Real const> const& sig,
                           Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];
    Real fxyz = facx + facy + facz;

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k)) {
            x(i,j,k) = x(i,j,k) /
                ((-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                             +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_jacobi_ha (Box const& bx, Array4<Real> const& sol, Array4<Real const> const& Ax,
                        Array4<Real const> const& rhs, Array4<Real const> const& sx,
                        Array4<Real const> const& sy, Array4<Real const> const& sz,
                        Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = -4.0 * (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = -4.0 * (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = -4.0 * (1.0/36.0)*dxinv[2]*dxinv[2];

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            sol(i,j,k) += (2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                        +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                  +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                        +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                  +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                        +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_jacobi_aa (Box const& bx, Array4<Real> const& sol, Array4<Real const> const& Ax,
                        Array4<Real const> const& rhs, Array4<Real const> const& sig,
                        Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real fxyz = -4.0 * (1.0/36.0)*(dxinv[0]*dxinv[0] +
                                   dxinv[1]*dxinv[1] +
                                   dxinv[2]*dxinv[2]);

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            sol(i,j,k) += (2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                        +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_gauss_seidel_ha (Box const& bx, Array4<Real> const& sol,
                              Array4<Real const> const& rhs, Array4<Real const> const& sx,
                              Array4<Real const> const& sy, Array4<Real const> const& sz,
                              Array4<int const> const& msk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];

    amrex::Loop(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            Real s0 = (-4.0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)
                                   +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                             +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)
                                   +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                             +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)
                                   +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));
            Real Ax = sol(i,j,k)*s0
                     + sol(i-1,j-1,k-1)*(facx*sx(i-1,j-1,k-1)
                                        +facy*sy(i-1,j-1,k-1)
                                        +facz*sz(i-1,j-1,k-1))
                     + sol(i+1,j-1,k-1)*(facx*sx(i  ,j-1,k-1)
                                        +facy*sy(i  ,j-1,k-1)
                                        +facz*sz(i  ,j-1,k-1))
                     + sol(i-1,j+1,k-1)*(facx*sx(i-1,j  ,k-1)
                                        +facy*sy(i-1,j  ,k-1)
                                        +facz*sz(i-1,j  ,k-1))
                     + sol(i+1,j+1,k-1)*(facx*sx(i  ,j  ,k-1)
                                        +facy*sy(i  ,j  ,k-1)
                                        +facz*sz(i  ,j  ,k-1))
                     + sol(i-1,j-1,k+1)*(facx*sx(i-1,j-1,k  )
                                        +facy*sy(i-1,j-1,k  )
                                        +facz*sz(i-1,j-1,k  ))
                     + sol(i+1,j-1,k+1)*(facx*sx(i  ,j-1,k  )
                                        +facy*sy(i  ,j-1,k  )
                                        +facz*sz(i  ,j-1,k  ))
                     + sol(i-1,j+1,k+1)*(facx*sx(i-1,j  ,k  )
                                        +facy*sy(i-1,j  ,k  )
                                        +facz*sz(i-1,j  ,k  ))
                     + sol(i+1,j+1,k+1)*(facx*sx(i  ,j  ,k  )
                                        +facy*sy(i  ,j  ,k  )
                                        +facz*sz(i  ,j  ,k  ))
                     +sol(i  ,j-1,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1))
                                        +2.0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1))
                                        +2.0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)))
                     +sol(i  ,j+1,k-1)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1))
                                        +2.0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1))
                                        +2.0*facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)))
                     +sol(i  ,j-1,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  ))
                                        +2.0*facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  ))
                                        +2.0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )))
                     +sol(i  ,j+1,k+1)*(    -facx*(sx(i-1,j  ,k  )+sx(i,j  ,k  ))
                                        +2.0*facy*(sy(i-1,j  ,k  )+sy(i,j  ,k  ))
                                        +2.0*facz*(sz(i-1,j  ,k  )+sz(i,j  ,k  )))
                     +sol(i-1,j  ,k-1)*( 2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1))
                                            -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1))
                                        +2.0*facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)))
                     +sol(i+1,j  ,k-1)*( 2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1))
                                            -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1))
                                        +2.0*facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)))
                     +sol(i-1,j  ,k+1)*( 2.0*facx*(sx(i-1,j-1,k  )+sx(i-1,j,k  ))
                                            -facy*(sy(i-1,j-1,k  )+sy(i-1,j,k  ))
                                        +2.0*facz*(sz(i-1,j-1,k  )+sz(i-1,j,k  )))
                     +sol(i+1,j  ,k+1)*( 2.0*facx*(sx(i  ,j-1,k  )+sx(i  ,j,k  ))
                                            -facy*(sy(i  ,j-1,k  )+sy(i  ,j,k  ))
                                        +2.0*facz*(sz(i  ,j-1,k  )+sz(i  ,j,k  )))
                     +sol(i-1,j-1,k  )*( 2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j-1,k))
                                        +2.0*facy*(sy(i-1,j-1,k-1)+sy(i-1,j-1,k))
                                            -facz*(sz(i-1,j-1,k-1)+sz(i-1,j-1,k)))
                     +sol(i+1,j-1,k  )*( 2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j-1,k))
                                        +2.0*facy*(sy(i  ,j-1,k-1)+sy(i  ,j-1,k))
                                            -facz*(sz(i  ,j-1,k-1)+sz(i  ,j-1,k)))
                     +sol(i-1,j+1,k  )*( 2.0*facx*(sx(i-1,j  ,k-1)+sx(i-1,j  ,k))
                                        +2.0*facy*(sy(i-1,j  ,k-1)+sy(i-1,j  ,k))
                                            -facz*(sz(i-1,j  ,k-1)+sz(i-1,j  ,k)))
                     +sol(i+1,j+1,k  )*( 2.0*facx*(sx(i  ,j  ,k-1)+sx(i  ,j  ,k))
                                        +2.0*facy*(sy(i  ,j  ,k-1)+sy(i  ,j  ,k))
                                            -facz*(sz(i  ,j  ,k-1)+sz(i  ,j  ,k)))
                     + 2.0*sol(i-1,j,k)*(2.0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)+sx(i-1,j-1,k)+sx(i-1,j,k))
                                            -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)+sy(i-1,j-1,k)+sy(i-1,j,k))
                                            -facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)+sz(i-1,j-1,k)+sz(i-1,j,k)))
                     + 2.0*sol(i+1,j,k)*(2.0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k))
                                            -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k))
                                            -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k)))
                     + 2.0*sol(i,j-1,k)*(   -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k))
                                        +2.0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k))
                                            -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k)))
                     + 2.0*sol(i,j+1,k)*(   -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k))
                                        +2.0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k))
                                            -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k)))
                     + 2.0*sol(i,j,k-1)*(   -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1))
                                            -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1))
                                        +2.0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1)))
                     + 2.0*sol(i,j,k+1)*(   -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  ))
                                            -facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  ))
                                        +2.0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )));

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
    Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];
    Real fxyz = facx + facy + facz;
    Real fmx2y2z = -facx + 2.0*facy + 2.0*facz;
    Real f2xmy2z = 2.0*facx - facy + 2.0*facz;
    Real f2x2ymz = 2.0*facx + 2.0*facy - facz;
    Real f4xm2ym2z = 4.0*facx - 2.0*facy - 2.0*facz;
    Real fm2x4ym2z = -2.0*facx + 4.0*facy - 2.0*facz;
    Real fm2xm2y4z = -2.0*facx - 2.0*facy + 4.0*facz;

    amrex::Loop(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            Real s0 = (-4.0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
                                  +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ));
            Real Ax = sol(i,j,k)*s0
                + fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                      + sol(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                      + sol(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                      + sol(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                      + sol(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                      + sol(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                      + sol(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                      + sol(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
                + fmx2y2z*(sol(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                         + sol(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                         + sol(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                         + sol(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
                + f2xmy2z*(sol(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                         + sol(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                         + sol(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                         + sol(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
                + f2x2ymz*(sol(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                         + sol(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                         + sol(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                         + sol(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
                + f4xm2ym2z*(sol(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                           + sol(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
                + fm2x4ym2z*(sol(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                           + sol(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
                + fm2xm2y4z*(sol(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                           + sol(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));

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
    int kk = k*2;
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = 0.0;
    } else {
        crse(i,j,k) = (1./64.)*(fine(ii-1,jj-1,kk-1)+fine(ii+1,jj-1,kk-1)
                               +fine(ii-1,jj+1,kk-1)+fine(ii+1,jj+1,kk-1)
                               +fine(ii-1,jj-1,kk+1)+fine(ii+1,jj-1,kk+1)
                               +fine(ii-1,jj+1,kk+1)+fine(ii+1,jj+1,kk+1))
                    + (1./32.)*(fine(ii  ,jj-1,kk-1)+fine(ii  ,jj+1,kk-1)
                               +fine(ii  ,jj-1,kk+1)+fine(ii  ,jj+1,kk+1)
                               +fine(ii-1,jj  ,kk-1)+fine(ii+1,jj  ,kk-1)
                               +fine(ii-1,jj  ,kk+1)+fine(ii+1,jj  ,kk+1)
                               +fine(ii-1,jj-1,kk  )+fine(ii+1,jj-1,kk  )
                               +fine(ii-1,jj+1,kk  )+fine(ii+1,jj+1,kk  ))
                    + (1./16.)*(fine(ii-1,jj,kk)+fine(ii+1,jj,kk)
                               +fine(ii,jj-1,kk)+fine(ii,jj+1,kk)
                               +fine(ii,jj,kk-1)+fine(ii,jj,kk+1))
                      + (1./8.)*fine(ii,jj,kk);
    }
}

//
// interpolation
//

namespace {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_x (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        return (w1*crse(ic,jc,kc)+w2*crse(ic+1,jc,kc))/(w1+w2);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_y (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        Real w2 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        return (w1*crse(ic,jc,kc)+w2*crse(ic,jc+1,kc))/(w1+w2);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_z (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        Real w2 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1*crse(ic,jc,kc)+w2*crse(ic,jc,kc+1))/(w1+w2);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_xy (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        Real w4 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        return (w1 * aa_interp_line_y(crse,sig,i-1,j  ,k,ic  ,jc  ,kc) +
                w2 * aa_interp_line_y(crse,sig,i+1,j  ,k,ic+1,jc  ,kc) +
                w3 * aa_interp_line_x(crse,sig,i  ,j-1,k,ic  ,jc  ,kc) +
                w4 * aa_interp_line_x(crse,sig,i  ,j+1,k,ic  ,jc+1,kc)) / (w1+w2+w3+w4);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_xz (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        Real w4 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sig,i-1,j,k  ,ic  ,jc,kc  ) +
                w2 * aa_interp_line_z(crse,sig,i+1,j,k  ,ic+1,jc,kc  ) +
                w3 * aa_interp_line_x(crse,sig,i  ,j,k-1,ic  ,jc,kc  ) +
                w4 * aa_interp_line_x(crse,sig,i  ,j,k+1,ic  ,jc,kc+1)) / (w1+w2+w3+w4);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_yz (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        Real w2 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        Real w4 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sig,i,j-1,k  ,ic,jc  ,kc  ) +
                w2 * aa_interp_line_z(crse,sig,i,j+1,k  ,ic,jc+1,kc  ) +
                w3 * aa_interp_line_y(crse,sig,i,j  ,k-1,ic,jc  ,kc  ) +
                w4 * aa_interp_line_y(crse,sig,i,j  ,k+1,ic,jc  ,kc+1)) / (w1+w2+w3+w4);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_interpadd_aa (int i, int j, int k, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sig,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        int kc = amrex::coarsen(k,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        bool k_is_odd = (kc*2 != k);
        if (i_is_odd and j_is_odd and k_is_odd) {
            // Fine node at center of cell
            Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
            Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
            Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
            Real w4 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
            Real w5 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
            Real w6 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
            fine(i,j,k) += (w1 * aa_interp_face_yz(crse,sig,i-1,j  ,k  ,ic  ,jc  ,kc  ) +
                            w2 * aa_interp_face_yz(crse,sig,i+1,j  ,k  ,ic+1,jc  ,kc  ) +
                            w3 * aa_interp_face_xz(crse,sig,i  ,j-1,k  ,ic  ,jc  ,kc  ) +
                            w4 * aa_interp_face_xz(crse,sig,i  ,j+1,k  ,ic  ,jc+1,kc  ) +
                            w5 * aa_interp_face_xy(crse,sig,i  ,j  ,k-1,ic  ,jc  ,kc  ) +
                            w6 * aa_interp_face_xy(crse,sig,i  ,j  ,k+1,ic  ,jc  ,kc+1))
                / (w1+w2+w3+w4+w5+w6);
        } else if (j_is_odd and k_is_odd) {
            // Node on a Y-Z face
            fine(i,j,k) += aa_interp_face_yz(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd and k_is_odd) {
            // Node on a Z-X face
            fine(i,j,k) += aa_interp_face_xz(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd and j_is_odd) {
            // Node on a X-Y face
            fine(i,j,k) += aa_interp_face_xy(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,k) += aa_interp_line_x(crse,sig,i,j,k,ic,jc,kc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,k) += aa_interp_line_y(crse,sig,i,j,k,ic,jc,kc);
        } else if (k_is_odd) {
            // Node on Z line
            fine(i,j,k) += aa_interp_line_z(crse,sig,i,j,k,ic,jc,kc);
        } else {
            // Node coincident with coarse node
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
}

namespace {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ha_interp_face_xy (Array4<Real const> const& crse,
                            Array4<Real const> const& sigx, Array4<Real const> const& sigy,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sigx(i-1,j-1,k-1) + sigx(i-1,j,k-1) + sigx(i-1,j-1,k) + sigx(i-1,j,k);
        Real w2 = sigx(i  ,j-1,k-1) + sigx(i  ,j,k-1) + sigx(i  ,j-1,k) + sigx(i  ,j,k);
        Real w3 = sigy(i-1,j-1,k-1) + sigy(i,j-1,k-1) + sigy(i-1,j-1,k) + sigy(i,j-1,k);
        Real w4 = sigy(i-1,j  ,k-1) + sigy(i,j  ,k-1) + sigy(i-1,j  ,k) + sigy(i,j  ,k);
        return (w1 * aa_interp_line_y(crse,sigy,i-1,j  ,k,ic  ,jc  ,kc) +
                w2 * aa_interp_line_y(crse,sigy,i+1,j  ,k,ic+1,jc  ,kc) +
                w3 * aa_interp_line_x(crse,sigx,i  ,j-1,k,ic  ,jc  ,kc) +
                w4 * aa_interp_line_x(crse,sigx,i  ,j+1,k,ic  ,jc+1,kc)) / (w1+w2+w3+w4);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ha_interp_face_xz (Array4<Real const> const& crse,
                            Array4<Real const> const& sigx, Array4<Real const> const& sigz,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sigx(i-1,j-1,k-1) + sigx(i-1,j,k-1) + sigx(i-1,j-1,k) + sigx(i-1,j,k);
        Real w2 = sigx(i  ,j-1,k-1) + sigx(i  ,j,k-1) + sigx(i  ,j-1,k) + sigx(i  ,j,k);
        Real w3 = sigz(i-1,j-1,k-1) + sigz(i,j-1,k-1) + sigz(i-1,j,k-1) + sigz(i,j,k-1);
        Real w4 = sigz(i-1,j-1,k  ) + sigz(i,j-1,k  ) + sigz(i-1,j,k  ) + sigz(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sigz,i-1,j,k  ,ic  ,jc,kc  ) +
                w2 * aa_interp_line_z(crse,sigz,i+1,j,k  ,ic+1,jc,kc  ) +
                w3 * aa_interp_line_x(crse,sigx,i  ,j,k-1,ic  ,jc,kc  ) +
                w4 * aa_interp_line_x(crse,sigx,i  ,j,k+1,ic  ,jc,kc+1)) / (w1+w2+w3+w4);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ha_interp_face_yz (Array4<Real const> const& crse,
                            Array4<Real const> const& sigy, Array4<Real const> const& sigz,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sigy(i-1,j-1,k-1) + sigy(i,j-1,k-1) + sigy(i-1,j-1,k) + sigy(i,j-1,k);
        Real w2 = sigy(i-1,j  ,k-1) + sigy(i,j  ,k-1) + sigy(i-1,j  ,k) + sigy(i,j  ,k);
        Real w3 = sigz(i-1,j-1,k-1) + sigz(i,j-1,k-1) + sigz(i-1,j,k-1) + sigz(i,j,k-1);
        Real w4 = sigz(i-1,j-1,k  ) + sigz(i,j-1,k  ) + sigz(i-1,j,k  ) + sigz(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sigz,i,j-1,k  ,ic,jc  ,kc  ) +
                w2 * aa_interp_line_z(crse,sigz,i,j+1,k  ,ic,jc+1,kc  ) +
                w3 * aa_interp_line_y(crse,sigy,i,j  ,k-1,ic,jc  ,kc  ) +
                w4 * aa_interp_line_y(crse,sigy,i,j  ,k+1,ic,jc  ,kc+1)) / (w1+w2+w3+w4);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_interpadd_ha (int i, int j, int k, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sigx,
                           Array4<Real const> const& sigy, Array4<Real const> const& sigz,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        int kc = amrex::coarsen(k,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        bool k_is_odd = (kc*2 != k);
        if (i_is_odd and j_is_odd and k_is_odd) {
            // Fine node at center of cell
            Real w1 = sigx(i-1,j-1,k-1) + sigx(i-1,j,k-1) + sigx(i-1,j-1,k) + sigx(i-1,j,k);
            Real w2 = sigx(i  ,j-1,k-1) + sigx(i  ,j,k-1) + sigx(i  ,j-1,k) + sigx(i  ,j,k);
            Real w3 = sigy(i-1,j-1,k-1) + sigy(i,j-1,k-1) + sigy(i-1,j-1,k) + sigy(i,j-1,k);
            Real w4 = sigy(i-1,j  ,k-1) + sigy(i,j  ,k-1) + sigy(i-1,j  ,k) + sigy(i,j  ,k);
            Real w5 = sigz(i-1,j-1,k-1) + sigz(i,j-1,k-1) + sigz(i-1,j,k-1) + sigz(i,j,k-1);
            Real w6 = sigz(i-1,j-1,k  ) + sigz(i,j-1,k  ) + sigz(i-1,j,k  ) + sigz(i,j,k  );
            fine(i,j,k) += (w1 * ha_interp_face_yz(crse,sigy,sigz,i-1,j  ,k  ,ic  ,jc  ,kc  ) +
                            w2 * ha_interp_face_yz(crse,sigy,sigz,i+1,j  ,k  ,ic+1,jc  ,kc  ) +
                            w3 * ha_interp_face_xz(crse,sigx,sigz,i  ,j-1,k  ,ic  ,jc  ,kc  ) +
                            w4 * ha_interp_face_xz(crse,sigx,sigz,i  ,j+1,k  ,ic  ,jc+1,kc  ) +
                            w5 * ha_interp_face_xy(crse,sigx,sigy,i  ,j  ,k-1,ic  ,jc  ,kc  ) +
                            w6 * ha_interp_face_xy(crse,sigx,sigy,i  ,j  ,k+1,ic  ,jc  ,kc+1))
                / (w1+w2+w3+w4+w5+w6);
        } else if (j_is_odd and k_is_odd) {
            // Node on a Y-Z face
            fine(i,j,k) += ha_interp_face_yz(crse,sigy,sigz,i,j,k,ic,jc,kc);
        } else if (i_is_odd and k_is_odd) {
            // Node on a Z-X face
            fine(i,j,k) += ha_interp_face_xz(crse,sigx,sigz,i,j,k,ic,jc,kc);
        } else if (i_is_odd and j_is_odd) {
            // Node on a X-Y face
            fine(i,j,k) += ha_interp_face_xy(crse,sigx,sigy,i,j,k,ic,jc,kc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,k) += aa_interp_line_x(crse,sigx,i,j,k,ic,jc,kc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,k) += aa_interp_line_y(crse,sigy,i,j,k,ic,jc,kc);
        } else if (k_is_odd) {
            // Node on Z line
            fine(i,j,k) += aa_interp_line_z(crse,sigz,i,j,k,ic,jc,kc);
        } else {
            // Node coincident with coarse node
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
}

//
// rhs & u
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_divu (int i, int j, int k, Array4<Real> const& rhs, Array4<Real const> const& vel,
                   Array4<int const> const& msk,
                   GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = 0.25*dxinv[0];
    Real facy = 0.25*dxinv[1];
    Real facz = 0.25*dxinv[2];

    if (msk(i,j,k)) {
        rhs(i,j,k) = 0.0;
    } else {
        rhs(i,j,k) = facx*(-vel(i-1,j-1,k-1,0)+vel(i,j-1,k-1,0)
                           -vel(i-1,j  ,k-1,0)+vel(i,j  ,k-1,0)
                           -vel(i-1,j-1,k  ,0)+vel(i,j-1,k  ,0)
                           -vel(i-1,j  ,k  ,0)+vel(i,j  ,k  ,0))
                   + facy*(-vel(i-1,j-1,k-1,1)-vel(i,j-1,k-1,1)
                           +vel(i-1,j  ,k-1,1)+vel(i,j  ,k-1,1)
                           -vel(i-1,j-1,k  ,1)-vel(i,j-1,k  ,1)
                           +vel(i-1,j  ,k  ,1)+vel(i,j  ,k  ,1))
                   + facz*(-vel(i-1,j-1,k-1,2)-vel(i,j-1,k-1,2)
                           -vel(i-1,j  ,k-1,2)-vel(i,j  ,k-1,2)
                           +vel(i-1,j-1,k  ,2)+vel(i,j-1,k  ,2)
                           +vel(i-1,j  ,k  ,2)+vel(i,j  ,k  ,2));
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
        r = 0.125 * (rhcc(i-1,j-1,k-1)+rhcc(i,j-1,k-1)+rhcc(i-1,j,k-1)+rhcc(i,j,k-1) +
                     rhcc(i-1,j-1,k  )+rhcc(i,j-1,k  )+rhcc(i-1,j,k  )+rhcc(i,j,k  ));
    }
    return r;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_mknewu (int i, int j, int k, Array4<Real> const& u, Array4<Real const> const& p,
                     Array4<Real const> const& sig, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = 0.25*dxinv[0];
    Real facy = 0.25*dxinv[1];
    Real facz = 0.25*dxinv[2];
    u(i,j,k,0) -= sig(i,j,k)*facx
        * (-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  )
           -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1));
    u(i,j,k,1) -= sig(i,j,k)*facy
        * (-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  )
           -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
    u(i,j,k,2) -= sig(i,j,k)*facz
        * (-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  )
           +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_divu_compute_fine_contrib (int i, int j, int k, Box const& fvbx,
                                        Array4<Real> const& frh, Array4<Real const> const& vel,
                                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    IntVect iv(i,j,k);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv))
    {
        frh(i,j,k) = 0.25*dxinv[0]*(-vel(i-1,j-1,k-1,0)+vel(i,j-1,k-1,0)
                                    -vel(i-1,j  ,k-1,0)+vel(i,j  ,k-1,0)
                                    -vel(i-1,j-1,k  ,0)+vel(i,j-1,k  ,0)
                                    -vel(i-1,j  ,k  ,0)+vel(i,j  ,k  ,0))
                   + 0.25*dxinv[1]*(-vel(i-1,j-1,k-1,1)-vel(i,j-1,k-1,1)
                                    +vel(i-1,j  ,k-1,1)+vel(i,j  ,k-1,1)
                                    -vel(i-1,j-1,k  ,1)-vel(i,j-1,k  ,1)
                                    +vel(i-1,j  ,k  ,1)+vel(i,j  ,k  ,1))
                   + 0.25*dxinv[2]*(-vel(i-1,j-1,k-1,2)-vel(i,j-1,k-1,2)
                                    -vel(i-1,j  ,k-1,2)-vel(i,j  ,k-1,2)
                                    +vel(i-1,j-1,k  ,2)+vel(i,j-1,k  ,2)
                                    +vel(i-1,j  ,k  ,2)+vel(i,j  ,k  ,2));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_divu_add_fine_contrib (int i, int j, int k, Box const& fvbx,
                                    Array4<Real> const& rhs, Array4<Real const> const& frh,
                                    Array4<int const> const& msk) noexcept
{
    constexpr Real rfd = 0.125_rt;
    constexpr Real chip = 0.5_rt;
    constexpr Real chip2 = 0.25_rt;
    constexpr Real chip3 = 0.125_rt;

    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    IntVect iv(ii,jj,kk);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv) and msk(ii,jj,kk))
    {
        rhs(i,j,k) +=
            rfd*(frh(ii,jj,kk)
                 + chip*(frh(ii,jj,kk-1)+frh(ii,jj,kk+1)
                         +frh(ii,jj-1,kk)+frh(ii,jj+1,kk)
                         +frh(ii-1,jj,kk)+frh(ii+1,jj,kk))
                 + chip2*(frh(ii,jj-1,kk-1)+frh(ii,jj+1,kk-1)+frh(ii,jj-1,kk+1)+frh(ii,jj+1,kk+1)
                          +frh(ii-1,jj,kk-1)+frh(ii+1,jj,kk-1)+frh(ii-1,jj,kk+1)+frh(ii+1,jj,kk+1)
                          +frh(ii-1,jj-1,kk)+frh(ii+1,jj-1,kk)+frh(ii-1,jj+1,kk)+frh(ii+1,jj+1,kk))
                 + chip3*(frh(ii-1,jj-1,kk-1)+frh(ii+1,jj-1,kk-1)
                          +frh(ii-1,jj+1,kk-1)+frh(ii+1,jj+1,kk-1)
                          +frh(ii-1,jj-1,kk+1)+frh(ii+1,jj-1,kk+1)
                          +frh(ii-1,jj+1,kk+1)+frh(ii+1,jj+1,kk+1)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_rhcc_fine_contrib (int i, int j, int k, Box const& fvbx,
                                Array4<Real> const& rhs, Array4<Real const> const& cc,
                                Array4<int const> const& msk) noexcept
{
    constexpr Real fac[] = {0.125_rt, 0.375_rt, 0.375_rt, 0.125_rt};
    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    IntVect iv(ii,jj,kk);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv) and msk(ii,jj,kk))
    {
        Real r = 0.0;
        for (int koff = -2; koff <= 1; ++koff) {
        for (int joff = -2; joff <= 1; ++joff) {
        for (int ioff = -2; ioff <= 1; ++ioff) {
            r += cc(ii+ioff,jj+joff,kk+koff) * fac[ioff+2]*fac[joff+2]*fac[koff+2];
        }}}
        rhs(i,j,k) += r;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_divu_cf_contrib (int i, int j, int k, Array4<Real> const& rhs,
                              Array4<Real const> const& vel, Array4<Real const> const& fc,
                              Array4<Real const> const& rhcc, Array4<int const> const& dmsk,
                              Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi) noexcept
{
    if (!dmsk(i,j,k) and ndmsk(i,j,k) == crse_fine_node) {
        Real facx = 0.25_rt * dxinv[0];
        Real facy = 0.25_rt * dxinv[1];
        Real facz = 0.25_rt * dxinv[2];
        const auto ndlo = amrex::lbound(nddom);
        const auto ndhi = amrex::ubound(nddom);
        Real r = fc(i,j,k);
        if (rhcc) {
            r += 0.125_rt*((1._rt-ccmsk(i-1,j-1,k-1)) * rhcc(i-1,j-1,k-1)
                         + (1._rt-ccmsk(i  ,j-1,k-1)) * rhcc(i  ,j-1,k-1)
                         + (1._rt-ccmsk(i-1,j  ,k-1)) * rhcc(i-1,j  ,k-1)
                         + (1._rt-ccmsk(i  ,j  ,k-1)) * rhcc(i  ,j  ,k-1)
                         + (1._rt-ccmsk(i-1,j-1,k  )) * rhcc(i-1,j-1,k  )
                         + (1._rt-ccmsk(i  ,j-1,k  )) * rhcc(i  ,j-1,k  )
                         + (1._rt-ccmsk(i-1,j  ,k  )) * rhcc(i-1,j  ,k  )
                         + (1._rt-ccmsk(i  ,j  ,k  )) * rhcc(i  ,j  ,k  ));
        }
        if (ccmsk(i-1,j-1,k-1) == crse_cell) {
            r += - facx*vel(i-1,j-1,k-1,0)
                 - facy*vel(i-1,j-1,k-1,1)
                 - facz*vel(i-1,j-1,k-1,2);
        }
        if (ccmsk(i,j-1,k-1) == crse_cell) {
            r += + facx*vel(i  ,j-1,k-1,0)
                 - facy*vel(i  ,j-1,k-1,1)
                 - facz*vel(i  ,j-1,k-1,2);
        }
        if (ccmsk(i-1,j,k-1) == crse_cell) {
            r += - facx*vel(i-1,j  ,k-1,0)
                 + facy*vel(i-1,j  ,k-1,1)
                 - facz*vel(i-1,j  ,k-1,2);
        }
        if (ccmsk(i,j,k-1) == crse_cell) {
            r += + facx*vel(i  ,j  ,k-1,0)
                 + facy*vel(i  ,j  ,k-1,1)
                 - facz*vel(i  ,j  ,k-1,2);
        }
        if (ccmsk(i-1,j-1,k) == crse_cell) {
            r += - facx*vel(i-1,j-1,k  ,0)
                 - facy*vel(i-1,j-1,k  ,1)
                 + facz*vel(i-1,j-1,k  ,2);
        }
        if (ccmsk(i,j-1,k) == crse_cell) {
            r += + facx*vel(i  ,j-1,k  ,0)
                 - facy*vel(i  ,j-1,k  ,1)
                 + facz*vel(i  ,j-1,k  ,2);
        }
        if (ccmsk(i-1,j,k) == crse_cell) {
            r += - facx*vel(i-1,j  ,k  ,0)
                 + facy*vel(i-1,j  ,k  ,1)
                 + facz*vel(i-1,j  ,k  ,2);
        }
        if (ccmsk(i,j,k) == crse_cell) {
            r += + facx*vel(i  ,j  ,k  ,0)
                 + facy*vel(i  ,j  ,k  ,1)
                 + facz*vel(i  ,j  ,k  ,2);
        }

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

        if (k == ndlo.z and ( bclo[2] == LinOpBCType::Neumann or
                              bclo[2] == LinOpBCType::inflow)) {
            r *= 2._rt;
        } else if (k == ndhi.z and ( bchi[2] == LinOpBCType::Neumann or
                                     bchi[2] == LinOpBCType::inflow)) {
            r *= 2._rt;
        }

        rhs(i,j,k) = r;
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
    if ((msk(i-1,j-1,k-1) == 0 or
         msk(i  ,j-1,k-1) == 0 or
         msk(i-1,j  ,k-1) == 0 or
         msk(i  ,j  ,k-1) == 0 or
         msk(i-1,j-1,k  ) == 0 or
         msk(i  ,j-1,k  ) == 0 or
         msk(i-1,j  ,k  ) == 0 or
         msk(i  ,j  ,k  ) == 0) and
        (msk(i-1,j-1,k-1) == 0 or
         msk(i  ,j-1,k-1) == 0 or
         msk(i-1,j  ,k-1) == 0 or
         msk(i  ,j  ,k-1) == 0 or
         msk(i-1,j-1,k  ) == 0 or
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
        if (k == ndlo.z and ( bclo[2] == LinOpBCType::Neumann or
                              bclo[2] == LinOpBCType::inflow)) {
            fac *= 2._rt;
        } else if (k == ndhi.z and ( bchi[2] == LinOpBCType::Neumann or
                                     bchi[2] == LinOpBCType::inflow)) {
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
void mlndhelm_res_fine_Ax (int i, int j, int k, Box const& fvbx, Array4<Real> const& Ax,
                          Array4<Real const> const& x, Array4<Real const> const& sig,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    IntVect iv(i,j,k);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv)) {
        Real facx = (1._rt/36._rt)*dxinv[0]*dxinv[0];
        Real facy = (1._rt/36._rt)*dxinv[1]*dxinv[1];
        Real facz = (1._rt/36._rt)*dxinv[2]*dxinv[2];
        Real fxyz = facx + facy + facz;
        Real fmx2y2z = -facx + 2._rt*facy + 2._rt*facz;
        Real f2xmy2z = 2._rt*facx - facy + 2._rt*facz;
        Real f2x2ymz = 2._rt*facx + 2._rt*facy - facz;
        Real f4xm2ym2z = 4._rt*facx - 2._rt*facy - 2._rt*facz;
        Real fm2x4ym2z = -2._rt*facx + 4._rt*facy - 2._rt*facz;
        Real fm2xm2y4z = -2._rt*facx - 2._rt*facy + 4._rt*facz;
        Ax(i,j,k) = x(i,j,k)*(-4._rt)*fxyz*
            (sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
            +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ))
            //
            + fxyz*(x(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                  + x(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                  + x(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                  + x(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                  + x(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                  + x(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                  + x(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                  + x(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
            + fmx2y2z*(x(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                     + x(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                     + x(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                     + x(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
            + f2xmy2z*(x(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                     + x(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                     + x(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                     + x(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
            + f2x2ymz*(x(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                     + x(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                     + x(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                     + x(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
            + f4xm2ym2z*(x(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                       + x(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
            + fm2x4ym2z*(x(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                       + x(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
            + fm2xm2y4z*(x(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                       + x(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_res_fine_contrib (int i, int j, int k, Array4<Real> const& f,
                               Array4<Real const> const& Ax, Array4<int const> const& msk) noexcept
{
    constexpr Real rfd = 0.125_rt;
    constexpr Real chip = 0.5_rt;
    constexpr Real chip2 = 0.25_rt;
    constexpr Real chip3 = 0.125_rt;

    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    if (msk(ii,jj,kk)) {
        f(i,j,k) += rfd*
            (Ax(ii,jj,kk)
             + chip*(Ax(ii,jj,kk-1)+Ax(ii,jj,kk+1)
                     +Ax(ii,jj-1,kk)+Ax(ii,jj+1,kk)
                     +Ax(ii-1,jj,kk)+Ax(ii+1,jj,kk))
             + chip2*(Ax(ii,jj-1,kk-1)+Ax(ii,jj+1,kk-1)+Ax(ii,jj-1,kk+1)+Ax(ii,jj+1,kk+1)
                      +Ax(ii-1,jj,kk-1)+Ax(ii+1,jj,kk-1)+Ax(ii-1,jj,kk+1)+Ax(ii+1,jj,kk+1)
                      +Ax(ii-1,jj-1,kk)+Ax(ii+1,jj-1,kk)+Ax(ii-1,jj+1,kk)+Ax(ii+1,jj+1,kk))
             + chip3*(Ax(ii-1,jj-1,kk-1)+Ax(ii+1,jj-1,kk-1)
                      +Ax(ii-1,jj+1,kk-1)+Ax(ii+1,jj+1,kk-1)
                      +Ax(ii-1,jj-1,kk+1)+Ax(ii+1,jj-1,kk+1)
                      +Ax(ii-1,jj+1,kk+1)+Ax(ii+1,jj+1,kk+1)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_res_cf_contrib (int i, int j, int k, Array4<Real> const& res,
                             Array4<Real const> const& phi, Array4<Real const> const& rhs,
                             Array4<Real const> const& sig, Array4<int const> const& dmsk,
                             Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                             Array4<Real const> const& fc,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv, Box const& nddom,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    if (!dmsk(i,j,k) and ndmsk(i,j,k) == crse_fine_node) {

        Real facx = (1._rt/36._rt)*dxinv[0]*dxinv[0];
        Real facy = (1._rt/36._rt)*dxinv[1]*dxinv[1];
        Real facz = (1._rt/36._rt)*dxinv[2]*dxinv[2];

        Real Ax = 0._rt;
        if (ccmsk(i-1,j-1,k-1) == crse_cell) {
            Ax += sig(i-1,j-1,k-1)*(facx*(4._rt*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                         +2._rt*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  ))
                                         +2._rt*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1))
                                          +     (phi(i-1,j-1,k-1)-phi(i  ,j-1,k-1)))
                                  + facy*(4._rt*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                         +2._rt*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  ))
                                         +2._rt*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1))
                                          +     (phi(i-1,j-1,k-1)-phi(i-1,j  ,k-1)))
                                  + facz*(4._rt*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                         +2._rt*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  ))
                                         +2._rt*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  ))
                                          +     (phi(i-1,j-1,k-1)-phi(i-1,j-1,k  ))));
        }
        if (ccmsk(i,j-1,k-1) == crse_cell) {
            Ax += sig(i,j-1,k-1)*(facx*(4._rt*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                       +2._rt*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  ))
                                       +2._rt*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1))
                                        +     (phi(i+1,j-1,k-1)-phi(i  ,j-1,k-1)))
                                + facy*(4._rt*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                       +2._rt*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  ))
                                       +2._rt*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1))
                                        +     (phi(i+1,j-1,k-1)-phi(i+1,j  ,k-1)))
                                + facz*(4._rt*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                       +2._rt*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  ))
                                       +2._rt*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  ))
                                        +     (phi(i+1,j-1,k-1)-phi(i+1,j-1,k  ))));
        }
        if (ccmsk(i-1,j,k-1) == crse_cell) {
            Ax += sig(i-1,j,k-1)*(facx*(4._rt*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                       +2._rt*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  ))
                                       +2._rt*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1))
                                        +     (phi(i-1,j+1,k-1)-phi(i  ,j+1,k-1)))
                                 + facy*(4._rt*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                        +2._rt*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  ))
                                        +2._rt*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1))
                                         +     (phi(i-1,j+1,k-1)-phi(i-1,j  ,k-1)))
                                 + facz*(4._rt*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                        +2._rt*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  ))
                                        +2._rt*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  ))
                                         +     (phi(i-1,j+1,k-1)-phi(i-1,j+1,k  ))));
        }
        if (ccmsk(i,j,k-1) == crse_cell) {
            Ax += sig(i,j,k-1)*(facx*(4._rt*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +2._rt*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  ))
                                     +2._rt*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1))
                                      +     (phi(i+1,j+1,k-1)-phi(i  ,j+1,k-1)))
                              + facy*(4._rt*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                     +2._rt*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  ))
                                     +2._rt*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1))
                                      +     (phi(i+1,j+1,k-1)-phi(i+1,j  ,k-1)))
                              + facz*(4._rt*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  ))
                                     +2._rt*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  ))
                                     +2._rt*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  ))
                                      +     (phi(i+1,j+1,k-1)-phi(i+1,j+1,k  ))));
        }
        if (ccmsk(i-1,j-1,k) == crse_cell) {
            Ax += sig(i-1,j-1,k)*(facx*(4._rt*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                       +2._rt*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  ))
                                       +2._rt*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1))
                                        +     (phi(i-1,j-1,k+1)-phi(i  ,j-1,k+1)))
                                + facy*(4._rt*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                       +2._rt*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  ))
                                       +2._rt*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1))
                                        +     (phi(i-1,j-1,k+1)-phi(i-1,j  ,k+1)))
                                + facz*(4._rt*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                       +2._rt*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  ))
                                       +2._rt*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  ))
                                        +     (phi(i-1,j-1,k+1)-phi(i-1,j-1,k  ))));
        }
        if (ccmsk(i,j-1,k) == crse_cell) {
            Ax += sig(i,j-1,k)*(facx*(4._rt*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +2._rt*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  ))
                                     +2._rt*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1))
                                      +     (phi(i+1,j-1,k+1)-phi(i  ,j-1,k+1)))
                              + facy*(4._rt*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  ))
                                     +2._rt*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  ))
                                     +2._rt*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1))
                                      +     (phi(i+1,j-1,k+1)-phi(i+1,j  ,k+1)))
                              + facz*(4._rt*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                     +2._rt*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  ))
                                     +2._rt*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  ))
                                      +     (phi(i+1,j-1,k+1)-phi(i+1,j-1,k  ))));
        }
        if (ccmsk(i-1,j,k) == crse_cell) {
            Ax += sig(i-1,j,k)*(facx*(4._rt*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  ))
                                     +2._rt*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  ))
                                     +2._rt*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1))
                                      +     (phi(i-1,j+1,k+1)-phi(i  ,j+1,k+1)))
                              + facy*(4._rt*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                     +2._rt*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  ))
                                     +2._rt*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1))
                                      +     (phi(i-1,j+1,k+1)-phi(i-1,j  ,k+1)))
                              + facz*(4._rt*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                     +2._rt*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  ))
                                     +2._rt*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  ))
                                      +     (phi(i-1,j+1,k+1)-phi(i-1,j+1,k  ))));
        }
        if (ccmsk(i,j,k) == crse_cell) {
            Ax += sig(i,j,k)*(facx*(4._rt*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  ))
                                   +2._rt*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  ))
                                   +2._rt*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1))
                                    +     (phi(i+1,j+1,k+1)-phi(i  ,j+1,k+1)))
                            + facy*(4._rt*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  ))
                                   +2._rt*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  ))
                                   +2._rt*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1))
                                    +     (phi(i+1,j+1,k+1)-phi(i+1,j  ,k+1)))
                            + facz*(4._rt*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  ))
                                   +2._rt*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  ))
                                   +2._rt*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  ))
                                    +     (phi(i+1,j+1,k+1)-phi(i+1,j+1,k  ))));
        }

        Real Axf = fc(i,j,k);
        const auto ndlo = amrex::lbound(nddom);
        const auto ndhi = amrex::ubound(nddom);

        if (i == ndlo.x and (bclo[0] == LinOpBCType::Neumann or
                             bclo[0] == LinOpBCType::inflow)) {
            Axf *= 2._rt;
        } else if (i== ndhi.x and (bchi[0] == LinOpBCType::Neumann or
                                   bchi[0] == LinOpBCType::inflow)) {
            Axf *= 2._rt;
        }

        if (j == ndlo.y and ( bclo[1] == LinOpBCType::Neumann or
                              bclo[1] == LinOpBCType::inflow)) {
            Axf *= 2._rt;
        } else if (j == ndhi.y and (bchi[1] == LinOpBCType::Neumann or
                                    bchi[1] == LinOpBCType::inflow)) {
            Axf *= 2._rt;
        }

        if (k == ndlo.z and (bclo[2] == LinOpBCType::Neumann or
                             bclo[2] == LinOpBCType::inflow)) {
            Axf *= 2._rt;
        } else if (k == ndhi.z and (bchi[2] == LinOpBCType::Neumann or
                                    bchi[2] == LinOpBCType::inflow)) {
            Axf *= 2._rt;
        }

        res(i,j,k) = rhs(i,j,k) - (Ax + Axf);
    }
}

//
// RAP
//

namespace {

    constexpr int ist_000 = 0;
    constexpr int ist_p00 = 1;
    constexpr int ist_0p0 = 2;
    constexpr int ist_00p = 3;
    constexpr int ist_pp0 = 4;
    constexpr int ist_p0p = 5;
    constexpr int ist_0pp = 6;
    constexpr int ist_ppp = 7;
    constexpr int ist_inv = 8;
    constexpr int n_sten  = 9;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_set_stencil (Box const& bx, Array4<Real> const& sten,
                          Array4<Real const> const& sig,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = (1.0/36.0)*dxinv[2]*dxinv[2];
    Real fxyz = facx + facy + facz;
    Real fmx2y2z = -facx + 2.0*facy + 2.0*facz;
    Real f2xmy2z = 2.0*facx - facy + 2.0*facz;
    Real f2x2ymz = 2.0*facx + 2.0*facy - facz;
    Real f4xm2ym2z = 4.0*facx - 2.0*facy - 2.0*facz;
    Real fm2x4ym2z = -2.0*facx + 4.0*facy - 2.0*facz;
    Real fm2xm2y4z = -2.0*facx - 2.0*facy + 4.0*facz;

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        // i+1,j,k
        sten(i,j,k,ist_p00) = f4xm2ym2z * (sig(i,j-1,k-1)+sig(i,j,k-1)+sig(i,j-1,k)+sig(i,j,k));
        // i-1,j,k: sten(i-1,j,k,ist_p00)

        // i,j+1,k
        sten(i,j,k,ist_0p0) = fm2x4ym2z * (sig(i-1,j,k-1)+sig(i,j,k-1)+sig(i-1,j,k)+sig(i,j,k));
        // i,j-1,k: sten(i,j-1,k,ist_0p0)

        // i,j,k+1
        sten(i,j,k,ist_00p) = fm2xm2y4z * (sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
        // i,j,k-1: sten(i,j,k-1,ist_00p)

        // i+1,j+1,k
        sten(i,j,k,ist_pp0) = f2x2ymz * (sig(i,j,k-1)+sig(i,j,k));
        // i-1,j-1,k: sten(i-1,j-1,k,ist_pp0)
        // i+1,j-1,k: sten(i  ,j-1,k,ist_pp0)
        // i-1,j+1,k: sten(i-1,j  ,k,ist_pp0)

        // i+1,j,k+1
        sten(i,j,k,ist_p0p) = f2xmy2z * (sig(i,j-1,k)+sig(i,j,k));
        // i-1,j,k-1: sten(i-1,j,k-1,ist_p0p)
        // i+1,j,k-1: sten(i  ,j,k-1,ist_p0p)
        // i-1,j,k+1: sten(i-1,j,k  ,ist_p0p)

        // i,j+1,k+1
        sten(i,j,k,ist_0pp) = fmx2y2z * (sig(i-1,j,k)+sig(i,j,k));
        // i,j-1,k-1: sten(i,j-1,k-1,ist_0pp)
        // i,j+1,k-1: sten(i,j  ,k-1,ist_0pp)
        // i,j-1,k+1: sten(i,j-1,k  ,ist_0pp)

        // i+1,j+1,k+1
        sten(i,j,k,ist_ppp) = fxyz * sig(i,j,k);
        // i-1,j-1,k-1: sten(i-1,j-1,k-1,ist_ppp)
        // i+1,j-1,k-1: sten(i  ,j-1,k-1,ist_ppp)
        // i-1,j+1,k-1: sten(i-1,j  ,k-1,ist_ppp)
        // i+1,j+1,k-1: sten(i  ,j  ,k-1,ist_ppp)
        // i-1,j-1,k+1: sten(i-1,j-1,k  ,ist_ppp)
        // i+1,j-1,k+1: sten(i  ,j-1,k  ,ist_ppp)
        // i-1,j+1,k+1: sten(i-1,j  ,k  ,ist_ppp)
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_set_stencil_s0 (int i, int j, int k, Array4<Real> const& sten) noexcept
{
    sten(i,j,k,ist_000) = -(sten(i-1,j,k,ist_p00) + sten(i,j,k,ist_p00)
                          + sten(i,j-1,k,ist_0p0) + sten(i,j,k,ist_0p0)
                          + sten(i,j,k-1,ist_00p) + sten(i,j,k,ist_00p)
                          + sten(i-1,j-1,k,ist_pp0) + sten(i,j-1,k,ist_pp0)
                          + sten(i-1,j,k,ist_pp0) + sten(i,j,k,ist_pp0)
                          + sten(i-1,j,k-1,ist_p0p) + sten(i,j,k-1,ist_p0p)
                          + sten(i-1,j,k,ist_p0p) + sten(i,j,k,ist_p0p)
                          + sten(i,j-1,k-1,ist_0pp) + sten(i,j,k-1,ist_0pp)
                          + sten(i,j-1,k,ist_0pp) + sten(i,j,k,ist_0pp)
                          + sten(i-1,j-1,k-1,ist_ppp) + sten(i,j-1,k-1,ist_ppp)
                          + sten(i-1,j,k-1,ist_ppp) + sten(i,j,k-1,ist_ppp)
                          + sten(i-1,j-1,k,ist_ppp) + sten(i,j-1,k,ist_ppp)
                          + sten(i-1,j,k,ist_ppp) + sten(i,j,k,ist_ppp));
    sten(i,j,k,ist_inv) = 1.0 /
        (  std::abs(sten(i-1,j,k,ist_p00)) + std::abs(sten(i,j,k,ist_p00))
         + std::abs(sten(i,j-1,k,ist_0p0)) + std::abs(sten(i,j,k,ist_0p0))
         + std::abs(sten(i,j,k-1,ist_00p)) + std::abs(sten(i,j,k,ist_00p))
         + std::abs(sten(i-1,j-1,k,ist_pp0)) + std::abs(sten(i,j-1,k,ist_pp0))
         + std::abs(sten(i-1,j,k,ist_pp0)) + std::abs(sten(i,j,k,ist_pp0))
         + std::abs(sten(i-1,j,k-1,ist_p0p)) + std::abs(sten(i,j,k-1,ist_p0p))
         + std::abs(sten(i-1,j,k,ist_p0p)) + std::abs(sten(i,j,k,ist_p0p))
         + std::abs(sten(i,j-1,k-1,ist_0pp)) + std::abs(sten(i,j,k-1,ist_0pp))
         + std::abs(sten(i,j-1,k,ist_0pp)) + std::abs(sten(i,j,k,ist_0pp))
         + std::abs(sten(i-1,j-1,k-1,ist_ppp)) + std::abs(sten(i,j-1,k-1,ist_ppp))
         + std::abs(sten(i-1,j,k-1,ist_ppp)) + std::abs(sten(i,j,k-1,ist_ppp))
         + std::abs(sten(i-1,j-1,k,ist_ppp)) + std::abs(sten(i,j-1,k,ist_ppp))
         + std::abs(sten(i-1,j,k,ist_ppp)) + std::abs(sten(i,j,k,ist_ppp)) + eps);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_stencil_rap (int i, int j, int k, Array4<Real> const& csten,
                          Array4<Real const> const& fsten) noexcept
{
    auto interp_from_mmm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = 1._rt;
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp)) + eps);
        p *= std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_pmm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = 1._rt;
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp)) + eps);
        p *= std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_mpm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = 1._rt;
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p *= std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_ppm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = 1._rt;
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) + eps);
        p *= std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_mmp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = 1._rt;
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_,ist_ppp)) + eps);
        p += std::abs(fsten(i_-1,j_-1,k_,ist_pp0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_-1,j_  ,k_,ist_p0p)) /
           ( std::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_-1,k_,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_,ist_ppp)) + eps);
        p *= std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_pmp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = 1._rt;
        p += std::abs(fsten(i_  ,j_  ,k_,ist_p00)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_-1,k_,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_-1,k_,ist_pp0)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_,ist_p0p)) /
           ( std::abs(fsten(i_  ,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_-1,k_,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_-1,k_,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_,ist_ppp)) + eps);
        p *= std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp)) * fsten(i_,j_  ,k_,ist_inv);
        return p;
    };

    auto interp_from_mpp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = 1._rt;
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_-1,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) /
           ( std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p *= std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_ppp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real p = 1._rt;
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) /
           ( std::abs(fsten(i_  ,j_-1,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) /
           ( std::abs(fsten(i_-1,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) /
           ( std::abs(fsten(i_-1,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) /
           ( std::abs(fsten(i_  ,j_  ,k_-1,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) /
           ( std::abs(fsten(i_  ,j_-1,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p += std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) /
           ( std::abs(fsten(i_-1,j_  ,k_  ,ist_ppp))
           + std::abs(fsten(i_  ,j_  ,k_  ,ist_ppp)) + eps);
        p *= std::abs(fsten(i_  ,j_  ,k_,ist_ppp)) * fsten(i_,j_,k_,ist_inv);
        return p;
    };

    auto interp_from_0mm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real wmm = std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (1._rt + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_0mp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real wmm = std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (1._rt + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_0pm_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real wmm = std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (1._rt + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_0pp_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) + eps);
        Real wmm = std::abs(fsten(i_  ,j_-1,k_-1,ist_0pp)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_0pp)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_  ,j_-1,k_  ,ist_0pp)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_0pp)) * (1._rt + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_m0m_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real wmm = std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (1._rt + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_p0m_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real wmm = std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (1._rt + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_m0p_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real wmm = std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (1._rt + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_p0p_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) + eps);
        Real wmm = std::abs(fsten(i_-1,j_  ,k_-1,ist_p0p)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_  ,k_-1,ist_p0p)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_p0p)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_p0p)) * (1._rt + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_mm0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps);
      Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                            +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real wmm = std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (1._rt + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_mp0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real wmm = std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (1._rt + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_pm0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real wmm = std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (1._rt + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_pp0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1m = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) + eps);
        Real w1p = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00)) / (std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real w2m = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) + eps);
        Real w2p = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0)) / (std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0))
                                                             +std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) + eps);
        Real wmm = std::abs(fsten(i_-1,j_-1,k_  ,ist_pp0)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(i_  ,j_-1,k_  ,ist_pp0)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(i_-1,j_  ,k_  ,ist_pp0)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(i_  ,j_  ,k_  ,ist_pp0)) * (1._rt + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto interp_from_00m_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p));
        Real w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w1 / (w1+w2);
        }
    };

    auto interp_from_00p_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = std::abs(fsten(i_  ,j_  ,k_-1,ist_00p));
        Real w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_00p));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w2 / (w1+w2);
        }
    };

    auto interp_from_0m0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0));
        Real w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w1 / (w1+w2);
        }
    };

    auto interp_from_0p0_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = std::abs(fsten(i_  ,j_-1,k_  ,ist_0p0));
        Real w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_0p0));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w2 / (w1+w2);
        }
    };

    auto interp_from_m00_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00));
        Real w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w1 / (w1+w2);
        }
    };

    auto interp_from_p00_to = [&fsten] (int i_, int j_, int k_) -> Real {
        Real w1 = std::abs(fsten(i_-1,j_  ,k_  ,ist_p00));
        Real w2 = std::abs(fsten(i_  ,j_  ,k_  ,ist_p00));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w2 / (w1+w2);
        }
    };

    auto Ammm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_-1,k_-1,ist_ppp);
    };

    auto A0mm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_-1,ist_0pp);
    };

    auto Apmm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_-1,ist_ppp);
    };

    auto Am0m = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_-1,ist_p0p);
    };

    auto A00m = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_-1,ist_00p);
    };

    auto Ap0m = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_-1,ist_p0p);
    };

    auto Ampm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_-1,ist_ppp);
    };

    auto A0pm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_-1,ist_0pp);
    };

    auto Appm = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_-1,ist_ppp);
    };

    auto Amm0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_-1,k_  ,ist_pp0);
    };

    auto A0m0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_  ,ist_0p0);
    };

    auto Apm0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_  ,ist_pp0);
    };

    auto Am00 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_  ,ist_p00);
    };

    auto A000 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_000);
    };

    auto Ap00 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_p00);
    };

    auto Amp0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_  ,ist_pp0);
    };

    auto A0p0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_0p0);
    };

    auto App0 = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_pp0);
    };

    auto Ammp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_-1,k_  ,ist_ppp);
    };

    auto A0mp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_  ,ist_0pp);
    };

    auto Apmp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_-1,k_  ,ist_ppp);
    };

    auto Am0p = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_  ,ist_p0p);
    };

    auto A00p = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_00p);
    };

    auto Ap0p = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_p0p);
    };

    auto Ampp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_-1,j_  ,k_  ,ist_ppp);
    };

    auto A0pp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_0pp);
    };

    auto Appp = [&fsten] (int i_, int j_, int k_) -> Real {
        return fsten(i_  ,j_  ,k_  ,ist_ppp);
    };

    auto restrict_from_mmm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = 1._rt;
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_p00)) /
           ( std::abs(fsten(ii_-1,jj_-2,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_0p0)) /
           ( std::abs(fsten(ii_-2,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-2,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_00p)) /
           ( std::abs(fsten(ii_-2,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-2,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_pp0)) /
           ( std::abs(fsten(ii_-1,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_p0p)) /
           ( std::abs(fsten(ii_-1,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_0pp)) /
           ( std::abs(fsten(ii_-2,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) + eps);
        r *= std::abs(fsten(ii_-1,jj_-1,kk_-1,ist_ppp)) * fsten(ii_-1,jj_-1,kk_-1,ist_inv);
        return r;
    };

    auto restrict_from_0mm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_,jj_-2,kk_-1,ist_0p0)) / (std::abs(fsten(ii_,jj_-2,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-2,kk_-1,ist_0pp)) + eps);
        Real w1p = std::abs(fsten(ii_,jj_-1,kk_-1,ist_0p0)) / (std::abs(fsten(ii_,jj_-1,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_-1,ist_0pp)) + eps);
        Real w2m = std::abs(fsten(ii_,jj_-1,kk_-2,ist_00p)) / (std::abs(fsten(ii_,jj_-2,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_-2,ist_0pp)) + eps);
        Real w2p = std::abs(fsten(ii_,jj_-1,kk_-1,ist_00p)) / (std::abs(fsten(ii_,jj_-2,kk_-1,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_-1,ist_0pp)) + eps);
        Real wmm = std::abs(fsten(ii_,jj_-2,kk_-2,ist_0pp)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_,jj_-1,kk_-2,ist_0pp)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_,jj_-2,kk_-1,ist_0pp)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_,jj_-1,kk_-1,ist_0pp)) * (1._rt + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_pmm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = 1._rt;
        r += std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_p00)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_0p0)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_00p)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_pp0)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_p0p)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_0pp)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_-1,ist_ppp)) + eps);
        r *= std::abs(fsten(ii_  ,jj_-1,kk_-1,ist_ppp)) * fsten(ii_+1,jj_-1,kk_-1,ist_inv);
        return r;
    };

    auto restrict_from_m0m_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_-2,jj_,kk_-1,ist_p00)) / (std::abs(fsten(ii_-2,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_-2,jj_,kk_-1,ist_p0p)) + eps);
        Real w1p = std::abs(fsten(ii_-1,jj_,kk_-1,ist_p00)) / (std::abs(fsten(ii_-1,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_-1,ist_p0p)) + eps);
        Real w2m = std::abs(fsten(ii_-1,jj_,kk_-2,ist_00p)) / (std::abs(fsten(ii_-2,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_-2,ist_p0p)) + eps);
        Real w2p = std::abs(fsten(ii_-1,jj_,kk_-1,ist_00p)) / (std::abs(fsten(ii_-2,jj_,kk_-1,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_-1,ist_p0p)) + eps);
        Real wmm = std::abs(fsten(ii_-2,jj_,kk_-2,ist_p0p)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_-1,jj_,kk_-2,ist_p0p)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_-2,jj_,kk_-1,ist_p0p)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_-1,jj_,kk_-1,ist_p0p)) * (1._rt + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_00m_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = std::abs(fsten(ii_,jj_,kk_-2,ist_00p));
        Real w2 = std::abs(fsten(ii_,jj_,kk_-1,ist_00p));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w2 / (w1+w2);
        }
    };

    auto restrict_from_p0m_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_  ,jj_,kk_-1,ist_p00)) / (std::abs(fsten(ii_  ,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_  ,jj_,kk_-1,ist_p0p)) + eps);
        Real w1p = std::abs(fsten(ii_+1,jj_,kk_-1,ist_p00)) / (std::abs(fsten(ii_+1,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_-1,ist_p0p)) + eps);
        Real w2m = std::abs(fsten(ii_+1,jj_,kk_-2,ist_00p)) / (std::abs(fsten(ii_  ,jj_,kk_-2,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_-2,ist_p0p)) + eps);
        Real w2p = std::abs(fsten(ii_+1,jj_,kk_-1,ist_00p)) / (std::abs(fsten(ii_  ,jj_,kk_-1,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_-1,ist_p0p)) + eps);
        Real wmm = std::abs(fsten(ii_  ,jj_,kk_-2,ist_p0p)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_+1,jj_,kk_-2,ist_p0p)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_  ,jj_,kk_-1,ist_p0p)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_+1,jj_,kk_-1,ist_p0p)) * (1._rt + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_mpm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = 1._rt;
        r += std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_p00)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_0p0)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-2,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_00p)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-2,jj_+1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_pp0)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_p0p)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_0pp)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) + eps);
        r *= std::abs(fsten(ii_-1,jj_  ,kk_-1,ist_ppp)) * fsten(ii_-1,jj_+1,kk_-1,ist_inv);
        return r;
    };

    auto restrict_from_0pm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_,jj_  ,kk_-1,ist_0p0)) / (std::abs(fsten(ii_,jj_  ,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_  ,kk_-1,ist_0pp)) + eps);
        Real w1p = std::abs(fsten(ii_,jj_+1,kk_-1,ist_0p0)) / (std::abs(fsten(ii_,jj_+1,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_-1,ist_0pp)) + eps);
        Real w2m = std::abs(fsten(ii_,jj_+1,kk_-2,ist_00p)) / (std::abs(fsten(ii_,jj_  ,kk_-2,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_-2,ist_0pp)) + eps);
        Real w2p = std::abs(fsten(ii_,jj_+1,kk_-1,ist_00p)) / (std::abs(fsten(ii_,jj_  ,kk_-1,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_-1,ist_0pp)) + eps);
        Real wmm = std::abs(fsten(ii_,jj_  ,kk_-2,ist_0pp)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_,jj_+1,kk_-2,ist_0pp)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_,jj_  ,kk_-1,ist_0pp)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_,jj_+1,kk_-1,ist_0pp)) * (1._rt + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_ppm_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = 1._rt;
        r += std::abs(fsten(ii_  ,jj_+1,kk_-1,ist_p00)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_  ,kk_-1,ist_0p0)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_+1,kk_-1,ist_00p)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_+1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_pp0)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-2,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_  ,jj_+1,kk_-1,ist_p0p)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_-1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_  ,kk_-1,ist_0pp)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_-1,ist_ppp)) + eps);
        r *= std::abs(fsten(ii_  ,jj_  ,kk_-1,ist_ppp)) * fsten(ii_+1,jj_+1,kk_-1,ist_inv);
        return r;
    };

    auto restrict_from_mm0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_-2,jj_-1,kk_,ist_p00)) / (std::abs(fsten(ii_-2,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-2,jj_-1,kk_,ist_pp0)) + eps);
        Real w1p = std::abs(fsten(ii_-1,jj_-1,kk_,ist_p00)) / (std::abs(fsten(ii_-1,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_-1,kk_,ist_pp0)) + eps);
        Real w2m = std::abs(fsten(ii_-1,jj_-2,kk_,ist_0p0)) / (std::abs(fsten(ii_-2,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_-2,kk_,ist_pp0)) + eps);
        Real w2p = std::abs(fsten(ii_-1,jj_-1,kk_,ist_0p0)) / (std::abs(fsten(ii_-2,jj_-1,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_-1,kk_,ist_pp0)) + eps);
        Real wmm = std::abs(fsten(ii_-2,jj_-2,kk_,ist_pp0)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_-1,jj_-2,kk_,ist_pp0)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_-2,jj_-1,kk_,ist_pp0)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_-1,jj_-1,kk_,ist_pp0)) * (1._rt + w1p + w2p);
        return wpp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_0m0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = std::abs(fsten(ii_,jj_-2,kk_,ist_0p0));
        Real w2 = std::abs(fsten(ii_,jj_-1,kk_,ist_0p0));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w2 / (w1+w2);
        }
    };

    auto restrict_from_pm0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_  ,jj_-1,kk_,ist_p00)) / (std::abs(fsten(ii_  ,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_  ,jj_-1,kk_,ist_pp0)) + eps);
        Real w1p = std::abs(fsten(ii_+1,jj_-1,kk_,ist_p00)) / (std::abs(fsten(ii_+1,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_-1,kk_,ist_pp0)) + eps);
        Real w2m = std::abs(fsten(ii_+1,jj_-2,kk_,ist_0p0)) / (std::abs(fsten(ii_  ,jj_-2,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_-2,kk_,ist_pp0)) + eps);
        Real w2p = std::abs(fsten(ii_+1,jj_-1,kk_,ist_0p0)) / (std::abs(fsten(ii_  ,jj_-1,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_-1,kk_,ist_pp0)) + eps);
        Real wmm = std::abs(fsten(ii_  ,jj_-2,kk_,ist_pp0)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_+1,jj_-2,kk_,ist_pp0)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_  ,jj_-1,kk_,ist_pp0)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_+1,jj_-1,kk_,ist_pp0)) * (1._rt + w1p + w2p);
        return wmp / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_m00_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = std::abs(fsten(ii_-2,jj_,kk_,ist_p00));
        Real w2 = std::abs(fsten(ii_-1,jj_,kk_,ist_p00));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w2 / (w1+w2);
        }
    };

    auto restrict_from_000_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        return 1._rt;
    };

    auto restrict_from_p00_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = std::abs(fsten(ii_  ,jj_,kk_,ist_p00));
        Real w2 = std::abs(fsten(ii_+1,jj_,kk_,ist_p00));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w1 / (w1+w2);
        }
    };

    auto restrict_from_mp0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_-2,jj_+1,kk_,ist_p00)) / (std::abs(fsten(ii_-2,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-2,jj_+1,kk_,ist_pp0)) + eps);
        Real w1p = std::abs(fsten(ii_-1,jj_+1,kk_,ist_p00)) / (std::abs(fsten(ii_-1,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_+1,kk_,ist_pp0)) + eps);
        Real w2m = std::abs(fsten(ii_-1,jj_  ,kk_,ist_0p0)) / (std::abs(fsten(ii_-2,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_  ,kk_,ist_pp0)) + eps);
        Real w2p = std::abs(fsten(ii_-1,jj_+1,kk_,ist_0p0)) / (std::abs(fsten(ii_-2,jj_+1,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_-1,jj_+1,kk_,ist_pp0)) + eps);
        Real wmm = std::abs(fsten(ii_-2,jj_  ,kk_,ist_pp0)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_-1,jj_  ,kk_,ist_pp0)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_-2,jj_+1,kk_,ist_pp0)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_-1,jj_+1,kk_,ist_pp0)) * (1._rt + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_0p0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = std::abs(fsten(ii_,jj_  ,kk_,ist_0p0));
        Real w2 = std::abs(fsten(ii_,jj_+1,kk_,ist_0p0));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w1 / (w1+w2);
        }
    };

    auto restrict_from_pp0_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_  ,jj_+1,kk_,ist_p00)) / (std::abs(fsten(ii_  ,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_  ,jj_+1,kk_,ist_pp0)) + eps);
        Real w1p = std::abs(fsten(ii_+1,jj_+1,kk_,ist_p00)) / (std::abs(fsten(ii_+1,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_+1,kk_,ist_pp0)) + eps);
        Real w2m = std::abs(fsten(ii_+1,jj_  ,kk_,ist_0p0)) / (std::abs(fsten(ii_  ,jj_  ,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_  ,kk_,ist_pp0)) + eps);
        Real w2p = std::abs(fsten(ii_+1,jj_+1,kk_,ist_0p0)) / (std::abs(fsten(ii_  ,jj_+1,kk_,ist_pp0))
                                                              +std::abs(fsten(ii_+1,jj_+1,kk_,ist_pp0)) + eps);
        Real wmm = std::abs(fsten(ii_  ,jj_  ,kk_,ist_pp0)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_+1,jj_  ,kk_,ist_pp0)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_  ,jj_+1,kk_,ist_pp0)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_+1,jj_+1,kk_,ist_pp0)) * (1._rt + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_mmp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = 1._rt;
        r += std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_p00)) /
           ( std::abs(fsten(ii_-1,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-2,kk_+1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_0p0)) /
           ( std::abs(fsten(ii_-2,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-2,jj_-1,kk_+1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_00p)) /
           ( std::abs(fsten(ii_-2,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-2,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_pp0)) /
           ( std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_p0p)) /
           ( std::abs(fsten(ii_-1,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_0pp)) /
           ( std::abs(fsten(ii_-2,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) + eps);
        r *= std::abs(fsten(ii_-1,jj_-1,kk_  ,ist_ppp)) * fsten(ii_-1,jj_-1,kk_+1,ist_inv);
        return r;
    };

    auto restrict_from_0mp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_,jj_-2,kk_+1,ist_0p0)) / (std::abs(fsten(ii_,jj_-2,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-2,kk_+1,ist_0pp)) + eps);
        Real w1p = std::abs(fsten(ii_,jj_-1,kk_+1,ist_0p0)) / (std::abs(fsten(ii_,jj_-1,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_+1,ist_0pp)) + eps);
        Real w2m = std::abs(fsten(ii_,jj_-1,kk_  ,ist_00p)) / (std::abs(fsten(ii_,jj_-2,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_  ,ist_0pp)) + eps);
        Real w2p = std::abs(fsten(ii_,jj_-1,kk_+1,ist_00p)) / (std::abs(fsten(ii_,jj_-2,kk_+1,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_-1,kk_+1,ist_0pp)) + eps);
        Real wmm = std::abs(fsten(ii_,jj_-2,kk_  ,ist_0pp)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_,jj_-1,kk_  ,ist_0pp)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_,jj_-2,kk_+1,ist_0pp)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_,jj_-1,kk_+1,ist_0pp)) * (1._rt + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_pmp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = 1._rt;
        r += std::abs(fsten(ii_  ,jj_-1,kk_+1,ist_p00)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-2,kk_+1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_-1,kk_+1,ist_0p0)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_+1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_-1,kk_  ,ist_00p)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_  ,ist_ppp)) + eps);
        r += std::abs(fsten(ii_  ,jj_-1,kk_+1,ist_pp0)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_p0p)) /
           ( std::abs(fsten(ii_  ,jj_-2,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_-1,kk_  ,ist_0pp)) /
           ( std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_-1,kk_  ,ist_ppp)) + eps);
        r *= std::abs(fsten(ii_  ,jj_-1,kk_  ,ist_ppp)) * fsten(ii_+1,jj_-1,kk_+1,ist_inv);
        return r;
    };

    auto restrict_from_m0p_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_-2,jj_,kk_+1,ist_p00)) / (std::abs(fsten(ii_-2,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_-2,jj_,kk_+1,ist_p0p)) + eps);
        Real w1p = std::abs(fsten(ii_-1,jj_,kk_+1,ist_p00)) / (std::abs(fsten(ii_-1,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_+1,ist_p0p)) + eps);
        Real w2m = std::abs(fsten(ii_-1,jj_,kk_  ,ist_00p)) / (std::abs(fsten(ii_-2,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_  ,ist_p0p)) + eps);
        Real w2p = std::abs(fsten(ii_-1,jj_,kk_+1,ist_00p)) / (std::abs(fsten(ii_-2,jj_,kk_+1,ist_p0p))
                                                              +std::abs(fsten(ii_-1,jj_,kk_+1,ist_p0p)) + eps);
        Real wmm = std::abs(fsten(ii_-2,jj_,kk_  ,ist_p0p)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_-1,jj_,kk_  ,ist_p0p)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_-2,jj_,kk_+1,ist_p0p)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_-1,jj_,kk_+1,ist_p0p)) * (1._rt + w1p + w2p);
        return wpm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_00p_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1 = std::abs(fsten(ii_,jj_,kk_  ,ist_00p));
        Real w2 = std::abs(fsten(ii_,jj_,kk_+1,ist_00p));
        if (w1 == 0._rt and w2 == 0._rt) {
            return 0.5_rt;
        } else {
            return w1 / (w1+w2);
        }
    };

    auto restrict_from_p0p_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_  ,jj_,kk_+1,ist_p00)) / (std::abs(fsten(ii_  ,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_  ,jj_,kk_+1,ist_p0p)) + eps);
        Real w1p = std::abs(fsten(ii_+1,jj_,kk_+1,ist_p00)) / (std::abs(fsten(ii_+1,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_+1,ist_p0p)) + eps);
        Real w2m = std::abs(fsten(ii_+1,jj_,kk_  ,ist_00p)) / (std::abs(fsten(ii_  ,jj_,kk_  ,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_  ,ist_p0p)) + eps);
        Real w2p = std::abs(fsten(ii_+1,jj_,kk_+1,ist_00p)) / (std::abs(fsten(ii_  ,jj_,kk_+1,ist_p0p))
                                                              +std::abs(fsten(ii_+1,jj_,kk_+1,ist_p0p)) + eps);
        Real wmm = std::abs(fsten(ii_  ,jj_,kk_  ,ist_p0p)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_+1,jj_,kk_  ,ist_p0p)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_  ,jj_,kk_+1,ist_p0p)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_+1,jj_,kk_+1,ist_p0p)) * (1._rt + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_mpp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = 1._rt;
        r += std::abs(fsten(ii_-1,jj_+1,kk_+1,ist_p00)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_+1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_  ,kk_+1,ist_0p0)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-2,jj_  ,kk_+1,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_+1,kk_  ,ist_00p)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-2,jj_+1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_  ,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_  ,kk_+1,ist_pp0)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_+1,kk_  ,ist_p0p)) /
           ( std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_+1,kk_  ,ist_ppp)) + eps);
        r += std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_0pp)) /
           ( std::abs(fsten(ii_-2,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp)) + eps);
        r *= std::abs(fsten(ii_-1,jj_  ,kk_  ,ist_ppp)) * fsten(ii_-1,jj_+1,kk_+1,ist_inv);
        return r;
    };

    auto restrict_from_0pp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real w1m = std::abs(fsten(ii_,jj_  ,kk_+1,ist_0p0)) / (std::abs(fsten(ii_,jj_  ,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_  ,kk_+1,ist_0pp)) + eps);
        Real w1p = std::abs(fsten(ii_,jj_+1,kk_+1,ist_0p0)) / (std::abs(fsten(ii_,jj_+1,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_+1,ist_0pp)) + eps);
        Real w2m = std::abs(fsten(ii_,jj_+1,kk_  ,ist_00p)) / (std::abs(fsten(ii_,jj_  ,kk_  ,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_  ,ist_0pp)) + eps);
        Real w2p = std::abs(fsten(ii_,jj_+1,kk_+1,ist_00p)) / (std::abs(fsten(ii_,jj_  ,kk_+1,ist_0pp))
                                                              +std::abs(fsten(ii_,jj_+1,kk_+1,ist_0pp)) + eps);
        Real wmm = std::abs(fsten(ii_,jj_  ,kk_  ,ist_0pp)) * (1._rt + w1m + w2m);
        Real wpm = std::abs(fsten(ii_,jj_+1,kk_  ,ist_0pp)) * (1._rt + w1p + w2m);
        Real wmp = std::abs(fsten(ii_,jj_  ,kk_+1,ist_0pp)) * (1._rt + w1m + w2p);
        Real wpp = std::abs(fsten(ii_,jj_+1,kk_+1,ist_0pp)) * (1._rt + w1p + w2p);
        return wmm / (wmm+wpm+wmp+wpp+eps);
    };

    auto restrict_from_ppp_to = [&fsten] (int ii_, int jj_, int kk_) -> Real {
        Real r = 1._rt;
        r += std::abs(fsten(ii_  ,jj_+1,kk_+1,ist_p00)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_+1,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_  ,kk_+1,ist_0p0)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_+1,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_+1,kk_  ,ist_00p)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_+1,kk_  ,ist_ppp)) + eps);
        r += std::abs(fsten(ii_  ,jj_  ,kk_+1,ist_pp0)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_  ,kk_+1,ist_ppp)) + eps);
        r += std::abs(fsten(ii_  ,jj_+1,kk_  ,ist_p0p)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_  ,jj_+1,kk_  ,ist_ppp)) + eps);
        r += std::abs(fsten(ii_+1,jj_  ,kk_  ,ist_0pp)) /
           ( std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp))
           + std::abs(fsten(ii_+1,jj_  ,kk_  ,ist_ppp)) + eps);
        r *= std::abs(fsten(ii_  ,jj_  ,kk_  ,ist_ppp)) * fsten(ii_+1,jj_+1,kk_+1,ist_inv);
        return r;
    };

    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    Array3D<Real,-1,1,-1,1,-1,1> p;
    Array3D<Real,-1,1,-1,1,-1,1> ap;
    Real cs1, cs2, cs3, cs4;

    // csten(i,j,k,ist_p00)
    int iii = ii;
    int jjj = jj;
    int kkk = kk;;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj-1,kkk-1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj-1,kkk-1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj  ,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj  ,kkk-1);
    p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj+1,kkk-1);
    p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj+1,kkk-1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj-1,kkk  );
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj-1,kkk  );
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj  ,kkk  );
    p( 0, 0, 0) = 1._rt;
    p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj+1,kkk  );
    p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj+1,kkk  );
    p(-1,-1,+1) = interp_from_ppm_to(iii+1,jjj-1,kkk+1);
    p( 0,-1,+1) = interp_from_0pm_to(iii+2,jjj-1,kkk+1);
    p(-1, 0,+1) = interp_from_p0m_to(iii+1,jjj  ,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii+2,jjj  ,kkk+1);
    p(-1,+1,+1) = interp_from_pmm_to(iii+1,jjj+1,kkk+1);
    p( 0,+1,+1) = interp_from_0mm_to(iii+2,jjj+1,kkk+1);
    ap(0,-1,-1) =
                     Ap00(iii,jjj-1,kkk-1) * p(-1,-1,-1)
      +              App0(iii,jjj-1,kkk-1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj-1,kkk-1) * p(-1,-1, 0)
      +              Appp(iii,jjj-1,kkk-1) * p(-1, 0, 0);
    ap(1,-1,-1) =
                     A000(iii+1,jjj-1,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj-1,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj-1,kkk-1) * p(-1, 0,-1)
      +              App0(iii+1,jjj-1,kkk-1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj-1,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj-1,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj-1,kkk-1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj-1,kkk-1) * p( 0, 0, 0);
    ap(0,0,-1) =
                     Apm0(iii,jjj,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii,jjj,kkk-1) * p(-1, 0,-1)
      +              App0(iii,jjj,kkk-1) * p(-1,+1,-1)
      +              Apmp(iii,jjj,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii,jjj,kkk-1) * p(-1, 0, 0)
      +              Appp(iii,jjj,kkk-1) * p(-1,+1, 0);
    ap(1,0,-1) =
                     A0m0(iii+1,jjj,kkk-1) * p(-1,-1,-1)
      +              Apm0(iii+1,jjj,kkk-1) * p( 0,-1,-1)
      +              A000(iii+1,jjj,kkk-1) * p(-1, 0,-1)
      +              Ap00(iii+1,jjj,kkk-1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj,kkk-1) * p(-1,+1,-1)
      +              App0(iii+1,jjj,kkk-1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk-1) * p(-1,-1, 0)
      +              Apmp(iii+1,jjj,kkk-1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj,kkk-1) * p(-1, 0, 0)
      +              Ap0p(iii+1,jjj,kkk-1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj,kkk-1) * p(-1,+1, 0)
      +              Appp(iii+1,jjj,kkk-1) * p( 0,+1, 0);
    ap(0,1,-1) =
                     Apm0(iii,jjj+1,kkk-1) * p(-1, 0,-1)
      +              Ap00(iii,jjj+1,kkk-1) * p(-1,+1,-1)
      +              Apmp(iii,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Ap0p(iii,jjj+1,kkk-1) * p(-1,+1, 0);
    ap(1,1,-1) =
                     A0m0(iii+1,jjj+1,kkk-1) * p(-1, 0,-1)
      +              Apm0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A000(iii+1,jjj+1,kkk-1) * p(-1,+1,-1)
      +              Ap00(iii+1,jjj+1,kkk-1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Apmp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj+1,kkk-1) * p(-1,+1, 0)
      +              Ap0p(iii+1,jjj+1,kkk-1) * p( 0,+1, 0);
    ap(0,-1,0) =
                     Ap0m(iii,jjj-1,kkk) * p(-1,-1,-1)
      +              Appm(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              Ap00(iii,jjj-1,kkk) * p(-1,-1, 0)
      +              App0(iii,jjj-1,kkk) * p(-1, 0, 0)
      +              Ap0p(iii,jjj-1,kkk) * p(-1,-1,+1)
      +              Appp(iii,jjj-1,kkk) * p(-1, 0,+1);
    ap(1,-1,0) =
                     A00m(iii+1,jjj-1,kkk) * p(-1,-1,-1)
      +              Ap0m(iii+1,jjj-1,kkk) * p( 0,-1,-1)
      +              A0pm(iii+1,jjj-1,kkk) * p(-1, 0,-1)
      +              Appm(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A000(iii+1,jjj-1,kkk) * p(-1,-1, 0)
      +              Ap00(iii+1,jjj-1,kkk) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj-1,kkk) * p(-1, 0, 0)
      +              App0(iii+1,jjj-1,kkk) * p( 0, 0, 0)
      +              A00p(iii+1,jjj-1,kkk) * p(-1,-1,+1)
      +              Ap0p(iii+1,jjj-1,kkk) * p( 0,-1,+1)
      +              A0pp(iii+1,jjj-1,kkk) * p(-1, 0,+1)
      +              Appp(iii+1,jjj-1,kkk) * p( 0, 0,+1);
    ap(0,0,0) =
                     Apmm(iii,jjj,kkk) * p(-1,-1,-1)
      +              Ap0m(iii,jjj,kkk) * p(-1, 0,-1)
      +              Appm(iii,jjj,kkk) * p(-1,+1,-1)
      +              Apm0(iii,jjj,kkk) * p(-1,-1, 0)
      +              Ap00(iii,jjj,kkk) * p(-1, 0, 0)
      +              App0(iii,jjj,kkk) * p(-1,+1, 0)
      +              Apmp(iii,jjj,kkk) * p(-1,-1,+1)
      +              Ap0p(iii,jjj,kkk) * p(-1, 0,+1)
      +              Appp(iii,jjj,kkk) * p(-1,+1,+1);
    ap(1,0,0) =
                     A0mm(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Apmm(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A00m(iii+1,jjj,kkk) * p(-1, 0,-1)
      +              Ap0m(iii+1,jjj,kkk) * p( 0, 0,-1)
      +              A0pm(iii+1,jjj,kkk) * p(-1,+1,-1)
      +              Appm(iii+1,jjj,kkk) * p( 0,+1,-1)
      +              A0m0(iii+1,jjj,kkk) * p(-1,-1, 0)
      +              Apm0(iii+1,jjj,kkk) * p( 0,-1, 0)
      +              A000(iii+1,jjj,kkk) * p(-1, 0, 0)
      +              Ap00(iii+1,jjj,kkk) * p( 0, 0, 0)
      +              A0p0(iii+1,jjj,kkk) * p(-1,+1, 0)
      +              App0(iii+1,jjj,kkk) * p( 0,+1, 0)
      +              A0mp(iii+1,jjj,kkk) * p(-1,-1,+1)
      +              Apmp(iii+1,jjj,kkk) * p( 0,-1,+1)
      +              A00p(iii+1,jjj,kkk) * p(-1, 0,+1)
      +              Ap0p(iii+1,jjj,kkk) * p( 0, 0,+1)
      +              A0pp(iii+1,jjj,kkk) * p(-1,+1,+1)
      +              Appp(iii+1,jjj,kkk) * p( 0,+1,+1);
    ap(0,1,0) =
                     Apmm(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              Ap0m(iii,jjj+1,kkk) * p(-1,+1,-1)
      +              Apm0(iii,jjj+1,kkk) * p(-1, 0, 0)
      +              Ap00(iii,jjj+1,kkk) * p(-1,+1, 0)
      +              Apmp(iii,jjj+1,kkk) * p(-1, 0,+1)
      +              Ap0p(iii,jjj+1,kkk) * p(-1,+1,+1);
    ap(1,1,0) =
                     A0mm(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Apmm(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A00m(iii+1,jjj+1,kkk) * p(-1,+1,-1)
      +              Ap0m(iii+1,jjj+1,kkk) * p( 0,+1,-1)
      +              A0m0(iii+1,jjj+1,kkk) * p(-1, 0, 0)
      +              Apm0(iii+1,jjj+1,kkk) * p( 0, 0, 0)
      +              A000(iii+1,jjj+1,kkk) * p(-1,+1, 0)
      +              Ap00(iii+1,jjj+1,kkk) * p( 0,+1, 0)
      +              A0mp(iii+1,jjj+1,kkk) * p(-1, 0,+1)
      +              Apmp(iii+1,jjj+1,kkk) * p( 0, 0,+1)
      +              A00p(iii+1,jjj+1,kkk) * p(-1,+1,+1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,+1,+1);
    ap(0,-1,1) =
                     Ap0m(iii,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Appm(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Ap00(iii,jjj-1,kkk+1) * p(-1,-1,+1)
      +              App0(iii,jjj-1,kkk+1) * p(-1, 0,+1);
    ap(1,-1,1) =
                     A00m(iii+1,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii+1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii+1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Appm(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A000(iii+1,jjj-1,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii+1,jjj-1,kkk+1) * p(-1, 0,+1)
      +              App0(iii+1,jjj-1,kkk+1) * p( 0, 0,+1);
    ap(0,0,1) =
                     Apmm(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii,jjj,kkk+1) * p(-1, 0, 0)
      +              Appm(iii,jjj,kkk+1) * p(-1,+1, 0)
      +              Apm0(iii,jjj,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii,jjj,kkk+1) * p(-1, 0,+1)
      +              App0(iii,jjj,kkk+1) * p(-1,+1,+1);
    ap(1,0,1) =
                     A0mm(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Apmm(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A00m(iii+1,jjj,kkk+1) * p(-1, 0, 0)
      +              Ap0m(iii+1,jjj,kkk+1) * p( 0, 0, 0)
      +              A0pm(iii+1,jjj,kkk+1) * p(-1,+1, 0)
      +              Appm(iii+1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0m0(iii+1,jjj,kkk+1) * p(-1,-1,+1)
      +              Apm0(iii+1,jjj,kkk+1) * p( 0,-1,+1)
      +              A000(iii+1,jjj,kkk+1) * p(-1, 0,+1)
      +              Ap00(iii+1,jjj,kkk+1) * p( 0, 0,+1)
      +              A0p0(iii+1,jjj,kkk+1) * p(-1,+1,+1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,+1,+1);
    ap(0,1,1) =
                     Apmm(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Ap0m(iii,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Apm0(iii,jjj+1,kkk+1) * p(-1, 0,+1)
      +              Ap00(iii,jjj+1,kkk+1) * p(-1,+1,+1);
    ap(1,1,1) =
                     A0mm(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Apmm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A00m(iii+1,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Ap0m(iii+1,jjj+1,kkk+1) * p( 0,+1, 0)
      +              A0m0(iii+1,jjj+1,kkk+1) * p(-1, 0,+1)
      +              Apm0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1)
      +              A000(iii+1,jjj+1,kkk+1) * p(-1,+1,+1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,+1,+1);
    csten(i,j,k,ist_p00) = 0.125_rt *
      ( restrict_from_0mm_to(iii,jjj,kkk) * ap( 0,-1,-1)
      + restrict_from_pmm_to(iii,jjj,kkk) * ap(+1,-1,-1)
      + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // csten(i,j,k,ist_0p0)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj+1,kkk-1);
    p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj+1,kkk-1);
    p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj+1,kkk-1);
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj+2,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj+2,kkk-1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj+2,kkk-1);
    p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj+1,kkk  );
    p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj+1,kkk  );
    p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj+1,kkk  );
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj+2,kkk  );
    p( 0, 0, 0) = 1._rt;
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj+2,kkk  );
    p(-1,-1,+1) = interp_from_ppm_to(iii-1,jjj+1,kkk+1);
    p( 0,-1,+1) = interp_from_0pm_to(iii  ,jjj+1,kkk+1);
    p(+1,-1,+1) = interp_from_mpm_to(iii+1,jjj+1,kkk+1);
    p(-1, 0,+1) = interp_from_p0m_to(iii-1,jjj+2,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii  ,jjj+2,kkk+1);
    p(+1, 0,+1) = interp_from_m0m_to(iii+1,jjj+2,kkk+1);
    ap(-1,0,-1) =
                     A0p0(iii-1,jjj,kkk-1) * p(-1,-1,-1)
      +              App0(iii-1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj,kkk-1) * p(-1,-1, 0)
      +              Appp(iii-1,jjj,kkk-1) * p( 0,-1, 0);
    ap(0,0,-1) =
                     Amp0(iii,jjj,kkk-1) * p(-1,-1,-1)
      +              A0p0(iii,jjj,kkk-1) * p( 0,-1,-1)
      +              App0(iii,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk-1) * p(-1,-1, 0)
      +              A0pp(iii,jjj,kkk-1) * p( 0,-1, 0)
      +              Appp(iii,jjj,kkk-1) * p(+1,-1, 0);
    ap(1,0,-1) =
                     Amp0(iii+1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj,kkk-1) * p(+1,-1, 0);
    ap(-1,1,-1) =
                     A000(iii-1,jjj+1,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii-1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj+1,kkk-1) * p(-1, 0,-1)
      +              App0(iii-1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A00p(iii-1,jjj+1,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii-1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Appp(iii-1,jjj+1,kkk-1) * p( 0, 0, 0);
    ap(0,1,-1) =
                     Am00(iii,jjj+1,kkk-1) * p(-1,-1,-1)
      +              A000(iii,jjj+1,kkk-1) * p( 0,-1,-1)
      +              Ap00(iii,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk-1) * p(-1, 0,-1)
      +              A0p0(iii,jjj+1,kkk-1) * p( 0, 0,-1)
      +              App0(iii,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk-1) * p(-1,-1, 0)
      +              A00p(iii,jjj+1,kkk-1) * p( 0,-1, 0)
      +              Ap0p(iii,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk-1) * p(-1, 0, 0)
      +              A0pp(iii,jjj+1,kkk-1) * p( 0, 0, 0)
      +              Appp(iii,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(1,1,-1) =
                     Am00(iii+1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A000(iii+1,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(-1,0,0) =
                     A0pm(iii-1,jjj,kkk) * p(-1,-1,-1)
      +              Appm(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj,kkk) * p(-1,-1, 0)
      +              App0(iii-1,jjj,kkk) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj,kkk) * p(-1,-1,+1)
      +              Appp(iii-1,jjj,kkk) * p( 0,-1,+1);
    ap(0,0,0) =
                     Ampm(iii,jjj,kkk) * p(-1,-1,-1)
      +              A0pm(iii,jjj,kkk) * p( 0,-1,-1)
      +              Appm(iii,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii,jjj,kkk) * p(-1,-1, 0)
      +              A0p0(iii,jjj,kkk) * p( 0,-1, 0)
      +              App0(iii,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii,jjj,kkk) * p(-1,-1,+1)
      +              A0pp(iii,jjj,kkk) * p( 0,-1,+1)
      +              Appp(iii,jjj,kkk) * p(+1,-1,+1);
    ap(1,0,0) =
                     Ampm(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0pm(iii+1,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj,kkk) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj,kkk) * p( 0,-1,+1)
      +              A0pp(iii+1,jjj,kkk) * p(+1,-1,+1);
    ap(-1,1,0) =
                     A00m(iii-1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0m(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pm(iii-1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appm(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A000(iii-1,jjj+1,kkk) * p(-1,-1, 0)
      +              Ap00(iii-1,jjj+1,kkk) * p( 0,-1, 0)
      +              A0p0(iii-1,jjj+1,kkk) * p(-1, 0, 0)
      +              App0(iii-1,jjj+1,kkk) * p( 0, 0, 0)
      +              A00p(iii-1,jjj+1,kkk) * p(-1,-1,+1)
      +              Ap0p(iii-1,jjj+1,kkk) * p( 0,-1,+1)
      +              A0pp(iii-1,jjj+1,kkk) * p(-1, 0,+1)
      +              Appp(iii-1,jjj+1,kkk) * p( 0, 0,+1);
    ap(0,1,0) =
                     Am0m(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              A00m(iii,jjj+1,kkk) * p( 0,-1,-1)
      +              Ap0m(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              A0pm(iii,jjj+1,kkk) * p( 0, 0,-1)
      +              Appm(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk) * p(-1,-1, 0)
      +              A000(iii,jjj+1,kkk) * p( 0,-1, 0)
      +              Ap00(iii,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii,jjj+1,kkk) * p(-1, 0, 0)
      +              A0p0(iii,jjj+1,kkk) * p( 0, 0, 0)
      +              App0(iii,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk) * p(-1,-1,+1)
      +              A00p(iii,jjj+1,kkk) * p( 0,-1,+1)
      +              Ap0p(iii,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii,jjj+1,kkk) * p(-1, 0,+1)
      +              A0pp(iii,jjj+1,kkk) * p( 0, 0,+1)
      +              Appp(iii,jjj+1,kkk) * p(+1, 0,+1);
    ap(1,1,0) =
                     Am0m(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00m(iii+1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pm(iii+1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii+1,jjj+1,kkk) * p( 0,-1, 0)
      +              A000(iii+1,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii+1,jjj+1,kkk) * p( 0, 0, 0)
      +              A0p0(iii+1,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii+1,jjj+1,kkk) * p( 0,-1,+1)
      +              A00p(iii+1,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii+1,jjj+1,kkk) * p( 0, 0,+1)
      +              A0pp(iii+1,jjj+1,kkk) * p(+1, 0,+1);
    ap(-1,0,1) =
                     A0pm(iii-1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appm(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0p0(iii-1,jjj,kkk+1) * p(-1,-1,+1)
      +              App0(iii-1,jjj,kkk+1) * p( 0,-1,+1);
    ap(0,0,1) =
                     Ampm(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              A0pm(iii,jjj,kkk+1) * p( 0,-1, 0)
      +              Appm(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii,jjj,kkk+1) * p(-1,-1,+1)
      +              A0p0(iii,jjj,kkk+1) * p( 0,-1,+1)
      +              App0(iii,jjj,kkk+1) * p(+1,-1,+1);
    ap(1,0,1) =
                     Ampm(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii+1,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii+1,jjj,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii+1,jjj,kkk+1) * p(+1,-1,+1);
    ap(-1,1,1) =
                     A00m(iii-1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii-1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appm(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A000(iii-1,jjj+1,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(-1, 0,+1)
      +              App0(iii-1,jjj+1,kkk+1) * p( 0, 0,+1);
    ap(0,1,1) =
                     Am0m(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              A00m(iii,jjj+1,kkk+1) * p( 0,-1, 0)
      +              Ap0m(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              A0pm(iii,jjj+1,kkk+1) * p( 0, 0, 0)
      +              Appm(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii,jjj+1,kkk+1) * p(-1,-1,+1)
      +              A000(iii,jjj+1,kkk+1) * p( 0,-1,+1)
      +              Ap00(iii,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii,jjj+1,kkk+1) * p(-1, 0,+1)
      +              A0p0(iii,jjj+1,kkk+1) * p( 0, 0,+1)
      +              App0(iii,jjj+1,kkk+1) * p(+1, 0,+1);
    ap(1,1,1) =
                     Am0m(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00m(iii+1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pm(iii+1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii+1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A000(iii+1,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(+1, 0,+1);
    csten(i,j,k,ist_0p0) = 0.125_rt *
      ( restrict_from_m0m_to(iii,jjj,kkk) * ap(-1, 0,-1)
      + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1)
      + restrict_from_mpm_to(iii,jjj,kkk) * ap(-1,+1,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // csten(i,j,k,ist_00p)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj-1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj-1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj-1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj  ,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj  ,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj  ,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii-1,jjj+1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii  ,jjj+1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii+1,jjj+1,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj-1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj-1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj-1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj  ,kkk+2);
    p( 0, 0, 0) = 1._rt;
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj  ,kkk+2);
    p(-1,+1, 0) = interp_from_pm0_to(iii-1,jjj+1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii  ,jjj+1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii+1,jjj+1,kkk+2);
    ap(-1,-1,0) =
                     A00p(iii-1,jjj-1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii-1,jjj-1,kkk) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj-1,kkk) * p(-1, 0,-1)
      +              Appp(iii-1,jjj-1,kkk) * p( 0, 0,-1);
    ap(0,-1,0) =
                     Am0p(iii,jjj-1,kkk) * p(-1,-1,-1)
      +              A00p(iii,jjj-1,kkk) * p( 0,-1,-1)
      +              Ap0p(iii,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              A0pp(iii,jjj-1,kkk) * p( 0, 0,-1)
      +              Appp(iii,jjj-1,kkk) * p(+1, 0,-1);
    ap(1,-1,0) =
                     Am0p(iii+1,jjj-1,kkk) * p( 0,-1,-1)
      +              A00p(iii+1,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0pp(iii+1,jjj-1,kkk) * p(+1, 0,-1);
    ap(-1,0,0) =
                     A0mp(iii-1,jjj,kkk) * p(-1,-1,-1)
      +              Apmp(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A00p(iii-1,jjj,kkk) * p(-1, 0,-1)
      +              Ap0p(iii-1,jjj,kkk) * p( 0, 0,-1)
      +              A0pp(iii-1,jjj,kkk) * p(-1,+1,-1)
      +              Appp(iii-1,jjj,kkk) * p( 0,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(-1,-1,-1)
      +              A0mp(iii,jjj,kkk) * p( 0,-1,-1)
      +              Apmp(iii,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii,jjj,kkk) * p(-1, 0,-1)
      +              A00p(iii,jjj,kkk) * p( 0, 0,-1)
      +              Ap0p(iii,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii,jjj,kkk) * p(-1,+1,-1)
      +              A0pp(iii,jjj,kkk) * p( 0,+1,-1)
      +              Appp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(1,0,0) =
                     Ammp(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0mp(iii+1,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii+1,jjj,kkk) * p( 0, 0,-1)
      +              A00p(iii+1,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii+1,jjj,kkk) * p( 0,+1,-1)
      +              A0pp(iii+1,jjj,kkk) * p(+1,+1,-1);
    ap(-1,1,0) =
                     A0mp(iii-1,jjj+1,kkk) * p(-1, 0,-1)
      +              Apmp(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A00p(iii-1,jjj+1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii-1,jjj+1,kkk) * p( 0,+1,-1);
    ap(0,1,0) =
                     Ammp(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              A0mp(iii,jjj+1,kkk) * p( 0, 0,-1)
      +              Apmp(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk) * p(-1,+1,-1)
      +              A00p(iii,jjj+1,kkk) * p( 0,+1,-1)
      +              Ap0p(iii,jjj+1,kkk) * p(+1,+1,-1);
    ap(1,1,0) =
                     Ammp(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0mp(iii+1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj+1,kkk) * p( 0,+1,-1)
      +              A00p(iii+1,jjj+1,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     A000(iii-1,jjj-1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii-1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              App0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii-1,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii-1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0);
    ap(0,-1,1) =
                     Am00(iii,jjj-1,kkk+1) * p(-1,-1,-1)
      +              A000(iii,jjj-1,kkk+1) * p( 0,-1,-1)
      +              Ap00(iii,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              A0p0(iii,jjj-1,kkk+1) * p( 0, 0,-1)
      +              App0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk+1) * p(-1,-1, 0)
      +              A00p(iii,jjj-1,kkk+1) * p( 0,-1, 0)
      +              Ap0p(iii,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              A0pp(iii,jjj-1,kkk+1) * p( 0, 0, 0)
      +              Appp(iii,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(1,-1,1) =
                     Am00(iii+1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A000(iii+1,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(-1,0,1) =
                     A0m0(iii-1,jjj,kkk+1) * p(-1,-1,-1)
      +              Apm0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A000(iii-1,jjj,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii-1,jjj,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj,kkk+1) * p(-1,+1,-1)
      +              App0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj,kkk+1) * p(-1,-1, 0)
      +              Apmp(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii-1,jjj,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj,kkk+1) * p(-1,+1, 0)
      +              Appp(iii-1,jjj,kkk+1) * p( 0,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              A0m0(iii,jjj,kkk+1) * p( 0,-1,-1)
      +              Apm0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii,jjj,kkk+1) * p(-1, 0,-1)
      +              A000(iii,jjj,kkk+1) * p( 0, 0,-1)
      +              Ap00(iii,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              A0p0(iii,jjj,kkk+1) * p( 0,+1,-1)
      +              App0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              A0mp(iii,jjj,kkk+1) * p( 0,-1, 0)
      +              Apmp(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii,jjj,kkk+1) * p(-1, 0, 0)
      +              A00p(iii,jjj,kkk+1) * p( 0, 0, 0)
      +              Ap0p(iii,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii,jjj,kkk+1) * p(-1,+1, 0)
      +              A0pp(iii,jjj,kkk+1) * p( 0,+1, 0)
      +              Appp(iii,jjj,kkk+1) * p(+1,+1, 0);
    ap(1,0,1) =
                     Amm0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0m0(iii+1,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii+1,jjj,kkk+1) * p( 0, 0,-1)
      +              A000(iii+1,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0p0(iii+1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0mp(iii+1,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii+1,jjj,kkk+1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii+1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0pp(iii+1,jjj,kkk+1) * p(+1,+1, 0);
    ap(-1,1,1) =
                     A0m0(iii-1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A000(iii-1,jjj+1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii-1,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii-1,jjj+1,kkk+1) * p( 0,+1, 0);
    ap(0,1,1) =
                     Amm0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              A0m0(iii,jjj+1,kkk+1) * p( 0, 0,-1)
      +              Apm0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk+1) * p(-1,+1,-1)
      +              A000(iii,jjj+1,kkk+1) * p( 0,+1,-1)
      +              Ap00(iii,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              A0mp(iii,jjj+1,kkk+1) * p( 0, 0, 0)
      +              Apmp(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk+1) * p(-1,+1, 0)
      +              A00p(iii,jjj+1,kkk+1) * p( 0,+1, 0)
      +              Ap0p(iii,jjj+1,kkk+1) * p(+1,+1, 0);
    ap(1,1,1) =
                     Amm0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii+1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii+1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A000(iii+1,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii+1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii+1,jjj+1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii+1,jjj+1,kkk+1) * p(+1,+1, 0);
    csten(i,j,k,ist_00p) = 0.125_rt *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // csten(i,j,k,ist_pp0)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj+1,kkk-1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj+1,kkk-1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj+2,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj+2,kkk-1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj+1,kkk  );
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj+1,kkk  );
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj+2,kkk  );
    p( 0, 0, 0) = 1._rt;
    p(-1,-1,+1) = interp_from_ppm_to(iii+1,jjj+1,kkk+1);
    p( 0,-1,+1) = interp_from_0pm_to(iii+2,jjj+1,kkk+1);
    p(-1, 0,+1) = interp_from_p0m_to(iii+1,jjj+2,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii+2,jjj+2,kkk+1);
    ap(0,0,-1) =
                     App0(iii,jjj,kkk-1) * p(-1,-1,-1)
      +              Appp(iii,jjj,kkk-1) * p(-1,-1, 0);
    ap(1,0,-1) =
                     A0p0(iii+1,jjj,kkk-1) * p(-1,-1,-1)
      +              App0(iii+1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj,kkk-1) * p(-1,-1, 0)
      +              Appp(iii+1,jjj,kkk-1) * p( 0,-1, 0);
    ap(0,1,-1) =
                     Ap00(iii,jjj+1,kkk-1) * p(-1,-1,-1)
      +              App0(iii,jjj+1,kkk-1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj+1,kkk-1) * p(-1,-1, 0)
      +              Appp(iii,jjj+1,kkk-1) * p(-1, 0, 0);
    ap(1,1,-1) =
                     A000(iii+1,jjj+1,kkk-1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj+1,kkk-1) * p(-1, 0,-1)
      +              App0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj+1,kkk-1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj+1,kkk-1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0);
    ap(0,0,0) =
                     Appm(iii,jjj,kkk) * p(-1,-1,-1)
      +              App0(iii,jjj,kkk) * p(-1,-1, 0)
      +              Appp(iii,jjj,kkk) * p(-1,-1,+1);
    ap(1,0,0) =
                     A0pm(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Appm(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj,kkk) * p(-1,-1, 0)
      +              App0(iii+1,jjj,kkk) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj,kkk) * p(-1,-1,+1)
      +              Appp(iii+1,jjj,kkk) * p( 0,-1,+1);
    ap(0,1,0) =
                     Ap0m(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              Appm(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              Ap00(iii,jjj+1,kkk) * p(-1,-1, 0)
      +              App0(iii,jjj+1,kkk) * p(-1, 0, 0)
      +              Ap0p(iii,jjj+1,kkk) * p(-1,-1,+1)
      +              Appp(iii,jjj+1,kkk) * p(-1, 0,+1);
    ap(1,1,0) =
                     A00m(iii+1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0m(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pm(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appm(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A000(iii+1,jjj+1,kkk) * p(-1,-1, 0)
      +              Ap00(iii+1,jjj+1,kkk) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj+1,kkk) * p(-1, 0, 0)
      +              App0(iii+1,jjj+1,kkk) * p( 0, 0, 0)
      +              A00p(iii+1,jjj+1,kkk) * p(-1,-1,+1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,-1,+1)
      +              A0pp(iii+1,jjj+1,kkk) * p(-1, 0,+1)
      +              Appp(iii+1,jjj+1,kkk) * p( 0, 0,+1);
    ap(0,0,1) =
                     Appm(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              App0(iii,jjj,kkk+1) * p(-1,-1,+1);
    ap(1,0,1) =
                     A0pm(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appm(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0p0(iii+1,jjj,kkk+1) * p(-1,-1,+1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,-1,+1);
    ap(0,1,1) =
                     Ap0m(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Appm(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Ap00(iii,jjj+1,kkk+1) * p(-1,-1,+1)
      +              App0(iii,jjj+1,kkk+1) * p(-1, 0,+1);
    ap(1,1,1) =
                     A00m(iii+1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0m(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A000(iii+1,jjj+1,kkk+1) * p(-1,-1,+1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(-1, 0,+1)
      +              App0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1);
    cs1 = 0.125_rt *
      ( restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i+1,j,k,ist_mp0)
    iii = ii+2;
    jjj = jj;
    kkk = kk;
    p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj+1,kkk-1);
    p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj+1,kkk-1);
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj+2,kkk-1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj+2,kkk-1);
    p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj+1,kkk  );
    p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj+1,kkk  );
    p( 0, 0, 0) = 1._rt;
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj+2,kkk  );
    p( 0,-1,+1) = interp_from_0pm_to(iii-2,jjj+1,kkk+1);
    p(+1,-1,+1) = interp_from_mpm_to(iii-1,jjj+1,kkk+1);
    p( 0, 0,+1) = interp_from_00m_to(iii-2,jjj+2,kkk+1);
    p(+1, 0,+1) = interp_from_m0m_to(iii-1,jjj+2,kkk+1);
    ap(-1,0,-1) =
                     Amp0(iii-1,jjj,kkk-1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj,kkk-1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj,kkk-1) * p(+1,-1, 0);
    ap(0,0,-1) =
                     Amp0(iii,jjj,kkk-1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk-1) * p(+1,-1, 0);
    ap(-1,1,-1) =
                     Am00(iii-1,jjj+1,kkk-1) * p( 0,-1,-1)
      +              A000(iii-1,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj+1,kkk-1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj+1,kkk-1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj+1,kkk-1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(0,1,-1) =
                     Am00(iii,jjj+1,kkk-1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk-1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk-1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk-1) * p(+1, 0, 0);
    ap(-1,0,0) =
                     Ampm(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0pm(iii-1,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj,kkk) * p( 0,-1, 0)
      +              A0p0(iii-1,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj,kkk) * p( 0,-1,+1)
      +              A0pp(iii-1,jjj,kkk) * p(+1,-1,+1);
    ap(0,0,0) =
                     Ampm(iii,jjj,kkk) * p(+1,-1,-1)
      +              Amp0(iii,jjj,kkk) * p(+1,-1, 0)
      +              Ampp(iii,jjj,kkk) * p(+1,-1,+1);
    ap(-1,1,0) =
                     Am0m(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00m(iii-1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pm(iii-1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii-1,jjj+1,kkk) * p( 0,-1, 0)
      +              A000(iii-1,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii-1,jjj+1,kkk) * p( 0, 0, 0)
      +              A0p0(iii-1,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii-1,jjj+1,kkk) * p( 0,-1,+1)
      +              A00p(iii-1,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii-1,jjj+1,kkk) * p( 0, 0,+1)
      +              A0pp(iii-1,jjj+1,kkk) * p(+1, 0,+1);
    ap(0,1,0) =
                     Am0m(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampm(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk) * p(+1,-1, 0)
      +              Amp0(iii,jjj+1,kkk) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk) * p(+1,-1,+1)
      +              Ampp(iii,jjj+1,kkk) * p(+1, 0,+1);
    ap(-1,0,1) =
                     Ampm(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pm(iii-1,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii-1,jjj,kkk+1) * p( 0,-1,+1)
      +              A0p0(iii-1,jjj,kkk+1) * p(+1,-1,+1);
    ap(0,0,1) =
                     Ampm(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Amp0(iii,jjj,kkk+1) * p(+1,-1,+1);
    ap(-1,1,1) =
                     Am0m(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00m(iii-1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pm(iii-1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii-1,jjj+1,kkk+1) * p( 0,-1,+1)
      +              A000(iii-1,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii-1,jjj+1,kkk+1) * p( 0, 0,+1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(+1, 0,+1);
    ap(0,1,1) =
                     Am0m(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampm(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am00(iii,jjj+1,kkk+1) * p(+1,-1,+1)
      +              Amp0(iii,jjj+1,kkk+1) * p(+1, 0,+1);
    cs2 = 0.125_rt *
      ( restrict_from_m0m_to(iii,jjj,kkk) * ap(-1, 0,-1)
      + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1)
      + restrict_from_mpm_to(iii,jjj,kkk) * ap(-1,+1,-1)
      + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1));

    csten(i,j,k,ist_pp0) = 0.5_rt*(cs1 + cs2);

    // csten(i,j,k,ist_p0p)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj-1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj-1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj  ,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj  ,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj+1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj+1,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj-1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj-1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj  ,kkk+2);
    p( 0, 0, 0) = 1._rt;
    p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj+1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj+1,kkk+2);
    ap(0,-1,0) =
                     Ap0p(iii,jjj-1,kkk) * p(-1,-1,-1)
      +              Appp(iii,jjj-1,kkk) * p(-1, 0,-1);
    ap(1,-1,0) =
                     A00p(iii+1,jjj-1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii+1,jjj-1,kkk) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj-1,kkk) * p(-1, 0,-1)
      +              Appp(iii+1,jjj-1,kkk) * p( 0, 0,-1);
    ap(0,0,0) =
                     Apmp(iii,jjj,kkk) * p(-1,-1,-1)
      +              Ap0p(iii,jjj,kkk) * p(-1, 0,-1)
      +              Appp(iii,jjj,kkk) * p(-1,+1,-1);
    ap(1,0,0) =
                     A0mp(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Apmp(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A00p(iii+1,jjj,kkk) * p(-1, 0,-1)
      +              Ap0p(iii+1,jjj,kkk) * p( 0, 0,-1)
      +              A0pp(iii+1,jjj,kkk) * p(-1,+1,-1)
      +              Appp(iii+1,jjj,kkk) * p( 0,+1,-1);
    ap(0,1,0) =
                     Apmp(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              Ap0p(iii,jjj+1,kkk) * p(-1,+1,-1);
    ap(1,1,0) =
                     A0mp(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Apmp(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A00p(iii+1,jjj+1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,+1,-1);
    ap(0,-1,1) =
                     Ap00(iii,jjj-1,kkk+1) * p(-1,-1,-1)
      +              App0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Appp(iii,jjj-1,kkk+1) * p(-1, 0, 0);
    ap(1,-1,1) =
                     A000(iii+1,jjj-1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              App0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj-1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0);
    ap(0,0,1) =
                     Apm0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii,jjj,kkk+1) * p(-1, 0,-1)
      +              App0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii,jjj,kkk+1) * p(-1, 0, 0)
      +              Appp(iii,jjj,kkk+1) * p(-1,+1, 0);
    ap(1,0,1) =
                     A0m0(iii+1,jjj,kkk+1) * p(-1,-1,-1)
      +              Apm0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A000(iii+1,jjj,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii+1,jjj,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj,kkk+1) * p(-1,+1,-1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Apmp(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii+1,jjj,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj,kkk+1) * p(-1,+1, 0)
      +              Appp(iii+1,jjj,kkk+1) * p( 0,+1, 0);
    ap(0,1,1) =
                     Apm0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii,jjj+1,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii,jjj+1,kkk+1) * p(-1,+1, 0);
    ap(1,1,1) =
                     A0m0(iii+1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A000(iii+1,jjj+1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj+1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii+1,jjj+1,kkk+1) * p( 0,+1, 0);
    cs1 = 0.125_rt *
      ( restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i+1,j,k,ist_m0p)
    iii = ii+2;
    jjj = jj;
    kkk = kk;
    p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj-1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj-1,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj  ,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj  ,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii-2,jjj+1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii-1,jjj+1,kkk+1);
    p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj-1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj-1,kkk+2);
    p( 0, 0, 0) = 1._rt;
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj  ,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii-2,jjj+1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii-1,jjj+1,kkk+2);

    ap(-1,-1,0) =
                     Am0p(iii-1,jjj-1,kkk) * p( 0,-1,-1)
      +              A00p(iii-1,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0pp(iii-1,jjj-1,kkk) * p(+1, 0,-1);
    ap(0,-1,0) =
                     Am0p(iii,jjj-1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj-1,kkk) * p(+1, 0,-1);
    ap(-1,0,0) =
                     Ammp(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0mp(iii-1,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii-1,jjj,kkk) * p( 0, 0,-1)
      +              A00p(iii-1,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii-1,jjj,kkk) * p( 0,+1,-1)
      +              A0pp(iii-1,jjj,kkk) * p(+1,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(+1,-1,-1)
      +              Am0p(iii,jjj,kkk) * p(+1, 0,-1)
      +              Ampp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(-1,1,0) =
                     Ammp(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0mp(iii-1,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj+1,kkk) * p( 0,+1,-1)
      +              A00p(iii-1,jjj+1,kkk) * p(+1,+1,-1);
    ap(0,1,0) =
                     Ammp(iii,jjj+1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     Am00(iii-1,jjj-1,kkk+1) * p( 0,-1,-1)
      +              A000(iii-1,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj-1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(0,-1,1) =
                     Am00(iii,jjj-1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj-1,kkk+1) * p(+1, 0, 0);
    ap(-1,0,1) =
                     Amm0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0m0(iii-1,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii-1,jjj,kkk+1) * p( 0, 0,-1)
      +              A000(iii-1,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0p0(iii-1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0mp(iii-1,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii-1,jjj,kkk+1) * p( 0, 0, 0)
      +              A00p(iii-1,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii-1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0pp(iii-1,jjj,kkk+1) * p(+1,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Am00(iii,jjj,kkk+1) * p(+1, 0,-1)
      +              Amp0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(+1,-1, 0)
      +              Am0p(iii,jjj,kkk+1) * p(+1, 0, 0)
      +              Ampp(iii,jjj,kkk+1) * p(+1,+1, 0);
    ap(-1,1,1) =
                     Amm0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii-1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii-1,jjj+1,kkk+1) * p( 0,+1,-1)
      +              A000(iii-1,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii-1,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii-1,jjj+1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii-1,jjj+1,kkk+1) * p(+1,+1, 0);
    ap(0,1,1) =
                     Amm0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj+1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj+1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj+1,kkk+1) * p(+1,+1, 0);
    cs2 = 0.125_rt *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1));

    csten(i,j,k,ist_p0p) = 0.5_rt*(cs1+cs2);

    // csten(i,j,k,ist_0pp)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj+1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj+1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj+1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj+2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj+2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj+2,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj+1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj+1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj+1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj+2,kkk+2);
    p( 0, 0, 0) = 1._rt;
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj+2,kkk+2);
    ap(-1,0,0) =
                     A0pp(iii-1,jjj,kkk) * p(-1,-1,-1)
      +              Appp(iii-1,jjj,kkk) * p( 0,-1,-1);
    ap(0,0,0) =
                     Ampp(iii,jjj,kkk) * p(-1,-1,-1)
      +              A0pp(iii,jjj,kkk) * p( 0,-1,-1)
      +              Appp(iii,jjj,kkk) * p(+1,-1,-1);
    ap(1,0,0) =
                     Ampp(iii+1,jjj,kkk) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj,kkk) * p(+1,-1,-1);
    ap(-1,1,0) =
                     A00p(iii-1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appp(iii-1,jjj+1,kkk) * p( 0, 0,-1);
    ap(0,1,0) =
                     Am0p(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              A00p(iii,jjj+1,kkk) * p( 0,-1,-1)
      +              Ap0p(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj+1,kkk) * p(-1, 0,-1)
      +              A0pp(iii,jjj+1,kkk) * p( 0, 0,-1)
      +              Appp(iii,jjj+1,kkk) * p(+1, 0,-1);
    ap(1,1,0) =
                     Am0p(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00p(iii+1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pp(iii+1,jjj+1,kkk) * p(+1, 0,-1);
    ap(-1,0,1) =
                     A0p0(iii-1,jjj,kkk+1) * p(-1,-1,-1)
      +              App0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appp(iii-1,jjj,kkk+1) * p( 0,-1, 0);
    ap(0,0,1) =
                     Amp0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              A0p0(iii,jjj,kkk+1) * p( 0,-1,-1)
      +              App0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk+1) * p(-1,-1, 0)
      +              A0pp(iii,jjj,kkk+1) * p( 0,-1, 0)
      +              Appp(iii,jjj,kkk+1) * p(+1,-1, 0);
    ap(1,0,1) =
                     Amp0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii+1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj,kkk+1) * p(+1,-1, 0);
    ap(-1,1,1) =
                     A000(iii-1,jjj+1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              App0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii-1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0);
    ap(0,1,1) =
                     Am00(iii,jjj+1,kkk+1) * p(-1,-1,-1)
      +              A000(iii,jjj+1,kkk+1) * p( 0,-1,-1)
      +              Ap00(iii,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              A0p0(iii,jjj+1,kkk+1) * p( 0, 0,-1)
      +              App0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              A00p(iii,jjj+1,kkk+1) * p( 0,-1, 0)
      +              Ap0p(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk+1) * p(-1, 0, 0)
      +              A0pp(iii,jjj+1,kkk+1) * p( 0, 0, 0)
      +              Appp(iii,jjj+1,kkk+1) * p(+1, 0, 0);
    ap(1,1,1) =
                     Am00(iii+1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A000(iii+1,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii+1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii+1,jjj+1,kkk+1) * p(+1, 0, 0);
    cs1 = 0.125_rt *
      ( restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i,j+1,k,ist_0mp)
    iii = ii;
    jjj = jj+2;
    kkk = kk;
    p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj-2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj-2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj-2,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii-1,jjj-1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii  ,jjj-1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii+1,jjj-1,kkk+1);
    p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj-2,kkk+2);
    p( 0, 0, 0) = 1._rt;
    p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj-2,kkk+2);
    p(-1,+1, 0) = interp_from_pm0_to(iii-1,jjj-1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii  ,jjj-1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii+1,jjj-1,kkk+2);
    ap(-1,-1,0) =
                     A0mp(iii-1,jjj-1,kkk) * p(-1, 0,-1)
      +              Apmp(iii-1,jjj-1,kkk) * p( 0, 0,-1)
      +              A00p(iii-1,jjj-1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii-1,jjj-1,kkk) * p( 0,+1,-1);
    ap(0,-1,0) =
                     Ammp(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              A0mp(iii,jjj-1,kkk) * p( 0, 0,-1)
      +              Apmp(iii,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk) * p(-1,+1,-1)
      +              A00p(iii,jjj-1,kkk) * p( 0,+1,-1)
      +              Ap0p(iii,jjj-1,kkk) * p(+1,+1,-1);
    ap(1,-1,0) =
                     Ammp(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0mp(iii+1,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii+1,jjj-1,kkk) * p( 0,+1,-1)
      +              A00p(iii+1,jjj-1,kkk) * p(+1,+1,-1);
    ap(-1,0,0) =
                     A0mp(iii-1,jjj,kkk) * p(-1,+1,-1)
      +              Apmp(iii-1,jjj,kkk) * p( 0,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(-1,+1,-1)
      +              A0mp(iii,jjj,kkk) * p( 0,+1,-1)
      +              Apmp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(1,0,0) =
                     Ammp(iii+1,jjj,kkk) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     A0m0(iii-1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A000(iii-1,jjj-1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii-1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii-1,jjj-1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii-1,jjj-1,kkk+1) * p( 0,+1, 0);
    ap(0,-1,1) =
                     Amm0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              A0m0(iii,jjj-1,kkk+1) * p( 0, 0,-1)
      +              Apm0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj-1,kkk+1) * p(-1,+1,-1)
      +              A000(iii,jjj-1,kkk+1) * p( 0,+1,-1)
      +              Ap00(iii,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              A0mp(iii,jjj-1,kkk+1) * p( 0, 0, 0)
      +              Apmp(iii,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj-1,kkk+1) * p(-1,+1, 0)
      +              A00p(iii,jjj-1,kkk+1) * p( 0,+1, 0)
      +              Ap0p(iii,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(1,-1,1) =
                     Amm0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii+1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii+1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A000(iii+1,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii+1,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii+1,jjj-1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii+1,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(-1,0,1) =
                     A0m0(iii-1,jjj,kkk+1) * p(-1,+1,-1)
      +              Apm0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj,kkk+1) * p(-1,+1, 0)
      +              Apmp(iii-1,jjj,kkk+1) * p( 0,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              A0m0(iii,jjj,kkk+1) * p( 0,+1,-1)
      +              Apm0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(-1,+1, 0)
      +              A0mp(iii,jjj,kkk+1) * p( 0,+1, 0)
      +              Apmp(iii,jjj,kkk+1) * p(+1,+1, 0);
    ap(1,0,1) =
                     Amm0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0m0(iii+1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii+1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0mp(iii+1,jjj,kkk+1) * p(+1,+1, 0);
    cs2 = 0.125_rt *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1));

    csten(i,j,k,ist_0pp) = 0.5_rt*(cs1+cs2);

    // csten(i,j,k,ist_ppp)
    iii = ii;
    jjj = jj;
    kkk = kk;
    p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj+1,kkk+1);
    p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj+1,kkk+1);
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj+2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj+2,kkk+1);
    p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj+1,kkk+2);
    p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj+1,kkk+2);
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj+2,kkk+2);
    p( 0, 0, 0) = 1._rt;
    ap(0,0,0) =
                     Appp(iii,jjj,kkk) * p(-1,-1,-1);
    ap(1,0,0) =
                     A0pp(iii+1,jjj,kkk) * p(-1,-1,-1)
      +              Appp(iii+1,jjj,kkk) * p( 0,-1,-1);
    ap(0,1,0) =
                     Ap0p(iii,jjj+1,kkk) * p(-1,-1,-1)
      +              Appp(iii,jjj+1,kkk) * p(-1, 0,-1);
    ap(1,1,0) =
                     A00p(iii+1,jjj+1,kkk) * p(-1,-1,-1)
      +              Ap0p(iii+1,jjj+1,kkk) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj+1,kkk) * p(-1, 0,-1)
      +              Appp(iii+1,jjj+1,kkk) * p( 0, 0,-1);
    ap(0,0,1) =
                     App0(iii,jjj,kkk+1) * p(-1,-1,-1)
      +              Appp(iii,jjj,kkk+1) * p(-1,-1, 0);
    ap(1,0,1) =
                     A0p0(iii+1,jjj,kkk+1) * p(-1,-1,-1)
      +              App0(iii+1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0pp(iii+1,jjj,kkk+1) * p(-1,-1, 0)
      +              Appp(iii+1,jjj,kkk+1) * p( 0,-1, 0);
    ap(0,1,1) =
                     Ap00(iii,jjj+1,kkk+1) * p(-1,-1,-1)
      +              App0(iii,jjj+1,kkk+1) * p(-1, 0,-1)
      +              Ap0p(iii,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Appp(iii,jjj+1,kkk+1) * p(-1, 0, 0);
    ap(1,1,1) =
                     A000(iii+1,jjj+1,kkk+1) * p(-1,-1,-1)
      +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii+1,jjj+1,kkk+1) * p(-1, 0,-1)
      +              App0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A00p(iii+1,jjj+1,kkk+1) * p(-1,-1, 0)
      +              Ap0p(iii+1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii+1,jjj+1,kkk+1) * p(-1, 0, 0)
      +              Appp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0);
    cs1 = 0.125_rt *
      ( restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1)
      + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1));

    // alternative: csten(i+1,j,k,ist_mpp)
    iii = ii+2;
    jjj = jj;
    kkk = kk;
    p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj+1,kkk+1);
    p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj+1,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj+2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj+2,kkk+1);
    p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj+1,kkk+2);
    p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj+1,kkk+2);
    p( 0, 0, 0) = 1._rt;
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj+2,kkk+2);
    ap(-1,0,0) =
                     Ampp(iii-1,jjj,kkk) * p( 0,-1,-1)
      +              A0pp(iii-1,jjj,kkk) * p(+1,-1,-1);
    ap(0,0,0) =
                     Ampp(iii,jjj,kkk) * p(+1,-1,-1);
    ap(-1,1,0) =
                     Am0p(iii-1,jjj+1,kkk) * p( 0,-1,-1)
      +              A00p(iii-1,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj+1,kkk) * p( 0, 0,-1)
      +              A0pp(iii-1,jjj+1,kkk) * p(+1, 0,-1);
    ap(0,1,0) =
                     Am0p(iii,jjj+1,kkk) * p(+1,-1,-1)
      +              Ampp(iii,jjj+1,kkk) * p(+1, 0,-1);
    ap(-1,0,1) =
                     Amp0(iii-1,jjj,kkk+1) * p( 0,-1,-1)
      +              A0p0(iii-1,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii-1,jjj,kkk+1) * p( 0,-1, 0)
      +              A0pp(iii-1,jjj,kkk+1) * p(+1,-1, 0);
    ap(0,0,1) =
                     Amp0(iii,jjj,kkk+1) * p(+1,-1,-1)
      +              Ampp(iii,jjj,kkk+1) * p(+1,-1, 0);
    ap(-1,1,1) =
                     Am00(iii-1,jjj+1,kkk+1) * p( 0,-1,-1)
      +              A000(iii-1,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1)
      +              A0p0(iii-1,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj+1,kkk+1) * p( 0,-1, 0)
      +              A00p(iii-1,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
      +              A0pp(iii-1,jjj+1,kkk+1) * p(+1, 0, 0);
    ap(0,1,1) =
                     Am00(iii,jjj+1,kkk+1) * p(+1,-1,-1)
      +              Amp0(iii,jjj+1,kkk+1) * p(+1, 0,-1)
      +              Am0p(iii,jjj+1,kkk+1) * p(+1,-1, 0)
      +              Ampp(iii,jjj+1,kkk+1) * p(+1, 0, 0);
    cs2 = 0.125_rt *
      ( restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0)
      + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1)
      + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1));

    // alternative: csten(i,j+1,k,ist_pmp)
    iii = ii;
    jjj = jj+2;
    kkk = kk;
    p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj-2,kkk+1);
    p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj-2,kkk+1);
    p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj-1,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj-1,kkk+1);
    p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj-2,kkk+2);
    p( 0, 0, 0) = 1._rt;
    p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj-1,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj-1,kkk+2);
    ap(0,-1,0) =
                     Apmp(iii,jjj-1,kkk) * p(-1, 0,-1)
      +              Ap0p(iii,jjj-1,kkk) * p(-1,+1,-1);
    ap(1,-1,0) =
                     A0mp(iii+1,jjj-1,kkk) * p(-1, 0,-1)
      +              Apmp(iii+1,jjj-1,kkk) * p( 0, 0,-1)
      +              A00p(iii+1,jjj-1,kkk) * p(-1,+1,-1)
      +              Ap0p(iii+1,jjj-1,kkk) * p( 0,+1,-1);
    ap(0,0,0) =
                     Apmp(iii,jjj,kkk) * p(-1,+1,-1);
    ap(1,0,0) =
                     A0mp(iii+1,jjj,kkk) * p(-1,+1,-1)
      +              Apmp(iii+1,jjj,kkk) * p( 0,+1,-1);
    ap(0,-1,1) =
                     Apm0(iii,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Ap00(iii,jjj-1,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Ap0p(iii,jjj-1,kkk+1) * p(-1,+1, 0);
    ap(1,-1,1) =
                     A0m0(iii+1,jjj-1,kkk+1) * p(-1, 0,-1)
      +              Apm0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A000(iii+1,jjj-1,kkk+1) * p(-1,+1,-1)
      +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj-1,kkk+1) * p(-1, 0, 0)
      +              Apmp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A00p(iii+1,jjj-1,kkk+1) * p(-1,+1, 0)
      +              Ap0p(iii+1,jjj-1,kkk+1) * p( 0,+1, 0);
    ap(0,0,1) =
                     Apm0(iii,jjj,kkk+1) * p(-1,+1,-1)
      +              Apmp(iii,jjj,kkk+1) * p(-1,+1, 0);
    ap(1,0,1) =
                     A0m0(iii+1,jjj,kkk+1) * p(-1,+1,-1)
      +              Apm0(iii+1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0mp(iii+1,jjj,kkk+1) * p(-1,+1, 0)
      +              Apmp(iii+1,jjj,kkk+1) * p( 0,+1, 0);
    cs3 = 0.125_rt *
      ( restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1)
      + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1));

    // alternative: csten(i+1,j+1,k,ist_mmp)
    iii = ii+2;
    jjj = jj+2;
    kkk = kk;
    p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj-2,kkk+1);
    p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj-2,kkk+1);
    p( 0,+1,-1) = interp_from_0mp_to(iii-2,jjj-1,kkk+1);
    p(+1,+1,-1) = interp_from_mmp_to(iii-1,jjj-1,kkk+1);
    p( 0, 0, 0) = 1._rt;
    p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj-2,kkk+2);
    p( 0,+1, 0) = interp_from_0m0_to(iii-2,jjj-1,kkk+2);
    p(+1,+1, 0) = interp_from_mm0_to(iii-1,jjj-1,kkk+2);
    ap(-1,-1,0) =
                     Ammp(iii-1,jjj-1,kkk) * p( 0, 0,-1)
      +              A0mp(iii-1,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii-1,jjj-1,kkk) * p( 0,+1,-1)
      +              A00p(iii-1,jjj-1,kkk) * p(+1,+1,-1);
    ap(0,-1,0) =
                     Ammp(iii,jjj-1,kkk) * p(+1, 0,-1)
      +              Am0p(iii,jjj-1,kkk) * p(+1,+1,-1);
    ap(-1,0,0) =
                     Ammp(iii-1,jjj,kkk) * p( 0,+1,-1)
      +              A0mp(iii-1,jjj,kkk) * p(+1,+1,-1);
    ap(0,0,0) =
                     Ammp(iii,jjj,kkk) * p(+1,+1,-1);
    ap(-1,-1,1) =
                     Amm0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1)
      +              A0m0(iii-1,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii-1,jjj-1,kkk+1) * p( 0,+1,-1)
      +              A000(iii-1,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0)
      +              A0mp(iii-1,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii-1,jjj-1,kkk+1) * p( 0,+1, 0)
      +              A00p(iii-1,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(0,-1,1) =
                     Amm0(iii,jjj-1,kkk+1) * p(+1, 0,-1)
      +              Am00(iii,jjj-1,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj-1,kkk+1) * p(+1, 0, 0)
      +              Am0p(iii,jjj-1,kkk+1) * p(+1,+1, 0);
    ap(-1,0,1) =
                     Amm0(iii-1,jjj,kkk+1) * p( 0,+1,-1)
      +              A0m0(iii-1,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii-1,jjj,kkk+1) * p( 0,+1, 0)
      +              A0mp(iii-1,jjj,kkk+1) * p(+1,+1, 0);
    ap(0,0,1) =
                     Amm0(iii,jjj,kkk+1) * p(+1,+1,-1)
      +              Ammp(iii,jjj,kkk+1) * p(+1,+1, 0);
    cs4 = 0.125_rt *
      ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0)
      + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0)
      + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0)
      + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0)
      + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1)
      + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1)
      + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1)
      + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1));

    csten(i,j,k,ist_ppp) = 0.25_rt*(cs1+cs2+cs3+cs4);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_adotx_sten (int i, int j, int k, Array4<Real> const& y, Array4<Real const> const& x,
                         Array4<Real const> const& sten, Array4<int const> const& msk) noexcept
{
    if (msk(i,j,k)) {
        y(i,j,k) = 0.0;
    } else {
        y(i,j,k) = x(i  ,j  ,k  ) * sten(i  ,j  ,k  ,ist_000)
            //
            +      x(i-1,j  ,k  ) * sten(i-1,j  ,k  ,ist_p00)
            +      x(i+1,j  ,k  ) * sten(i  ,j  ,k  ,ist_p00)
            //
            +      x(i  ,j-1,k  ) * sten(i  ,j-1,k  ,ist_0p0)
            +      x(i  ,j+1,k  ) * sten(i  ,j  ,k  ,ist_0p0)
            //
            +      x(i  ,j  ,k-1) * sten(i  ,j  ,k-1,ist_00p)
            +      x(i  ,j  ,k+1) * sten(i  ,j  ,k  ,ist_00p)
            //
            +      x(i-1,j-1,k  ) * sten(i-1,j-1,k  ,ist_pp0)
            +      x(i+1,j-1,k  ) * sten(i  ,j-1,k  ,ist_pp0)
            +      x(i-1,j+1,k  ) * sten(i-1,j  ,k  ,ist_pp0)
            +      x(i+1,j+1,k  ) * sten(i  ,j  ,k  ,ist_pp0)
            //
            +      x(i-1,j  ,k-1) * sten(i-1,j  ,k-1,ist_p0p)
            +      x(i+1,j  ,k-1) * sten(i  ,j  ,k-1,ist_p0p)
            +      x(i-1,j  ,k+1) * sten(i-1,j  ,k  ,ist_p0p)
            +      x(i+1,j  ,k+1) * sten(i  ,j  ,k  ,ist_p0p)
            //
            +      x(i  ,j-1,k-1) * sten(i  ,j-1,k-1,ist_0pp)
            +      x(i  ,j+1,k-1) * sten(i  ,j  ,k-1,ist_0pp)
            +      x(i  ,j-1,k+1) * sten(i  ,j-1,k  ,ist_0pp)
            +      x(i  ,j+1,k+1) * sten(i  ,j  ,k  ,ist_0pp)
            //
            +      x(i-1,j-1,k-1) * sten(i-1,j-1,k-1,ist_ppp)
            +      x(i+1,j-1,k-1) * sten(i  ,j-1,k-1,ist_ppp)
            +      x(i-1,j+1,k-1) * sten(i-1,j  ,k-1,ist_ppp)
            +      x(i+1,j+1,k-1) * sten(i  ,j  ,k-1,ist_ppp)
            +      x(i-1,j-1,k+1) * sten(i-1,j-1,k  ,ist_ppp)
            +      x(i+1,j-1,k+1) * sten(i  ,j-1,k  ,ist_ppp)
            +      x(i-1,j+1,k+1) * sten(i-1,j  ,k  ,ist_ppp)
            +      x(i+1,j+1,k+1) * sten(i  ,j  ,k  ,ist_ppp);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_gauss_seidel_sten (Box const& bx, Array4<Real> const& sol,
                                Array4<Real const> const& rhs,
                                Array4<Real const> const& sten,
                                Array4<int const> const& msk) noexcept
{
    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else if (sten(i,j,k,ist_000) != 0.0) {
            Real Ax  = sol(i  ,j  ,k  ) * sten(i  ,j  ,k  ,ist_000)
                //
                +      sol(i-1,j  ,k  ) * sten(i-1,j  ,k  ,ist_p00)
                +      sol(i+1,j  ,k  ) * sten(i  ,j  ,k  ,ist_p00)
                //
                +      sol(i  ,j-1,k  ) * sten(i  ,j-1,k  ,ist_0p0)
                +      sol(i  ,j+1,k  ) * sten(i  ,j  ,k  ,ist_0p0)
                //
                +      sol(i  ,j  ,k-1) * sten(i  ,j  ,k-1,ist_00p)
                +      sol(i  ,j  ,k+1) * sten(i  ,j  ,k  ,ist_00p)
                //
                +      sol(i-1,j-1,k  ) * sten(i-1,j-1,k  ,ist_pp0)
                +      sol(i+1,j-1,k  ) * sten(i  ,j-1,k  ,ist_pp0)
                +      sol(i-1,j+1,k  ) * sten(i-1,j  ,k  ,ist_pp0)
                +      sol(i+1,j+1,k  ) * sten(i  ,j  ,k  ,ist_pp0)
                //
                +      sol(i-1,j  ,k-1) * sten(i-1,j  ,k-1,ist_p0p)
                +      sol(i+1,j  ,k-1) * sten(i  ,j  ,k-1,ist_p0p)
                +      sol(i-1,j  ,k+1) * sten(i-1,j  ,k  ,ist_p0p)
                +      sol(i+1,j  ,k+1) * sten(i  ,j  ,k  ,ist_p0p)
                //
                +      sol(i  ,j-1,k-1) * sten(i  ,j-1,k-1,ist_0pp)
                +      sol(i  ,j+1,k-1) * sten(i  ,j  ,k-1,ist_0pp)
                +      sol(i  ,j-1,k+1) * sten(i  ,j-1,k  ,ist_0pp)
                +      sol(i  ,j+1,k+1) * sten(i  ,j  ,k  ,ist_0pp)
                //
                +      sol(i-1,j-1,k-1) * sten(i-1,j-1,k-1,ist_ppp)
                +      sol(i+1,j-1,k-1) * sten(i  ,j-1,k-1,ist_ppp)
                +      sol(i-1,j+1,k-1) * sten(i-1,j  ,k-1,ist_ppp)
                +      sol(i+1,j+1,k-1) * sten(i  ,j  ,k-1,ist_ppp)
                +      sol(i-1,j-1,k+1) * sten(i-1,j-1,k  ,ist_ppp)
                +      sol(i+1,j-1,k+1) * sten(i  ,j-1,k  ,ist_ppp)
                +      sol(i-1,j+1,k+1) * sten(i-1,j  ,k  ,ist_ppp)
                +      sol(i+1,j+1,k+1) * sten(i  ,j  ,k  ,ist_ppp);

            sol(i,j,k) += (rhs(i,j,k) - Ax) / sten(i,j,k,ist_000);
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_interpadd_rap (int i, int j, int k, Array4<Real> const& fine,
                            Array4<Real const> const& crse, Array4<Real const> const& sten,
                            Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k) and sten(i,j,k,ist_000) != 0.0) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        int kc = amrex::coarsen(k,2);
        bool ieven = ic*2 == i;
        bool jeven = jc*2 == j;
        bool keven = kc*2 == k;
        Real fv;
        if (ieven and jeven and keven) {
            fv = crse(ic,jc,kc);
        } else if (ieven and jeven) {
            Real w1 = std::abs(sten(i,j,k-1,ist_00p));
            Real w2 = std::abs(sten(i,j,k  ,ist_00p));
            if (w1 == 0.0 and w2 == 0.0) {
                fv = 0.5*(crse(ic,jc,kc)+crse(ic,jc,kc+1));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic,jc,kc+1)) / (w1+w2);
            }
        } else if (ieven and keven) {
            Real w1 = std::abs(sten(i,j-1,k,ist_0p0));
            Real w2 = std::abs(sten(i,j  ,k,ist_0p0));
            if (w1 == 0.0 and w2 == 0.0) {
                fv = 0.5*(crse(ic,jc,kc)+crse(ic,jc+1,kc));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic,jc+1,kc)) / (w1+w2);
            }
        } else if (jeven and keven) {
            Real w1 = std::abs(sten(i-1,j,k,ist_p00));
            Real w2 = std::abs(sten(i  ,j,k,ist_p00));
            if (w1 == 0.0 and w2 == 0.0) {
                fv = 0.5*(crse(ic,jc,kc)+crse(ic+1,jc,kc));
            } else {
                fv = (w1*crse(ic,jc,kc) + w2*crse(ic+1,jc,kc)) / (w1+w2);
            }
        } else if (ieven) {
            Real w1m = std::abs(sten(i,j-1,k,ist_0p0)) / (std::abs(sten(i,j-1,k-1,ist_0pp))
                                                         +std::abs(sten(i,j-1,k  ,ist_0pp)) + eps);
            Real w1p = std::abs(sten(i,j  ,k,ist_0p0)) / (std::abs(sten(i,j  ,k-1,ist_0pp))
                                                         +std::abs(sten(i,j  ,k  ,ist_0pp)) + eps);
            Real w2m = std::abs(sten(i,j,k-1,ist_00p)) / (std::abs(sten(i,j-1,k-1,ist_0pp))
                                                         +std::abs(sten(i,j  ,k-1,ist_0pp)) + eps);
            Real w2p = std::abs(sten(i,j,k  ,ist_00p)) / (std::abs(sten(i,j-1,k  ,ist_0pp))
                                                         +std::abs(sten(i,j  ,k  ,ist_0pp)) + eps);
            Real wmm = std::abs(sten(i,j-1,k-1,ist_0pp)) * (1.0 + w1m + w2m);
            Real wpm = std::abs(sten(i,j  ,k-1,ist_0pp)) * (1.0 + w1p + w2m);
            Real wmp = std::abs(sten(i,j-1,k  ,ist_0pp)) * (1.0 + w1m + w2p);
            Real wpp = std::abs(sten(i,j  ,k  ,ist_0pp)) * (1.0 + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic,jc+1,kc)
                  + wmp*crse(ic,jc,kc+1) + wpp*crse(ic,jc+1,kc+1))
                / (wmm+wpm+wmp+wpp+eps);
        } else if (jeven) {
            Real w1m = std::abs(sten(i-1,j,k,ist_p00)) / (std::abs(sten(i-1,j,k-1,ist_p0p))
                                                         +std::abs(sten(i-1,j,k  ,ist_p0p)) + eps);
            Real w1p = std::abs(sten(i  ,j,k,ist_p00)) / (std::abs(sten(i  ,j,k-1,ist_p0p))
                                                         +std::abs(sten(i  ,j,k  ,ist_p0p)) + eps);
            Real w2m = std::abs(sten(i,j,k-1,ist_00p)) / (std::abs(sten(i-1,j,k-1,ist_p0p))
                                                         +std::abs(sten(i  ,j,k-1,ist_p0p)) + eps);
            Real w2p = std::abs(sten(i,j,k  ,ist_00p)) / (std::abs(sten(i-1,j,k  ,ist_p0p))
                                                         +std::abs(sten(i  ,j,k  ,ist_p0p)) + eps);
            Real wmm = std::abs(sten(i-1,j,k-1,ist_p0p)) * (1.0 + w1m + w2m);
            Real wpm = std::abs(sten(i  ,j,k-1,ist_p0p)) * (1.0 + w1p + w2m);
            Real wmp = std::abs(sten(i-1,j,k  ,ist_p0p)) * (1.0 + w1m + w2p);
            Real wpp = std::abs(sten(i  ,j,k  ,ist_p0p)) * (1.0 + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic+1,jc,kc)
                  + wmp*crse(ic,jc,kc+1) + wpp*crse(ic+1,jc,kc+1))
                / (wmm+wpm+wmp+wpp+eps);
        } else if (keven) {
            Real w1m = std::abs(sten(i-1,j,k,ist_p00)) / (std::abs(sten(i-1,j-1,k,ist_pp0))
                                                         +std::abs(sten(i-1,j  ,k,ist_pp0)) + eps);
            Real w1p = std::abs(sten(i  ,j,k,ist_p00)) / (std::abs(sten(i  ,j-1,k,ist_pp0))
                                                         +std::abs(sten(i  ,j  ,k,ist_pp0)) + eps);
            Real w2m = std::abs(sten(i,j-1,k,ist_0p0)) / (std::abs(sten(i-1,j-1,k,ist_pp0))
                                                         +std::abs(sten(i  ,j-1,k,ist_pp0)) + eps);
            Real w2p = std::abs(sten(i,j  ,k,ist_0p0)) / (std::abs(sten(i-1,j  ,k,ist_pp0))
                                                         +std::abs(sten(i  ,j  ,k,ist_pp0)) + eps);
            Real wmm = std::abs(sten(i-1,j-1,k,ist_pp0)) * (1.0 + w1m + w2m);
            Real wpm = std::abs(sten(i  ,j-1,k,ist_pp0)) * (1.0 + w1p + w2m);
            Real wmp = std::abs(sten(i-1,j  ,k,ist_pp0)) * (1.0 + w1m + w2p);
            Real wpp = std::abs(sten(i  ,j  ,k,ist_pp0)) * (1.0 + w1p + w2p);
            fv = (wmm*crse(ic,jc,kc) + wpm*crse(ic+1,jc,kc)
                  + wmp*crse(ic,jc+1,kc) + wpp*crse(ic+1,jc+1,kc))
                / (wmm+wpm+wmp+wpp+eps);
        } else {
            Real wmmm = 1.0;
            Real wpmm = 1.0;
            Real wmpm = 1.0;
            Real wppm = 1.0;
            Real wmmp = 1.0;
            Real wpmp = 1.0;
            Real wmpp = 1.0;
            Real wppp = 1.0;

            Real wtmp = std::abs(sten(i-1,j,k,ist_p00)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j  ,k-1,ist_ppp))
                + std::abs(sten(i-1,j-1,k  ,ist_ppp))
                + std::abs(sten(i-1,j  ,k  ,ist_ppp)) + eps);
            wmmm += wtmp;
            wmpm += wtmp;
            wmmp += wtmp;
            wmpp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_p00)) /
                ( std::abs(sten(i,j-1,k-1,ist_ppp))
                + std::abs(sten(i,j  ,k-1,ist_ppp))
                + std::abs(sten(i,j-1,k  ,ist_ppp))
                + std::abs(sten(i,j  ,k  ,ist_ppp)) + eps);
            wpmm += wtmp;
            wppm += wtmp;
            wpmp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i,j-1,k,ist_0p0)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i  ,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j-1,k  ,ist_ppp))
                + std::abs(sten(i  ,j-1,k  ,ist_ppp)) + eps);
            wmmm += wtmp;
            wpmm += wtmp;
            wmmp += wtmp;
            wpmp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_0p0)) /
                ( std::abs(sten(i-1,j,k-1,ist_ppp))
                + std::abs(sten(i  ,j,k-1,ist_ppp))
                + std::abs(sten(i-1,j,k  ,ist_ppp))
                + std::abs(sten(i  ,j,k  ,ist_ppp)) + eps);
            wmpm += wtmp;
            wppm += wtmp;
            wmpp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i,j,k-1,ist_00p)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i  ,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j  ,k-1,ist_ppp))
                + std::abs(sten(i  ,j  ,k-1,ist_ppp)) + eps);
            wmmm += wtmp;
            wpmm += wtmp;
            wmpm += wtmp;
            wppm += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_00p)) /
                ( std::abs(sten(i-1,j-1,k,ist_ppp))
                + std::abs(sten(i  ,j-1,k,ist_ppp))
                + std::abs(sten(i-1,j  ,k,ist_ppp))
                + std::abs(sten(i  ,j  ,k,ist_ppp)) + eps);
            wmmp += wtmp;
            wpmp += wtmp;
            wmpp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i-1,j-1,k,ist_pp0)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j-1,k  ,ist_ppp)) + eps);
            wmmm += wtmp;
            wmmp += wtmp;

            wtmp = std::abs(sten(i,j-1,k,ist_pp0)) /
                ( std::abs(sten(i,j-1,k-1,ist_ppp))
                + std::abs(sten(i,j-1,k  ,ist_ppp)) + eps);
            wpmm += wtmp;
            wpmp += wtmp;

            wtmp = std::abs(sten(i-1,j,k,ist_pp0)) /
                ( std::abs(sten(i-1,j,k-1,ist_ppp))
                + std::abs(sten(i-1,j,k  ,ist_ppp)) + eps);
            wmpm += wtmp;
            wmpp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_pp0)) /
                ( std::abs(sten(i,j,k-1,ist_ppp))
                + std::abs(sten(i,j,k  ,ist_ppp)) + eps);
            wppm += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i-1,j,k-1,ist_p0p)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i-1,j  ,k-1,ist_ppp)) + eps);
            wmmm += wtmp;
            wmpm += wtmp;

            wtmp = std::abs(sten(i,j,k-1,ist_p0p)) /
                ( std::abs(sten(i,j-1,k-1,ist_ppp))
                + std::abs(sten(i,j  ,k-1,ist_ppp)) + eps);
            wpmm += wtmp;
            wppm += wtmp;

            wtmp = std::abs(sten(i-1,j,k,ist_p0p)) /
                ( std::abs(sten(i-1,j-1,k,ist_ppp))
                + std::abs(sten(i-1,j  ,k,ist_ppp)) + eps);
            wmmp += wtmp;
            wmpp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_p0p)) /
                ( std::abs(sten(i,j-1,k,ist_ppp))
                + std::abs(sten(i,j  ,k,ist_ppp)) + eps);
            wpmp += wtmp;
            wppp += wtmp;

            wtmp = std::abs(sten(i,j-1,k-1,ist_0pp)) /
                ( std::abs(sten(i-1,j-1,k-1,ist_ppp))
                + std::abs(sten(i  ,j-1,k-1,ist_ppp)) + eps);
            wmmm += wtmp;
            wpmm += wtmp;

            wtmp = std::abs(sten(i,j,k-1,ist_0pp)) /
                ( std::abs(sten(i-1,j,k-1,ist_ppp))
                + std::abs(sten(i  ,j,k-1,ist_ppp)) + eps);
            wmpm += wtmp;
            wppm += wtmp;

            wtmp = std::abs(sten(i,j-1,k,ist_0pp)) /
                ( std::abs(sten(i-1,j-1,k,ist_ppp))
                + std::abs(sten(i  ,j-1,k,ist_ppp)) + eps);
            wmmp += wtmp;
            wpmp += wtmp;

            wtmp = std::abs(sten(i,j,k,ist_0pp)) /
                ( std::abs(sten(i-1,j,k,ist_ppp))
                + std::abs(sten(i  ,j,k,ist_ppp)) + eps);
            wmpp += wtmp;
            wppp += wtmp;

            wmmm *= std::abs(sten(i-1,j-1,k-1,ist_ppp));
            wpmm *= std::abs(sten(i  ,j-1,k-1,ist_ppp));
            wmpm *= std::abs(sten(i-1,j  ,k-1,ist_ppp));
            wppm *= std::abs(sten(i  ,j  ,k-1,ist_ppp));
            wmmp *= std::abs(sten(i-1,j-1,k  ,ist_ppp));
            wpmp *= std::abs(sten(i  ,j-1,k  ,ist_ppp));
            wmpp *= std::abs(sten(i-1,j  ,k  ,ist_ppp));
            wppp *= std::abs(sten(i  ,j  ,k  ,ist_ppp));
            fv = (wmmm*crse(ic,jc  ,kc  ) + wpmm*crse(ic+1,jc  ,kc  )
                  + wmpm*crse(ic,jc+1,kc  ) + wppm*crse(ic+1,jc+1,kc  )
                  + wmmp*crse(ic,jc  ,kc+1) + wpmp*crse(ic+1,jc  ,kc+1)
                  + wmpp*crse(ic,jc+1,kc+1) + wppp*crse(ic+1,jc+1,kc+1))
                / (wmmm + wpmm + wmpm + wppm + wmmp + wpmp + wmpp + wppp + eps);
        }

        fine(i,j,k) += fv;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_restriction_rap (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine, Array4<Real const> const& sten,
                              Array4<int const> const& msk) noexcept
{
    int ii = i*2;
    int jj = j*2;
    int kk = k*2;
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = 0.0;
    } else {

        Real cv = fine(ii,jj,kk);

        // ************************************
        // Adding fine(ii-1,jj,kk)
        // ************************************

        Real sten_lo = std::abs(sten(ii-2,jj,kk,ist_p00));
        Real sten_hi = std::abs(sten(ii-1,jj,kk,ist_p00));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii-1,jj,kk);
        } else {
            cv += fine(ii-1,jj,kk) * sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii+1,jj,kk)
        // ************************************

        sten_lo = std::abs(sten(ii  ,jj,kk,ist_p00));
        sten_hi = std::abs(sten(ii+1,jj,kk,ist_p00));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii+1,jj,kk);
        } else {
            cv += fine(ii+1,jj,kk) * sten_lo / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj-1,kk)
        // ************************************

        sten_lo = std::abs(sten(ii,jj-2,kk,ist_0p0));
        sten_hi = std::abs(sten(ii,jj-1,kk,ist_0p0));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj-1,kk);
        } else {
            cv += fine(ii,jj-1,kk) * sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj+1,kk)
        // ************************************

        sten_lo = std::abs(sten(ii,jj  ,kk,ist_0p0));
        sten_hi = std::abs(sten(ii,jj+1,kk,ist_0p0));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj+1,kk);
        } else {
            cv += fine(ii,jj+1,kk) * sten_lo / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj,kk-1)
        // ************************************

        sten_lo = std::abs(sten(ii,jj,kk-2,ist_00p));
        sten_hi = std::abs(sten(ii,jj,kk-1,ist_00p));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj,kk-1);
        } else {
            cv += fine(ii,jj,kk-1)*sten_hi / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii,jj,kk+1)
        // ************************************

        sten_lo = std::abs(sten(ii,jj,kk  ,ist_00p));
        sten_hi = std::abs(sten(ii,jj,kk+1,ist_00p));

        if (sten_lo == 0.0 and sten_hi == 0.0) {
            cv += 0.5*fine(ii,jj,kk+1);
        } else {
            cv += fine(ii,jj,kk+1)*sten_lo  / (sten_lo + sten_hi);
        }

        // ************************************
        // Adding fine(ii-1,jj-1,kk)
        // ************************************

        // keven
        Real w1m = std::abs(sten(ii-2,jj-1,kk,ist_p00))
            / (    std::abs(sten(ii-2,jj-2,kk,ist_pp0))
                  +std::abs(sten(ii-2,jj-1,kk,ist_pp0)) + eps);
        Real w1p = std::abs(sten(ii-1,jj-1,kk,ist_p00))
            / (    std::abs(sten(ii-1,jj-2,kk,ist_pp0))
                  +std::abs(sten(ii-1,jj-1,kk,ist_pp0)) + eps);
        Real w2m = std::abs(sten(ii-1,jj-2,kk,ist_0p0))
            / (    std::abs(sten(ii-2,jj-2,kk,ist_pp0))
                  +std::abs(sten(ii-1,jj-2,kk,ist_pp0)) + eps);
        Real w2p = std::abs(sten(ii-1,jj-1,kk,ist_0p0))
            / (    std::abs(sten(ii-2,jj-1,kk,ist_pp0))
                  +std::abs(sten(ii-1,jj-1,kk,ist_pp0)) + eps);
        Real wmm = std::abs(sten(ii-2,jj-2,kk,ist_pp0)) * (1.0 + w1m + w2m);
        Real wpm = std::abs(sten(ii-1,jj-2,kk,ist_pp0)) * (1.0 + w1p + w2m);
        Real wmp = std::abs(sten(ii-2,jj-1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        Real wpp = std::abs(sten(ii-1,jj-1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj-1,kk)*wpp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj-1,kk)
        // ************************************

        w1m = std::abs(sten(ii  ,jj-1,kk,ist_p00))
           / (std::abs(sten(ii  ,jj-2,kk,ist_pp0))
             +std::abs(sten(ii  ,jj-1,kk,ist_pp0)) + eps);
        w1p = std::abs(sten(ii+1,jj-1,kk,ist_p00))
           / (std::abs(sten(ii+1,jj-2,kk,ist_pp0))
             +std::abs(sten(ii+1,jj-1,kk,ist_pp0)) + eps);
        w2m = std::abs(sten(ii+1,jj-2,kk,ist_0p0))
           / (std::abs(sten(ii  ,jj-2,kk,ist_pp0))
             +std::abs(sten(ii+1,jj-2,kk,ist_pp0)) + eps);
        w2p = std::abs(sten(ii+1,jj-1,kk,ist_0p0))
           / (std::abs(sten(ii  ,jj-1,kk,ist_pp0))
             +std::abs(sten(ii+1,jj-1,kk,ist_pp0)) + eps);
        wmm = std::abs(sten(ii  ,jj-2,kk,ist_pp0)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj-2,kk,ist_pp0)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj-1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj-1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj-1,kk)*wmp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii-1,jj+1,kk)
        // ************************************

        w1m = std::abs(sten(ii-2,jj+1,kk,ist_p00))
           / (std::abs(sten(ii-2,jj  ,kk,ist_pp0))
             +std::abs(sten(ii-2,jj+1,kk,ist_pp0)) + eps);
        w1p = std::abs(sten(ii-1,jj+1,kk,ist_p00))
           / (std::abs(sten(ii-1,jj  ,kk,ist_pp0))
             +std::abs(sten(ii-1,jj+1,kk,ist_pp0)) + eps);
        w2m = std::abs(sten(ii-1,jj  ,kk,ist_0p0))
           / (std::abs(sten(ii-2,jj  ,kk,ist_pp0))
             +std::abs(sten(ii-1,jj  ,kk,ist_pp0)) + eps);
        w2p = std::abs(sten(ii-1,jj+1,kk,ist_0p0))
           / (std::abs(sten(ii-2,jj+1,kk,ist_pp0))
             +std::abs(sten(ii-1,jj+1,kk,ist_pp0)) + eps);
        wmm = std::abs(sten(ii-2,jj  ,kk,ist_pp0)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii-1,jj  ,kk,ist_pp0)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii-2,jj+1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii-1,jj+1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj+1,kk)*wpm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj+1,kk)
        // ************************************

        w1m = std::abs(sten(ii  ,jj+1,kk,ist_p00))
           / (std::abs(sten(ii  ,jj+1,kk,ist_pp0))
             +std::abs(sten(ii  ,jj  ,kk,ist_pp0)) + eps);
        w1p = std::abs(sten(ii+1,jj+1,kk,ist_p00))
           / (std::abs(sten(ii+1,jj+1,kk,ist_pp0))
             +std::abs(sten(ii+1,jj  ,kk,ist_pp0)) + eps);
        w2m = std::abs(sten(ii+1,jj  ,kk,ist_0p0))
           / (std::abs(sten(ii+1,jj  ,kk,ist_pp0))
             +std::abs(sten(ii  ,jj  ,kk,ist_pp0)) + eps);
        w2p = std::abs(sten(ii+1,jj+1,kk,ist_0p0))
           / (std::abs(sten(ii+1,jj+1,kk,ist_pp0))
             +std::abs(sten(ii  ,jj+1,kk,ist_pp0)) + eps);
        wmm = std::abs(sten(ii  ,jj  ,kk,ist_pp0)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj  ,kk,ist_pp0)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj+1,kk,ist_pp0)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj+1,kk,ist_pp0)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj+1,kk)*wmm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii-1,jj,kk-1)
        // ************************************

        // jeven
        w1m = std::abs(sten(ii-2,jj,kk-1,ist_p00))
           / (std::abs(sten(ii-2,jj,kk-2,ist_p0p))
             +std::abs(sten(ii-2,jj,kk-1,ist_p0p)) + eps);
        w1p = std::abs(sten(ii-1,jj,kk-1,ist_p00))
           / (std::abs(sten(ii-1,jj,kk-2,ist_p0p))
             +std::abs(sten(ii-1,jj,kk-1,ist_p0p)) + eps);
        w2m = std::abs(sten(ii-1,jj,kk-2,ist_00p))
           / (std::abs(sten(ii-2,jj,kk-2,ist_p0p))
             +std::abs(sten(ii-1,jj,kk-2,ist_p0p)) + eps);
        w2p = std::abs(sten(ii-1,jj,kk-1,ist_00p))
           / (std::abs(sten(ii-2,jj,kk-1,ist_p0p))
             +std::abs(sten(ii-1,jj,kk-1,ist_p0p)) + eps);
        wmm = std::abs(sten(ii-2,jj,kk-2,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii-1,jj,kk-2,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii-2,jj,kk-1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii-1,jj,kk-1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj,kk-1)*wpp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj,kk-1)
        // ************************************

        w1m = std::abs(sten(ii  ,jj,kk-1,ist_p00))
           / (std::abs(sten(ii  ,jj,kk-2,ist_p0p))
             +std::abs(sten(ii  ,jj,kk-1,ist_p0p)) + eps);
        w1p = std::abs(sten(ii+1,jj,kk-1,ist_p00))
           / (std::abs(sten(ii+1,jj,kk-2,ist_p0p))
             +std::abs(sten(ii+1,jj,kk-1,ist_p0p)) + eps);
        w2m = std::abs(sten(ii+1,jj,kk-2,ist_00p))
           / (std::abs(sten(ii+1,jj,kk-2,ist_p0p))
             +std::abs(sten(ii  ,jj,kk-2,ist_p0p)) + eps);
        w2p = std::abs(sten(ii+1,jj,kk-1,ist_00p))
           / (std::abs(sten(ii+1,jj,kk-1,ist_p0p))
             +std::abs(sten(ii  ,jj,kk-1,ist_p0p)) + eps);
        wmm = std::abs(sten(ii  ,jj,kk-2,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj,kk-2,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj,kk-1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj,kk-1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj,kk-1)*wmp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii-1,jj,kk+1)
        // ************************************

        w1m = std::abs(sten(ii-2,jj,kk+1,ist_p00))
           / (std::abs(sten(ii-2,jj,kk+1,ist_p0p))
             +std::abs(sten(ii-2,jj,kk  ,ist_p0p)) + eps);
        w1p = std::abs(sten(ii-1,jj,kk+1,ist_p00))
           / (std::abs(sten(ii-1,jj,kk+1,ist_p0p))
             +std::abs(sten(ii-1,jj,kk  ,ist_p0p)) + eps);
        w2m = std::abs(sten(ii-1,jj,kk  ,ist_00p))
           / (std::abs(sten(ii-2,jj,kk  ,ist_p0p))
             +std::abs(sten(ii-1,jj,kk  ,ist_p0p)) + eps);
        w2p = std::abs(sten(ii-1,jj,kk+1,ist_00p))
           / (std::abs(sten(ii-2,jj,kk+1,ist_p0p))
             +std::abs(sten(ii-1,jj,kk+1,ist_p0p)) + eps);
        wmm = std::abs(sten(ii-2,jj,kk  ,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii-1,jj,kk  ,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii-2,jj,kk+1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii-1,jj,kk+1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii-1,jj,kk+1)*wpm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii+1,jj,kk+1)
        // ************************************

        w1m = std::abs(sten(ii  ,jj,kk+1,ist_p00))
           / (std::abs(sten(ii  ,jj,kk+1,ist_p0p))
             +std::abs(sten(ii  ,jj,kk  ,ist_p0p)) + eps);
        w1p = std::abs(sten(ii+1,jj,kk+1,ist_p00))
           / (std::abs(sten(ii+1,jj,kk+1,ist_p0p))
              +std::abs(sten(ii+1,jj,kk  ,ist_p0p)) + eps);
        w2m = std::abs(sten(ii+1,jj,kk  ,ist_00p))
           / (std::abs(sten(ii+1,jj,kk  ,ist_p0p))
             +std::abs(sten(ii  ,jj,kk  ,ist_p0p)) + eps);
        w2p = std::abs(sten(ii+1,jj,kk+1,ist_00p))
           / (std::abs(sten(ii+1,jj,kk+1,ist_p0p))
             +std::abs(sten(ii  ,jj,kk+1,ist_p0p)) + eps);
        wmm = std::abs(sten(ii  ,jj,kk  ,ist_p0p)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii+1,jj,kk  ,ist_p0p)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii  ,jj,kk+1,ist_p0p)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii+1,jj,kk+1,ist_p0p)) * (1.0 + w1p + w2p);
        cv += fine(ii+1,jj,kk+1)*wmm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj-1,kk-1)
        // ************************************

        // ieven
        w1m = std::abs(sten(ii,jj-2,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj-2,kk-2,ist_0pp))
             +std::abs(sten(ii,jj-2,kk-1,ist_0pp)) + eps);
        w2m = std::abs(sten(ii,jj-1,kk-2,ist_00p))
           / (std::abs(sten(ii,jj-2,kk-2,ist_0pp))
             +std::abs(sten(ii,jj-1,kk-2,ist_0pp)) + eps);
        w1p = std::abs(sten(ii,jj-1,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj-1,kk-2,ist_0pp))
             +std::abs(sten(ii,jj-1,kk-1,ist_0pp)) + eps);
        w2p = std::abs(sten(ii,jj-1,kk-1,ist_00p))
           / (std::abs(sten(ii,jj-2,kk-1,ist_0pp))
             +std::abs(sten(ii,jj-1,kk-1,ist_0pp)) + eps);
        wmm = std::abs(sten(ii,jj-2,kk-2,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii,jj-1,kk-2,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii,jj-2,kk-1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii,jj-1,kk-1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj-1,kk-1)*wpp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj+1,kk-1)
        // ************************************

        w1m = std::abs(sten(ii,jj  ,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj  ,kk-2,ist_0pp))
             +std::abs(sten(ii,jj  ,kk-1,ist_0pp)) + eps);
        w1p = std::abs(sten(ii,jj+1,kk-1,ist_0p0))
           / (std::abs(sten(ii,jj+1,kk-2,ist_0pp))
             +std::abs(sten(ii,jj+1,kk-1,ist_0pp)) + eps);
        w2m = std::abs(sten(ii,jj+1,kk-2,ist_00p))
           / (std::abs(sten(ii,jj+1,kk-2,ist_0pp))
             +std::abs(sten(ii,jj  ,kk-2,ist_0pp)) + eps);
        w2p = std::abs(sten(ii,jj+1,kk-1,ist_00p))
           / (std::abs(sten(ii,jj+1,kk-1,ist_0pp))
             +std::abs(sten(ii,jj  ,kk-1,ist_0pp)) + eps);
        wmm = std::abs(sten(ii,jj  ,kk-2,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii,jj+1,kk-2,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii,jj  ,kk-1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii,jj+1,kk-1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj+1,kk-1)*wmp/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj-1,kk+1)
        // ************************************

        w1m = std::abs(sten(ii,jj-2,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj-2,kk+1,ist_0pp))
             +std::abs(sten(ii,jj-2,kk  ,ist_0pp)) + eps);
        w1p = std::abs(sten(ii,jj-1,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj-1,kk+1,ist_0pp))
             +std::abs(sten(ii,jj-1,kk  ,ist_0pp)) + eps);
        w2m = std::abs(sten(ii,jj-1,kk  ,ist_00p))
           / (std::abs(sten(ii,jj-2,kk  ,ist_0pp))
             +std::abs(sten(ii,jj-1,kk  ,ist_0pp)) + eps);
        w2p = std::abs(sten(ii,jj-1,kk+1,ist_00p))
           / (std::abs(sten(ii,jj-2,kk+1,ist_0pp))
             +std::abs(sten(ii,jj-1,kk+1,ist_0pp)) + eps);
        wmm = std::abs(sten(ii,jj-2,kk  ,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii,jj-1,kk  ,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii,jj-2,kk+1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii,jj-1,kk+1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj-1,kk+1)*wpm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine(ii,jj+1,kk+1)
        // ************************************

        w1m = std::abs(sten(ii,jj  ,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj  ,kk+1,ist_0pp))
             +std::abs(sten(ii,jj  ,kk  ,ist_0pp)) + eps);
        w1p = std::abs(sten(ii,jj+1,kk+1,ist_0p0))
           / (std::abs(sten(ii,jj+1,kk+1,ist_0pp))
             +std::abs(sten(ii,jj+1,kk  ,ist_0pp)) + eps);
        w2m = std::abs(sten(ii,jj+1,kk  ,ist_00p))
           / (std::abs(sten(ii,jj+1,kk  ,ist_0pp))
             +std::abs(sten(ii,jj  ,kk  ,ist_0pp)) + eps);
        w2p = std::abs(sten(ii,jj+1,kk+1,ist_00p))
           / (std::abs(sten(ii,jj+1,kk+1,ist_0pp))
             +std::abs(sten(ii,jj  ,kk+1,ist_0pp)) + eps);
        wmm = std::abs(sten(ii,jj  ,kk  ,ist_0pp)) * (1.0 + w1m + w2m);
        wpm = std::abs(sten(ii,jj+1,kk  ,ist_0pp)) * (1.0 + w1p + w2m);
        wmp = std::abs(sten(ii,jj  ,kk+1,ist_0pp)) * (1.0 + w1m + w2p);
        wpp = std::abs(sten(ii,jj+1,kk+1,ist_0pp)) * (1.0 + w1p + w2p);
        cv += fine(ii,jj+1,kk+1)*wmm/(wmm+wpm+wmp+wpp+eps);

        // ************************************
        // Adding fine at corners
        // ************************************

        Real wmmm = 1.0
            + std::abs(sten(ii  ,jj+1,kk+1,ist_p00)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj  ,kk+1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj+1,kk  ,ist_00p)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj+1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj  ,kk+1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj+1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj  ,kk  ,ist_0pp)) /
            ( std::abs(sten(ii  ,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk  ,ist_ppp)) + eps);
        wmmm *= std::abs(sten(ii,jj,kk,ist_ppp));
        cv += wmmm*fine(ii+1,jj+1,kk+1)*sten(ii+1,jj+1,kk+1,ist_inv);

        Real wpmm = 1.0
            + std::abs(sten(ii-1,jj+1,kk+1,ist_p00)) /
            ( std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk+1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj  ,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj+1,kk  ,ist_00p)) /
            ( std::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj+1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk+1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj+1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii-1,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk  ,ist_0pp)) /
            ( std::abs(sten(ii-2,jj  ,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk  ,ist_ppp)) + eps);
        wpmm *= std::abs(sten(ii-1,jj,kk,ist_ppp));
        cv += wpmm*fine(ii-1,jj+1,kk+1)*sten(ii-1,jj+1,kk+1,ist_inv);

        Real wmpm = 1.0
            + std::abs(sten(ii  ,jj-1,kk+1,ist_p00)) /
            ( std::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-2,kk+1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk+1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk+1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk  ,ist_00p)) /
            ( std::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj-1,kk+1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj-1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii  ,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk  ,ist_0pp)) /
            ( std::abs(sten(ii  ,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk  ,ist_ppp)) + eps);
        wmpm *= std::abs(sten(ii  ,jj-1,kk  ,ist_ppp));
        cv += wmpm*fine(ii+1,jj-1,kk+1)*sten(ii+1,jj-1,kk+1,ist_inv);

        Real wppm = 1.0
            + std::abs(sten(ii-1,jj-1,kk+1,ist_p00)) /
            ( std::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk+1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk+1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk  ,ist_00p)) /
            ( std::abs(sten(ii-2,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk+1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk  ,ist_p0p)) /
            ( std::abs(sten(ii-1,jj-2,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk  ,ist_0pp)) /
            ( std::abs(sten(ii-2,jj-1,kk  ,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps);
        wppm *= std::abs(sten(ii-1,jj-1,kk  ,ist_ppp));
        cv += wppm*fine(ii-1,jj-1,kk+1)*sten(ii-1,jj-1,kk+1,ist_inv);

        Real wmmp = 1.0
            + std::abs(sten(ii  ,jj+1,kk-1,ist_p00)) /
            ( std::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj  ,kk-1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj+1,kk-1,ist_00p)) /
            ( std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj  ,kk-1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj  ,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj+1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj  ,kk-1,ist_0pp)) /
            ( std::abs(sten(ii  ,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj  ,kk-1,ist_ppp)) + eps);
        wmmp *= std::abs(sten(ii  ,jj  ,kk-1,ist_ppp));
        cv += wmmp*fine(ii+1,jj+1,kk-1)*sten(ii+1,jj+1,kk-1,ist_inv);

        Real wpmp = 1.0
            + std::abs(sten(ii-1,jj+1,kk-1,ist_p00)) /
            ( std::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk-1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj+1,kk-1,ist_00p)) /
            ( std::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-2,jj+1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk-1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj  ,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj+1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii-1,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj  ,kk-1,ist_0pp)) /
            ( std::abs(sten(ii-2,jj  ,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps);
        wpmp *= std::abs(sten(ii-1,jj  ,kk-1,ist_ppp));
        cv += wpmp*fine(ii-1,jj+1,kk-1)*sten(ii-1,jj+1,kk-1,ist_inv);

        Real wmpp = 1.0
            + std::abs(sten(ii  ,jj-1,kk-1,ist_p00)) /
            ( std::abs(sten(ii  ,jj-2,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk-1,ist_0p0)) /
            ( std::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk-1,ist_00p)) /
            ( std::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj-1,kk-1,ist_pp0)) /
            ( std::abs(sten(ii  ,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii  ,jj-1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii  ,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii+1,jj-1,kk-1,ist_0pp)) /
            ( std::abs(sten(ii  ,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps);
        wmpp *= std::abs(sten(ii  ,jj-1,kk-1,ist_ppp));
        cv += wmpp*fine(ii+1,jj-1,kk-1)*sten(ii+1,jj-1,kk-1,ist_inv);

        Real wppp = 1.0
            + std::abs(sten(ii-1,jj-1,kk-1,ist_p00)) /
            ( std::abs(sten(ii-1,jj-2,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_0p0)) /
            ( std::abs(sten(ii-2,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_00p)) /
            ( std::abs(sten(ii-2,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_pp0)) /
            ( std::abs(sten(ii-1,jj-1,kk-2,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_p0p)) /
            ( std::abs(sten(ii-1,jj-2,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
            + std::abs(sten(ii-1,jj-1,kk-1,ist_0pp)) /
            ( std::abs(sten(ii-2,jj-1,kk-1,ist_ppp))
            + std::abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps);
        wppp *= std::abs(sten(ii-1,jj-1,kk-1,ist_ppp));
        cv += wppp*fine(ii-1,jj-1,kk-1)*sten(ii-1,jj-1,kk-1,ist_inv);

        crse(i,j,k) = cv * 0.125;
    }
}

#ifdef AMREX_USE_EB

namespace {

    constexpr int i_S_x     = 0;
    constexpr int i_S_y     = 1;
    constexpr int i_S_z     = 2;
    constexpr int i_S_x2    = 3;
    constexpr int i_S_y2    = 4;
    constexpr int i_S_z2    = 5;
    constexpr int i_S_x_y   = 6;
    constexpr int i_S_x_z   = 7;
    constexpr int i_S_y_z   = 8;
    constexpr int i_S_x2_y  = 9;
    constexpr int i_S_x2_z  = 10;
    constexpr int i_S_x_y2  = 11;
    constexpr int i_S_y2_z  = 12;
    constexpr int i_S_x_z2  = 13;
    constexpr int i_S_y_z2  = 14;
    constexpr int i_S_x2_y2 = 15;
    constexpr int i_S_x2_z2 = 16;
    constexpr int i_S_y2_z2 = 17;
    constexpr int i_S_xyz   = 18;
    constexpr int n_Sintg   = 19;

    constexpr int i_c_xmym = 0;
    constexpr int i_c_xmyb = 1;
    constexpr int i_c_xmyp = 2;
    constexpr int i_c_xbym = 3;
    constexpr int i_c_xbyb = 4;
    constexpr int i_c_xbyp = 5;
    constexpr int i_c_xpym = 6;
    constexpr int i_c_xpyb = 7;
    constexpr int i_c_xpyp = 8;
    constexpr int i_c_xmzm = 9;
    constexpr int i_c_xmzb = 10;
    constexpr int i_c_xmzp = 11;
    constexpr int i_c_xbzm = 12;
    constexpr int i_c_xbzb = 13;
    constexpr int i_c_xbzp = 14;
    constexpr int i_c_xpzm = 15;
    constexpr int i_c_xpzb = 16;
    constexpr int i_c_xpzp = 17;
    constexpr int i_c_ymzm = 18;
    constexpr int i_c_ymzb = 19;
    constexpr int i_c_ymzp = 20;
    constexpr int i_c_ybzm = 21;
    constexpr int i_c_ybzb = 22;
    constexpr int i_c_ybzp = 23;
    constexpr int i_c_ypzm = 24;
    constexpr int i_c_ypzb = 25;
    constexpr int i_c_ypzp = 26;
    constexpr int n_conn = 27;

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_set_connection (int i, int j, int k, Array4<Real> const& conn,
                             Array4<Real const> const& intg, Array4<Real const> const& vol,
                             Array4<EBCellFlag const> const& flag) noexcept
{
    if (flag(i,j,k).isCovered()) {
        for (int n = 0; n < n_conn; ++n) conn(i,j,k,n) = 0._rt;
    } else if (flag(i,j,k).isRegular() or vol(i,j,k) >= almostone) {
        for (int n = 0; n < n_conn; ++n) conn(i,j,k,n) = 1._rt;
    } else {
        // Scaled by 9
        conn(i,j,k,i_c_xmym) = 0.5625_rt*vol(i,j,k)
            + 2.25_rt*(-intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_y)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + 9._rt*( intg(i,j,k,i_S_x_y ) - intg(i,j,k,i_S_x2_y)
                     -intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 18
        conn(i,j,k,i_c_xmyb) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*(-intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2))
            + 18._rt*( intg(i,j,k,i_S_x_y2) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xmyp) =  0.5625_rt*vol(i,j,k)
            + 2.25_rt*(-intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_y)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + 9._rt*(-intg(i,j,k,i_S_x_y ) + intg(i,j,k,i_S_x2_y)
                     -intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 18
        conn(i,j,k,i_c_xbym) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*(-intg(i,j,k,i_S_y) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + 18._rt*(intg(i,j,k,i_S_x2_y) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 36
        conn(i,j,k,i_c_xbyb) = 2.25_rt*vol(i,j,k)
            + 9._rt*(-intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2))
            + 36._rt*intg(i,j,k,i_S_x2_y2);

        // Scaled by 18
        conn(i,j,k,i_c_xbyp) =  1.125_rt*vol(i,j,k)
            + 4.5_rt*( intg(i,j,k,i_S_y) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + 18._rt*(-intg(i,j,k,i_S_x2_y) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xpym) = 0.5625_rt*vol(i,j,k)
            + 2.25_rt*( intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_y)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + 9._rt*(-intg(i,j,k,i_S_x_y ) - intg(i,j,k,i_S_x2_y)
                     +intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 18
        conn(i,j,k,i_c_xpyb) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*( intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2))
            + 18._rt*(-intg(i,j,k,i_S_x_y2) - intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xpyp) =  0.5625_rt*vol(i,j,k)
            + 2.25_rt*( intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_y)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2))
            + 9._rt*( intg(i,j,k,i_S_x_y ) + intg(i,j,k,i_S_x2_y)
                      +intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2));

        // Scaled by 9
        conn(i,j,k,i_c_xmzm) = 0.5625_rt*vol(i,j,k)
            + 2.25_rt*(-intg(i,j,k,i_S_x) - intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + 9._rt*(intg(i,j,k,i_S_x_z) - intg(i,j,k,i_S_x2_z)
                     -intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xmzb) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*(-intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_z2))
            + 18._rt*(intg(i,j,k,i_S_x_z2) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_xmzp) = 0.5625_rt*vol(i,j,k)
            + 2.25_rt*(-intg(i,j,k,i_S_x  ) + intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + 9._rt*(-intg(i,j,k,i_S_x_z  ) + intg(i,j,k,i_S_x2_z)
                     -intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xbzm) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*(-intg(i,j,k,i_S_z) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + 18._rt*(intg(i,j,k,i_S_x2_z) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xbzb) = 2.25_rt*vol(i,j,k)
            + 9._rt*(-intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_z2))
            + 36._rt*intg(i,j,k,i_S_x2_z2);

        // Scaled by 18
        conn(i,j,k,i_c_xbzp) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*( intg(i,j,k,i_S_z) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + 18._rt*(-intg(i,j,k,i_S_x2_z) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_xpzm) = 0.5625_rt*vol(i,j,k)
            + 2.25_rt*( intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + 9._rt*(-intg(i,j,k,i_S_x_z ) - intg(i,j,k,i_S_x2_z)
                     +intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_xpzb) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*( intg(i,j,k,i_S_x   ) + intg(i,j,k,i_S_x2   ) - intg(i,j,k,i_S_z2))
            + 18._rt*(-intg(i,j,k,i_S_x_z2) - intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_xpzp) = 0.5625_rt*vol(i,j,k)
            + 2.25_rt*( intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2))
            + 9._rt*( intg(i,j,k,i_S_x_z ) + intg(i,j,k,i_S_x2_z)
                      +intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ymzm) = 0.5625_rt*vol(i,j,k)
            + 2.25_rt*(-intg(i,j,k,i_S_y) - intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + 9._rt*(intg(i,j,k,i_S_y_z) - intg(i,j,k,i_S_y2_z)
                     -intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_ymzb) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*(-intg(i,j,k,i_S_y) + intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2))
            + 18._rt*(intg(i,j,k,i_S_y_z2) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ymzp) = 0.5625_rt*vol(i,j,k)
            + 2.25_rt*(-intg(i,j,k,i_S_y ) + intg(i,j,k,i_S_z)
                       +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + 9._rt*(-intg(i,j,k,i_S_y_z ) + intg(i,j,k,i_S_y2_z)
                     -intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_ybzm) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*(-intg(i,j,k,i_S_z) - intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + 18._rt*(intg(i,j,k,i_S_y2_z) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 36
        conn(i,j,k,i_c_ybzb) = 2.25_rt*vol(i,j,k)
            + 9._rt*(-intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2))
                     + 36._rt*intg(i,j,k,i_S_y2_z2);

        // Scaled by 18
        conn(i,j,k,i_c_ybzp) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*( intg(i,j,k,i_S_z) - intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + 18._rt*(-intg(i,j,k,i_S_y2_z) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ypzm) = 0.5625_rt*vol(i,j,k)
            + 2.25_rt*( intg(i,j,k,i_S_y ) - intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + 9._rt*(-intg(i,j,k,i_S_y_z ) - intg(i,j,k,i_S_y2_z)
                     +intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));

        // Scaled by 18
        conn(i,j,k,i_c_ypzb) = 1.125_rt*vol(i,j,k)
            + 4.5_rt*( intg(i,j,k,i_S_y   ) + intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2))
            + 18._rt*(-intg(i,j,k,i_S_y_z2) - intg(i,j,k,i_S_y2_z2));

        // Scaled by 9
        conn(i,j,k,i_c_ypzp) = 0.5625_rt*vol(i,j,k)
            + 2.25_rt*( intg(i,j,k,i_S_y ) + intg(i,j,k,i_S_z)
                        +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2))
            + 9._rt*( intg(i,j,k,i_S_y_z ) + intg(i,j,k,i_S_y2_z)
                      +intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_set_stencil_eb (int i, int j, int k, Array4<Real> const& sten,
                             Array4<Real const> const& sig, Array4<Real const> const& conn,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1._rt/36._rt)*dxinv[0]*dxinv[0];
    Real facy = (1._rt/36._rt)*dxinv[1]*dxinv[1];
    Real facz = (1._rt/36._rt)*dxinv[2]*dxinv[2];

    // i+1,j,k
    sten(i,j,k,ist_p00) = (
        sig(i,j  ,k  )*(4._rt*facx*conn(i,j  ,k  ,i_c_ymzm) - 2._rt*facy*conn(i,j  ,k  ,i_c_xbzm) - 2._rt*facz*conn(i,j  ,k  ,i_c_xbym) ) +
        sig(i,j-1,k  )*(4._rt*facx*conn(i,j-1,k  ,i_c_ypzm) - 2._rt*facy*conn(i,j-1,k  ,i_c_xbzm) - 2._rt*facz*conn(i,j-1,k  ,i_c_xbyp) ) +
        sig(i,j  ,k-1)*(4._rt*facx*conn(i,j  ,k-1,i_c_ymzp) - 2._rt*facy*conn(i,j  ,k-1,i_c_xbzp) - 2._rt*facz*conn(i,j  ,k-1,i_c_xbym) ) +
        sig(i,j-1,k-1)*(4._rt*facx*conn(i,j-1,k-1,i_c_ypzp) - 2._rt*facy*conn(i,j-1,k-1,i_c_xbzp) - 2._rt*facz*conn(i,j-1,k-1,i_c_xbyp) ) );

    // i,j+1,k
    sten(i,j,k,ist_0p0) = (
        sig(i  ,j,k  )*(-2._rt*facx*conn(i  ,j,k  ,i_c_ybzm) + 4._rt*facy*conn(i  ,j,k  ,i_c_xmzm) - 2._rt*facz*conn(i  ,j,k  ,i_c_xmyb) ) +
        sig(i-1,j,k  )*(-2._rt*facx*conn(i-1,j,k  ,i_c_ybzm) + 4._rt*facy*conn(i-1,j,k  ,i_c_xpzm) - 2._rt*facz*conn(i-1,j,k  ,i_c_xpyb) ) +
        sig(i  ,j,k-1)*(-2._rt*facx*conn(i  ,j,k-1,i_c_ybzp) + 4._rt*facy*conn(i  ,j,k-1,i_c_xmzp) - 2._rt*facz*conn(i  ,j,k-1,i_c_xmyb) ) +
        sig(i-1,j,k-1)*(-2._rt*facx*conn(i-1,j,k-1,i_c_ybzp) + 4._rt*facy*conn(i-1,j,k-1,i_c_xpzp) - 2._rt*facz*conn(i-1,j,k-1,i_c_xpyb) ) );

    // i,j,k+1
    sten(i,j,k,ist_00p) = (
        sig(i  ,j  ,k)*(-2._rt*facx*conn(i  ,j  ,k,i_c_ymzb) - 2._rt*facy*conn(i  ,j  ,k,i_c_xmzb) + 4._rt*facz*conn(i  ,j  ,k,i_c_xmym) ) +
        sig(i-1,j  ,k)*(-2._rt*facx*conn(i-1,j  ,k,i_c_ymzb) - 2._rt*facy*conn(i-1,j  ,k,i_c_xpzb) + 4._rt*facz*conn(i-1,j  ,k,i_c_xpym) ) +
        sig(i  ,j-1,k)*(-2._rt*facx*conn(i  ,j-1,k,i_c_ypzb) - 2._rt*facy*conn(i  ,j-1,k,i_c_xmzb) + 4._rt*facz*conn(i  ,j-1,k,i_c_xmyp) ) +
        sig(i-1,j-1,k)*(-2._rt*facx*conn(i-1,j-1,k,i_c_ypzb) - 2._rt*facy*conn(i-1,j-1,k,i_c_xpzb) + 4._rt*facz*conn(i-1,j-1,k,i_c_xpyp) ) );

    // i+1,j+1,k
    sten(i,j,k,ist_pp0) = (
        sig(i,j,k  )*(2._rt*facx*conn(i,j,k  ,i_c_ybzm) + 2._rt*facy*conn(i,j,k  ,i_c_xbzm) - facz*conn(i,j,k  ,i_c_xbyb) ) +
        sig(i,j,k-1)*(2._rt*facx*conn(i,j,k-1,i_c_ybzp) + 2._rt*facy*conn(i,j,k-1,i_c_xbzp) - facz*conn(i,j,k-1,i_c_xbyb) ) );

    // i+1,j,k+1
    sten(i,j,k,ist_p0p) = (
        sig(i,j,k  )*(2._rt*facx*conn(i,j,k  ,i_c_ymzb) - facy*conn(i,j,k  ,i_c_xbzb) + 2._rt*facz*conn(i,j,k  ,i_c_xbym) ) +
        sig(i,j-1,k)*(2._rt*facx*conn(i,j-1,k,i_c_ypzb) - facy*conn(i,j-1,k,i_c_xbzb) + 2._rt*facz*conn(i,j-1,k,i_c_xbyp) ) );

    // i,j+1,k+1
    sten(i,j,k,ist_0pp) = (
        sig(i  ,j,k)*(-facx*conn(i  ,j,k,i_c_ybzb) + 2._rt*facy*conn(i  ,j,k,i_c_xmzb) + 2._rt*facz*conn(i  ,j,k,i_c_xmyb) ) +
        sig(i-1,j,k)*(-facx*conn(i-1,j,k,i_c_ybzb) + 2._rt*facy*conn(i-1,j,k,i_c_xpzb) + 2._rt*facz*conn(i-1,j,k,i_c_xpyb) ) );

    // i+1,j+1,k+1
    sten(i,j,k,ist_ppp) = sig(i,j,k) * (facx*conn(i,j,k,i_c_ybzb) + facy*conn(i,j,k,i_c_xbzb) + facz*conn(i,j,k,i_c_xbyb) );
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_divu_eb (int i, int j, int k, Array4<Real> const& rhs, Array4<Real const> const& vel,
                      Array4<Real const> const& vfrac, Array4<Real const> const& intg,
                      Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = 0.25_rt*dxinv[0];
    Real facy = 0.25_rt*dxinv[1];
    Real facz = 0.25_rt*dxinv[2];

    if (!msk(i,j,k)) {
        rhs(i,j,k) = facx*(
            vel(i-1,j-1,k  ,0)*(    -vfrac(i-1,j-1,k  )
                               -2._rt*intg(i-1,j-1,k  ,i_S_y)
                               +2._rt*intg(i-1,j-1,k  ,i_S_z)
                               +4._rt*intg(i-1,j-1,k  ,i_S_y_z))
           +vel(i  ,j-1,k  ,0)*(     vfrac(i  ,j-1,k  )
                               +2._rt*intg(i  ,j-1,k  ,i_S_y)
                               -2._rt*intg(i  ,j-1,k  ,i_S_z)
                               -4._rt*intg(i  ,j-1,k  ,i_S_y_z))
           +vel(i-1,j  ,k  ,0)*(    -vfrac(i-1,j  ,k  )
                               +2._rt*intg(i-1,j  ,k  ,i_S_y)
                               +2._rt*intg(i-1,j  ,k  ,i_S_z)
                               -4._rt*intg(i-1,j  ,k  ,i_S_y_z))
           +vel(i  ,j  ,k  ,0)*(     vfrac(i  ,j  ,k  )
                               -2._rt*intg(i  ,j  ,k  ,i_S_y)
                               -2._rt*intg(i  ,j  ,k  ,i_S_z)
                               +4._rt*intg(i  ,j  ,k  ,i_S_y_z))
           +vel(i-1,j-1,k-1,0)*(    -vfrac(i-1,j-1,k-1)
                               -2._rt*intg(i-1,j-1,k-1,i_S_y)
                               -2._rt*intg(i-1,j-1,k-1,i_S_z)
                               -4._rt*intg(i-1,j-1,k-1,i_S_y_z))
           +vel(i  ,j-1,k-1,0)*(     vfrac(i  ,j-1,k-1)
                               +2._rt*intg(i  ,j-1,k-1,i_S_y)
                               +2._rt*intg(i  ,j-1,k-1,i_S_z)
                               +4._rt*intg(i  ,j-1,k-1,i_S_y_z))
           +vel(i-1,j  ,k-1,0)*(    -vfrac(i-1,j  ,k-1)
                               +2._rt*intg(i-1,j  ,k-1,i_S_y)
                               -2._rt*intg(i-1,j  ,k-1,i_S_z)
                               +4._rt*intg(i-1,j  ,k-1,i_S_y_z))
           +vel(i  ,j  ,k-1,0)*(     vfrac(i  ,j  ,k-1)
                               -2._rt*intg(i  ,j  ,k-1,i_S_y)
                               +2._rt*intg(i  ,j  ,k-1,i_S_z)
                               -4._rt*intg(i  ,j  ,k-1,i_S_y_z)) )
            + facy*(
            vel(i-1,j-1,k  ,1)*(    -vfrac(i-1,j-1,k  )
                               -2._rt*intg(i-1,j-1,k  ,i_S_x)
                               +2._rt*intg(i-1,j-1,k  ,i_S_z)
                               +4._rt*intg(i-1,j-1,k  ,i_S_x_z))
           +vel(i  ,j-1,k  ,1)*(    -vfrac(i  ,j-1,k  )
                               +2._rt*intg(i  ,j-1,k  ,i_S_x)
                               +2._rt*intg(i  ,j-1,k  ,i_S_z)
                               -4._rt*intg(i  ,j-1,k  ,i_S_x_z))
           +vel(i-1,j  ,k  ,1)*(     vfrac(i-1,j  ,k  )
                               +2._rt*intg(i-1,j  ,k  ,i_S_x)
                               -2._rt*intg(i-1,j  ,k  ,i_S_z)
                               -4._rt*intg(i-1,j  ,k  ,i_S_x_z))
           +vel(i  ,j  ,k  ,1)*(     vfrac(i  ,j  ,k  )
                               -2._rt*intg(i  ,j  ,k  ,i_S_x)
                               -2._rt*intg(i  ,j  ,k  ,i_S_z)
                               +4._rt*intg(i  ,j  ,k  ,i_S_x_z))
           +vel(i-1,j-1,k-1,1)*(    -vfrac(i-1,j-1,k-1)
                               -2._rt*intg(i-1,j-1,k-1,i_S_x)
                               -2._rt*intg(i-1,j-1,k-1,i_S_z)
                               -4._rt*intg(i-1,j-1,k-1,i_S_x_z))
           +vel(i  ,j-1,k-1,1)*(    -vfrac(i  ,j-1,k-1)
                               +2._rt*intg(i  ,j-1,k-1,i_S_x)
                               -2._rt*intg(i  ,j-1,k-1,i_S_z)
                               +4._rt*intg(i  ,j-1,k-1,i_S_x_z))
           +vel(i-1,j  ,k-1,1)*(     vfrac(i-1,j  ,k-1)
                               +2._rt*intg(i-1,j  ,k-1,i_S_x)
                               +2._rt*intg(i-1,j  ,k-1,i_S_z)
                               +4._rt*intg(i-1,j  ,k-1,i_S_x_z))
           +vel(i  ,j  ,k-1,1)*(     vfrac(i  ,j  ,k-1)
                               -2._rt*intg(i  ,j  ,k-1,i_S_x)
                               +2._rt*intg(i  ,j  ,k-1,i_S_z)
                               -4._rt*intg(i  ,j  ,k-1,i_S_x_z)) )
            + facz*(
            vel(i-1,j-1,k  ,2)*(     vfrac(i-1,j-1,k  )
                               +2._rt*intg(i-1,j-1,k  ,i_S_x)
                               +2._rt*intg(i-1,j-1,k  ,i_S_y)
                               +4._rt*intg(i-1,j-1,k  ,i_S_x_y))
           +vel(i  ,j-1,k  ,2)*(     vfrac(i  ,j-1,k  )
                               -2._rt*intg(i  ,j-1,k  ,i_S_x)
                               +2._rt*intg(i  ,j-1,k  ,i_S_y)
                               -4._rt*intg(i  ,j-1,k  ,i_S_x_y))
           +vel(i-1,j  ,k  ,2)*(     vfrac(i-1,j  ,k  )
                               +2._rt*intg(i-1,j  ,k  ,i_S_x)
                               -2._rt*intg(i-1,j  ,k  ,i_S_y)
                               -4._rt*intg(i-1,j  ,k  ,i_S_x_y))
           +vel(i  ,j  ,k  ,2)*(     vfrac(i  ,j  ,k  )
                               -2._rt*intg(i  ,j  ,k  ,i_S_x)
                               -2._rt*intg(i  ,j  ,k  ,i_S_y)
                               +4._rt*intg(i  ,j  ,k  ,i_S_x_y))
           +vel(i-1,j-1,k-1,2)*(    -vfrac(i-1,j-1,k-1)
                               -2._rt*intg(i-1,j-1,k-1,i_S_x)
                               -2._rt*intg(i-1,j-1,k-1,i_S_y)
                               -4._rt*intg(i-1,j-1,k-1,i_S_x_y))
           +vel(i  ,j-1,k-1,2)*(    -vfrac(i  ,j-1,k-1)
                               +2._rt*intg(i  ,j-1,k-1,i_S_x)
                               -2._rt*intg(i  ,j-1,k-1,i_S_y)
                               +4._rt*intg(i  ,j-1,k-1,i_S_x_y))
           +vel(i-1,j  ,k-1,2)*(    -vfrac(i-1,j  ,k-1)
                               -2._rt*intg(i-1,j  ,k-1,i_S_x)
                               +2._rt*intg(i-1,j  ,k-1,i_S_y)
                               +4._rt*intg(i-1,j  ,k-1,i_S_x_y))
           +vel(i  ,j  ,k-1,2)*(    -vfrac(i  ,j  ,k-1)
                               +2._rt*intg(i  ,j  ,k-1,i_S_x)
                               +2._rt*intg(i  ,j  ,k-1,i_S_y)
                               -4._rt*intg(i  ,j  ,k-1,i_S_x_y)) );
    } else {
        rhs(i,j,k) = 0._rt;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndhelm_mknewu_eb (int i, int j, int k, Array4<Real> const& u, Array4<Real const> const& p,
                        Array4<Real const> const& sig, Array4<Real const> const& vfrac,
                        Array4<Real const> const& intg, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (vfrac(i,j,k) == 0._rt) {
        u(i,j,k,0) = u(i,j,k,1) = u(i,j,k,2) = 0._rt;
    } else {
        Real dpdx = 0.25_rt*(-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  )
                             -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1));
        Real dpdy = 0.25_rt*(-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  )
                             -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));
        Real dpdz = 0.25_rt*(-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  )
                             +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1));

        Real dpp_xy = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1)
                      +p(i+1,j+1,k  ) - p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        Real dpp_xz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) + p(i+1,j,k+1) - p(i,j,k+1)
                      -p(i+1,j+1,k  ) + p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        Real dpp_yz = (p(i+1,j+1,k+1) + p(i,j+1,k+1) - p(i+1,j,k+1) - p(i,j,k+1)
                      -p(i+1,j+1,k  ) - p(i,j+1,k  ) + p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k);

        Real dpp_xyz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1)
                       -p(i+1,j+1,k  ) + p(i,j+1,k  ) + p(i+1,j,k  ) - p(i,j,k  ) ) / vfrac(i,j,k);

        u(i,j,k,0) -= sig(i,j,k)*dxinv[0]*(dpdx + .5_rt*intg(i,j,k,i_S_y  )*dpp_xy +
                                                  .5_rt*intg(i,j,k,i_S_z  )*dpp_xz +
                                                        intg(i,j,k,i_S_y_z)*dpp_xyz );
        u(i,j,k,1) -= sig(i,j,k)*dxinv[1]*(dpdy + .5_rt*intg(i,j,k,i_S_x  )*dpp_xy +
                                                  .5_rt*intg(i,j,k,i_S_z  )*dpp_yz +
                                                        intg(i,j,k,i_S_x_z)*dpp_xyz );
        u(i,j,k,2) -= sig(i,j,k)*dxinv[2]*(dpdz + .5_rt*intg(i,j,k,i_S_x  )*dpp_xz +
                                                  .5_rt*intg(i,j,k,i_S_y  )*dpp_yz +
                                                        intg(i,j,k,i_S_x_y)*dpp_xyz );
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndhelm_rhcc_eb (int i, int j, int k, Array4<Real const> const& rhcc,
                      Array4<Real const> const& vfrac, Array4<Real const> const& intg,
                      Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k)) {
        return
                          rhcc(i  ,j  ,k  ) *
            ( 0.125_rt * vfrac(i  ,j  ,k  )
            + 0.25_rt * (-intg(i  ,j  ,k  ,i_S_x)
                         -intg(i  ,j  ,k  ,i_S_y)
                         -intg(i  ,j  ,k  ,i_S_z))
            + 0.5_rt * (  intg(i  ,j  ,k  ,i_S_x_y)
                         +intg(i  ,j  ,k  ,i_S_x_z)
                         +intg(i  ,j  ,k  ,i_S_y_z))
            +          ( -intg(i  ,j  ,k  ,i_S_xyz)))
            //
            +             rhcc(i-1,j  ,k  ) *
            ( 0.125_rt * vfrac(i-1,j  ,k  )
            + 0.25_rt * ( intg(i-1,j  ,k  ,i_S_x)
                         -intg(i-1,j  ,k  ,i_S_y)
                         -intg(i-1,j  ,k  ,i_S_z))
            + 0.5_rt * ( -intg(i-1,j  ,k  ,i_S_x_y)
                         -intg(i-1,j  ,k  ,i_S_x_z)
                         +intg(i-1,j  ,k  ,i_S_y_z))
            +          (  intg(i-1,j  ,k  ,i_S_xyz)))
            //
            +             rhcc(i  ,j-1,k  ) *
            ( 0.125_rt * vfrac(i  ,j-1,k  )
            + 0.25_rt * (-intg(i  ,j-1,k  ,i_S_x)
                         +intg(i  ,j-1,k  ,i_S_y)
                         -intg(i  ,j-1,k  ,i_S_z))
            + 0.5_rt * ( -intg(i  ,j-1,k  ,i_S_x_y)
                         +intg(i  ,j-1,k  ,i_S_x_z)
                         -intg(i  ,j-1,k  ,i_S_y_z))
            +          (  intg(i  ,j-1,k  ,i_S_xyz)))
            //
            +             rhcc(i-1,j-1,k  ) *
            ( 0.125_rt * vfrac(i-1,j-1,k  )
            + 0.25_rt * ( intg(i-1,j-1,k  ,i_S_x)
                         +intg(i-1,j-1,k  ,i_S_y)
                         -intg(i-1,j-1,k  ,i_S_z))
            + 0.5_rt * (  intg(i-1,j-1,k  ,i_S_x_y)
                         -intg(i-1,j-1,k  ,i_S_x_z)
                         -intg(i-1,j-1,k  ,i_S_y_z))
            +          ( -intg(i-1,j-1,k  ,i_S_xyz)))
            //
            +             rhcc(i  ,j  ,k-1) *
            ( 0.125_rt * vfrac(i  ,j  ,k-1)
            + 0.25_rt * (-intg(i  ,j  ,k-1,i_S_x)
                         -intg(i  ,j  ,k-1,i_S_y)
                         +intg(i  ,j  ,k-1,i_S_z))
            + 0.5_rt * (  intg(i  ,j  ,k-1,i_S_x_y)
                         -intg(i  ,j  ,k-1,i_S_x_z)
                         -intg(i  ,j  ,k-1,i_S_y_z))
            +          (  intg(i  ,j  ,k-1,i_S_xyz)))
            //
            +             rhcc(i-1,j  ,k-1) *
            ( 0.125_rt * vfrac(i-1,j  ,k-1)
            + 0.25_rt * ( intg(i-1,j  ,k-1,i_S_x)
                         -intg(i-1,j  ,k-1,i_S_y)
                         +intg(i-1,j  ,k-1,i_S_z))
            + 0.5_rt * ( -intg(i-1,j  ,k-1,i_S_x_y)
                         +intg(i-1,j  ,k-1,i_S_x_z)
                         -intg(i-1,j  ,k-1,i_S_y_z))
            +          ( -intg(i-1,j  ,k-1,i_S_xyz)))
            //
            +             rhcc(i  ,j-1,k-1) *
            ( 0.125_rt * vfrac(i  ,j-1,k-1)
            + 0.25_rt * (-intg(i  ,j-1,k-1,i_S_x)
                         +intg(i  ,j-1,k-1,i_S_y)
                         +intg(i  ,j-1,k-1,i_S_z))
            + 0.5_rt * ( -intg(i  ,j-1,k-1,i_S_x_y)
                         -intg(i  ,j-1,k-1,i_S_x_z)
                         +intg(i  ,j-1,k-1,i_S_y_z))
            +          ( -intg(i  ,j-1,k-1,i_S_xyz)))
            //
            +             rhcc(i-1,j-1,k-1) *
            ( 0.125_rt * vfrac(i-1,j-1,k-1)
            + 0.25_rt * ( intg(i-1,j-1,k-1,i_S_x)
                         +intg(i-1,j-1,k-1,i_S_y)
                         +intg(i-1,j-1,k-1,i_S_z))
            + 0.5_rt * (  intg(i-1,j-1,k-1,i_S_x_y)
                         +intg(i-1,j-1,k-1,i_S_x_z)
                         +intg(i-1,j-1,k-1,i_S_y_z))
             +         (  intg(i-1,j-1,k-1,i_S_xyz)));
    } else {
        return 0._rt;
    }
}

#endif

}
#endif
