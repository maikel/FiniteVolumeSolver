// Copyright (c) 2019 Maikel Nadolski
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

#include "fub/equations/CompressibleAdvection.hpp"
#include "fub/AMReX/ForEachFab.hpp"

namespace fub {
namespace {
double LimitSlopes(double qL, double qM, double qR) {
  const double sL = qM - qL;
  const double sR = qR - qM;
  double r = 0.0;
  if (sL * sR > 0.0) {
    r = sL / sR;
  }
  if (r < 0.0) {
    return 0.0;
  } else {
    return 0.5 * std::min(2 * r / (1 + r), 2 / (1 + r)) * (sL + sR);
  }
}
} // namespace

template <int SpaceDimension>
Duration CompressibleAdvectionFluxMethod<SpaceDimension>::ComputeStableDt(
    amrex::IntegratorContext& context, int level, Direction dir) {
  CompressibleAdvection<SpaceDimension> equation{};
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  const double dx = geom.CellSize(int(dir));
  ::amrex::MultiFab& fluxes = context.GetFluxes(level, dir);
  Duration time_point = context.GetTimePoint(level);
  amrex::ForEachFab(fluxes, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = fluxes[mfi];
    const ::amrex::Box box = mfi.growntilebox();
    View<Conservative> flux = amrex::MakeView<Conservative>(fab, equation, box);
    ForEachIndex(Box<0>(flux), [&](auto... is) {
      ::amrex::IntVect face{int(is)...};
      std::array<double, AMREX_SPACEDIM> x{};
      geom.LoFace(face, int(dir), x.data());
      flux.PTdensity(is...) = Pv_function_(x, time_point, dir);
    });
  });

  const ::amrex::MultiFab& scratch = context.GetScratch(level);
  Duration min_dt(std::numeric_limits<double>::max());
  const int stencil = GetStencilWidth();
  amrex::ForEachFab(scratch, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::Box cell_box = mfi.growntilebox();
    const ::amrex::Box face_box =
        amrex::GetFacesInStencilRange(cell_box, stencil, dir);
    View<const Complete> cells =
        amrex::MakeView<const Complete>(scratch[mfi], equation, cell_box);
    View<const Conservative> flux =
        amrex::MakeView<const Conservative>(fluxes[mfi], equation, face_box);
    StridedDataView<const double, SpaceDimension> Pv = flux.PTdensity;
    Duration local_dt = ComputeStableDt(cells, Pv, dx, dir);
    min_dt = std::min(local_dt, min_dt);
  });
  return min_dt;
}

template <int SpaceDimension>
Duration CompressibleAdvectionFluxMethod<SpaceDimension>::ComputeStableDt(
    const View<const Complete>& states,
    const StridedDataView<const double, SpaceDimension> Pv, double dx,
    Direction dir) {
  double max_signal = std::numeric_limits<double>::lowest();
  ForEachIndex(Box<0>(states), [&](auto... is) {
    using Index = std::array<std::ptrdiff_t, SpaceDimension>;
    Index cell{is...};
    Index faceL = cell;
    Index faceR = Shift(cell, dir, 1);
    double v_advect = 0.5 * (Pv(faceL) + Pv(faceR)) / states.PTdensity(cell);
    max_signal = std::max(max_signal, std::abs(v_advect));
  });
  // FUB_ASSERT(max_signal > 0.0);
  return max_signal > 0 ? Duration(dx / max_signal)
                        : Duration(std::numeric_limits<double>::max());
}

template <int SpaceDimension>
typename CompressibleAdvection<SpaceDimension>::Conservative
CompressibleAdvectionFluxMethod<SpaceDimension>::ComputeNumericFluxes(
    const std::array<Complete, 4>& stencil, const std::array<double, 5> Pvs,
    Duration dt, double dx, Direction) {
  // Reconstruction
  double slope_chi_L = LimitSlopes(stencil[0].PTinverse, stencil[1].PTinverse,
                                   stencil[2].PTinverse) /
                       dx;
  double slope_chi_R = LimitSlopes(stencil[1].PTinverse, stencil[2].PTinverse,
                                   stencil[3].PTinverse) /
                       dx;

  std::array<double, SpaceDimension> slope_velocity_L{};
  std::array<double, SpaceDimension> slope_velocity_R{};

  for (int dim = 0; dim < SpaceDimension; ++dim) {
    slope_velocity_L[dim] =
        LimitSlopes(stencil[0].velocity[dim], stencil[1].velocity[dim],
                    stencil[2].velocity[dim]) /
        dx;
    slope_velocity_R[dim] =
        LimitSlopes(stencil[1].velocity[dim], stencil[2].velocity[dim],
                    stencil[3].velocity[dim]) /
        dx;
  }

  std::array<double, 4> v_advect{};
  // for (int i = 0; i < 4; ++i) {
  // We only use the inner two v_advects
  for (int i = 1; i < 2; ++i) {
    v_advect[i] = 0.5 * (Pvs[i] + Pvs[i + 1]) / stencil[i].PTdensity;
  }

  const double lambda = dt.count() / dx;

  double rec_chi_L = stencil[1].PTinverse +
                     0.5 * dx * slope_chi_L * (1.0 - v_advect[1] * lambda);
  double rec_chi_R = stencil[2].PTinverse -
                     0.5 * dx * slope_chi_R * (1.0 + v_advect[2] * lambda);

  std::array<double, SpaceDimension> rec_velocity_L{};
  std::array<double, SpaceDimension> rec_velocity_R{};

  for (int dim = 0; dim < SpaceDimension; ++dim) {
    rec_velocity_L[dim] =
        stencil[1].velocity[dim] +
        0.5 * dx * slope_velocity_L[dim] * (1.0 - v_advect[1] * lambda);
    rec_velocity_R[dim] =
        stencil[2].velocity[dim] -
        0.5 * dx * slope_velocity_R[dim] * (1.0 + v_advect[2] * lambda);
  }

  int upwind = (Pvs[2] > 0.0) - (Pvs[2] < 0.0);

  Conservative flux{};
  flux.density =
      Pvs[2] * 0.5 * ((1 + upwind) * rec_chi_L + (1 - upwind) * rec_chi_R);
  for (int dim = 0; dim < SpaceDimension; ++dim) {
    flux.momentum[dim] = Pvs[2] * 0.5 *
                         ((1.0 + upwind) * rec_velocity_L[dim] +
                          (1.0 - upwind) * rec_velocity_R[dim]);
  }
  flux.PTdensity = Pvs[2];
  return flux;
}

template <int SpaceDimension>
void CompressibleAdvectionFluxMethod<SpaceDimension>::ComputeNumericFluxes(
    const View<Conservative>& fluxes, const View<const Complete>& states,
    const StridedDataView<const double, SpaceDimension>& Pv_array, Duration dt,
    double dx, Direction dir) {
  std::array<Complete, 4> stencil{};
  std::array<double, 5> Pvs{};

  ForEachIndex(Box<0>(fluxes), [&](auto... is) {
    using Index = std::array<std::ptrdiff_t, SpaceDimension>;
    Index face{is...};
    Index cell_LL = Shift(face, dir, -2);
    Index cell_L = Shift(face, dir, -1);
    Index cell_R = Shift(face, dir, 0);
    Index cell_RR = Shift(face, dir, 1);
    Load(stencil[0], states, cell_LL);
    Load(stencil[1], states, cell_L);
    Load(stencil[2], states, cell_R);
    Load(stencil[3], states, cell_RR);
    // Pvs[0] = Pv_array(Shift(face, dir, -2));
    Pvs[1] = Pv_array(Shift(face, dir, -1));
    Pvs[2] = Pv_array(Shift(face, dir, +0));
    Pvs[3] = Pv_array(Shift(face, dir, +1));
    // Pvs[4] = Pv_array(Shift(face, dir, +2));
    Conservative flux = ComputeNumericFluxes(stencil, Pvs, dt, dx, dir);
    Store(fluxes, flux, face);
  });
}

template <int SpaceDimension>
void CompressibleAdvectionFluxMethod<SpaceDimension>::ComputeNumericFluxes(
    const View<Conservative>& fluxes, const View<const Complete>& states,
    Duration dt, double dx, Direction dir) {
  ComputeNumericFluxes(fluxes, states, fluxes.PTdensity, dt, dx, dir);
}

template struct CompressibleAdvectionFluxMethod<2>;

} // namespace fub