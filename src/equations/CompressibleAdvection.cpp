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
#include "fub/AMReX/CompressibleAdvectionIntegratorContext.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

namespace fub {
namespace {
double LimitSlopes(double qL, double qM, double qR) {
  const double sL = qM - qL;
  const double sR = qR - qM;
  // double r = 0.0;
  return 0.5 * (sL + sR);
  //   if (sL * sR > 0.0) {
  //     r = sL / sR;
  //   }
  //   if (r < 0.0) {
  //     return 0.0;
  //   } else {
  //     return 0.5 * std::min(2 * r / (1 + r), 2 / (1 + r)) * (sL + sR);
  //   }
}
} // namespace

template <int SpaceDimension, int VelocityDimension>
Duration CompressibleAdvectionFluxMethod<SpaceDimension, VelocityDimension>::
    ComputeStableDt(amrex::IntegratorContext& base_context, int level,
                    Direction dir) {
  // try to convert the IntegratorContext to
  // CompressibleAdvectionIntegratorContext, else terminate the program.
  amrex::CompressibleAdvectionIntegratorContext* pointer_to_context =
      dynamic_cast<amrex::CompressibleAdvectionIntegratorContext*>(
          &base_context);
  FUB_ASSERT(pointer_to_context);
  amrex::CompressibleAdvectionIntegratorContext& context = *pointer_to_context;
  CompressibleAdvection<SpaceDimension> equation{};
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  const std::size_t dir_v = static_cast<std::size_t>(dir);
  const double dx = geom.CellSize(static_cast<int>(dir));
  const ::amrex::MultiFab& scratch = context.GetScratch(level);
  const amrex::CompressibleAdvectionAdvectiveFluxes& Pvs =
      context.GetAdvectiveFluxes(level);
  Duration min_dt(std::numeric_limits<double>::max());
  // const int stencil = GetStencilWidth();
  amrex::ForEachFab(
      execution::openmp, Pvs.on_faces[dir_v], [&](const ::amrex::MFIter& mfi) {
        const ::amrex::Box face_tilebox = mfi.growntilebox();
        ::amrex::Box cell_tilebox = ::amrex::enclosedCells(face_tilebox);
        View<const Complete> cells = amrex::MakeView<const Complete>(
            scratch[mfi], equation, cell_tilebox);
        StridedDataView<const double, SpaceDimension> Pv;
        if constexpr (SpaceDimension == AMREX_SPACEDIM) {
          Pv = amrex::MakePatchDataView(Pvs.on_faces[dir_v][mfi], 0)
                   .Subview(amrex::AsIndexBox<SpaceDimension>(face_tilebox));
        } else if constexpr (SpaceDimension + 1 == AMREX_SPACEDIM) {
          Pv = SliceLast(
              amrex::MakePatchDataView(Pvs.on_faces[dir_v][mfi], 0)
                  .Subview(amrex::AsIndexBox<AMREX_SPACEDIM>(face_tilebox)),
              0);
        }
        Duration local_dt = ComputeStableDt(cells, Pv, dx, dir);
        min_dt = std::min(local_dt, min_dt);
      });
  return min_dt;
}

template <int SpaceDimension, int VelocityDimension>
Duration CompressibleAdvectionFluxMethod<SpaceDimension, VelocityDimension>::
    ComputeStableDt(const View<const Complete>& states,
                    const StridedDataView<const double, SpaceDimension> Pv,
                    double dx, Direction dir) {
  double max_signal = std::numeric_limits<double>::lowest();
  ForEachIndex(Shrink(Pv.Box(), dir, {0, 1}), [&](auto... is) {
    using Index =
        std::array<std::ptrdiff_t, static_cast<std::size_t>(SpaceDimension)>;
    Index faceL{is...};
    Index faceR = Shift(faceL, dir, 1);
    Index cell = faceL;
    FUB_ASSERT(states.PTdensity(cell) > 0.0);
    double v_advect = 0.5 * (Pv(faceL) + Pv(faceR)) / states.PTdensity(cell);
    max_signal = std::max(max_signal, std::abs(v_advect));
  });
  // FUB_ASSERT(max_signal > 0.0);
  return max_signal > 0 ? Duration(dx / max_signal)
                        : Duration(std::numeric_limits<double>::max());
}

template <int SpaceDimension, int VelocityDimension>
typename CompressibleAdvection<SpaceDimension, VelocityDimension>::Conservative
CompressibleAdvectionFluxMethod<SpaceDimension, VelocityDimension>::
    ComputeNumericFluxes(const std::array<Complete, 4>& stencil,
                         const std::array<double, 5> Pvs, Duration dt,
                         double dx, Direction) {
  const std::size_t VelocityDimensionli =
      static_cast<std::size_t>(VelocityDimension);

  // Reconstruction
  double slope_chi_L = LimitSlopes(stencil[0].PTinverse, stencil[1].PTinverse,
                                   stencil[2].PTinverse) /
                       dx;
  double slope_chi_R = LimitSlopes(stencil[1].PTinverse, stencil[2].PTinverse,
                                   stencil[3].PTinverse) /
                       dx;

  std::array<double, VelocityDimension> slope_velocity_L{};
  std::array<double, VelocityDimension> slope_velocity_R{};

  for (int dim = 0; dim < VelocityDimension; ++dim) {
    const std::size_t dimli = static_cast<std::size_t>(dim);
    std::array<double, 4> velocity;
    velocity[0] = stencil[0].momentum[dim] / stencil[0].density;
    velocity[1] = stencil[1].momentum[dim] / stencil[1].density;
    velocity[2] = stencil[2].momentum[dim] / stencil[2].density;
    velocity[3] = stencil[3].momentum[dim] / stencil[3].density;
    slope_velocity_L[dimli] =
        LimitSlopes(velocity[0], velocity[1], velocity[2]) / dx;
    slope_velocity_R[dimli] =
        LimitSlopes(velocity[1], velocity[2], velocity[3]) / dx;
  }

  double slope_chi_fast_L =
      LimitSlopes(stencil[0].chi_fast / stencil[0].density,
                  stencil[1].chi_fast / stencil[1].density,
                  stencil[2].chi_fast / stencil[2].density);
  double slope_chi_fast_R =
      LimitSlopes(stencil[1].chi_fast / stencil[1].density,
                  stencil[2].chi_fast / stencil[2].density,
                  stencil[3].chi_fast / stencil[3].density);

  std::array<double, 4> v_advect{};
  // for (int i = 0; i < 4; ++i) {
  // We only use the inner two v_advects
  for (std::size_t i = 1; i < 3; ++i) {
    v_advect[i] = 0.5 * (Pvs[i] + Pvs[i + 1]) / stencil[i].PTdensity;
  }

  const double lambda = dt.count() / dx;

  double rec_chi_L = stencil[1].PTinverse +
                     0.5 * dx * slope_chi_L * (1.0 - v_advect[1] * lambda);
  double rec_chi_R = stencil[2].PTinverse -
                     0.5 * dx * slope_chi_R * (1.0 + v_advect[2] * lambda);

  std::array<double, VelocityDimensionli> rec_velocity_L{};
  std::array<double, VelocityDimensionli> rec_velocity_R{};

  for (int dim = 0; dim < VelocityDimension; ++dim) {
    const std::size_t dimli = static_cast<std::size_t>(dim);
    rec_velocity_L[dimli] =
        stencil[1].momentum[dim] / stencil[1].density +
        0.5 * dx * slope_velocity_L[dimli] * (1.0 - v_advect[1] * lambda);
    rec_velocity_R[dimli] =
        stencil[2].momentum[dim] / stencil[2].density -
        0.5 * dx * slope_velocity_R[dimli] * (1.0 + v_advect[2] * lambda);
  }

  double rec_chi_fast_L =
      stencil[1].chi_fast / stencil[1].density +
      0.5 * dx * slope_chi_fast_L[dimli] * (1.0 - v_advect[1] * lambda);
  double rec_chi_fast_R =
      stencil[2].chi_fast / stencil[2].density -
      0.5 * dx * slope_velocity_R[dimli] * (1.0 + v_advect[2] * lambda);

  int upwind = (Pvs[2] > 0.0) - (Pvs[2] < 0.0);

  Conservative flux{};
  flux.density =
      Pvs[2] * 0.5 * ((1 + upwind) * rec_chi_L + (1 - upwind) * rec_chi_R);
  for (int dim = 0; dim < VelocityDimension; ++dim) {
    const std::size_t dimli = static_cast<std::size_t>(dim);
    flux.momentum[dim] = Pvs[2] * 0.5 *
                         ((1.0 + upwind) * rec_chi_L * rec_velocity_L[dimli] +
                          (1.0 - upwind) * rec_chi_R * rec_velocity_R[dimli]);
  }
  flux.PTdensity = Pvs[2];
  flux.chi_fast = Pvs[2] * 0.5 *
                  ((1.0 + upwind) * rec_chi_L * rec_chi_fast_L +
                   (1.0 - upwind) * rec_chi_R * rec_chi_fast_R);
  return flux;
}

template <int SpaceDimension, int VelocityDimension>
void CompressibleAdvectionFluxMethod<SpaceDimension, VelocityDimension>::
    ComputeNumericFluxes(
        const View<Conservative>& fluxes, const View<const Complete>& states,
        const StridedDataView<const double, SpaceDimension>& Pv_array,
        Duration dt, double dx, Direction dir) {
  std::array<Complete, 4> stencil{};
  std::array<double, 5> Pvs{};

  ForEachIndex(Box<0>(fluxes), [&](auto... is) {
    using Index =
        std::array<std::ptrdiff_t, static_cast<std::size_t>(SpaceDimension)>;
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

template <int SpaceDimension, int VelocityDimension>
void CompressibleAdvectionFluxMethod<SpaceDimension, VelocityDimension>::
    ComputeNumericFluxes(amrex::IntegratorContext& base_context, int level,
                         Duration dt, Direction dir) {
  // try to convert the IntegratorContext to
  // CompressibleAdvectionIntegratorContext, else terminate the program.
  amrex::CompressibleAdvectionIntegratorContext* pointer_to_context =
      dynamic_cast<amrex::CompressibleAdvectionIntegratorContext*>(
          &base_context);
  FUB_ASSERT(pointer_to_context);
  amrex::CompressibleAdvectionIntegratorContext& context = *pointer_to_context;
  CompressibleAdvection<SpaceDimension> equation{};
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  const auto dir_v = static_cast<int>(dir);
  const auto dir_s = static_cast<std::size_t>(dir);
  const double dx = geom.CellSize(dir_v);
  const ::amrex::MultiFab& scratch = context.GetScratch(level);
  ::amrex::MultiFab& fluxes = context.GetFluxes(level, dir);
  const amrex::CompressibleAdvectionAdvectiveFluxes& Pvs =
      context.GetAdvectiveFluxes(level);
  const int stencil = GetStencilWidth();
  amrex::ForEachFab(execution::openmp, fluxes, [&](const ::amrex::MFIter& mfi) {
    // Get a view of all complete state variables
    const ::amrex::Box face_tilebox = mfi.growntilebox();
    const ::amrex::Box cell_validbox = scratch[mfi].box();
    const auto [cell_box, face_box] = amrex::GetCellsAndFacesInStencilRange(
        cell_validbox, face_tilebox, stencil, dir);
    const ::amrex::FArrayBox& sfab = scratch[mfi];
    ::amrex::FArrayBox& ffab = fluxes[mfi];
    View<const Complete> cells =
        amrex::MakeView<const Complete>(sfab, equation, cell_box);
    View<Conservative> fluxes =
        amrex::MakeView<Conservative>(ffab, equation, face_box);
    StridedDataView<const double, SpaceDimension> Pv;
    if constexpr (SpaceDimension == AMREX_SPACEDIM) {
      Pv = amrex::MakePatchDataView(Pvs.on_faces[dir_s][mfi], 0)
               .Subview(amrex::AsIndexBox<SpaceDimension>(
                   ::amrex::grow(face_box, dir_v, 1)));
    } else if constexpr (SpaceDimension + 1 == AMREX_SPACEDIM) {
      Pv = SliceLast(amrex::MakePatchDataView(Pvs.on_faces[dir_s][mfi], 0)
                         .Subview(amrex::AsIndexBox<AMREX_SPACEDIM>(
                             ::amrex::grow(face_box, dir_v, 1))),
                     0);
    }
    ComputeNumericFluxes(fluxes, cells, Pv, dt, dx, dir);
  });
}

template struct CompressibleAdvectionFluxMethod<2>;

void Reflect(Complete<CompressibleAdvection<2>>& reflected,
             const Complete<CompressibleAdvection<2>>& state,
             const Eigen::Vector2d& normal, const CompressibleAdvection<2>&) {
  reflected.density = state.density;
  reflected.PTdensity = state.PTdensity;
  reflected.PTinverse = state.PTinverse;
  reflected.momentum =
      state.momentum -
      2 * (state.momentum.matrix().dot(normal) * normal).array();
}

} // namespace fub
