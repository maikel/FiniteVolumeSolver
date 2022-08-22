// Copyright (c) 2021 Maikel Nadolski
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

#ifndef FUB_AMREX_AXIAL_HYPERBOLIC_SPLIT_TIME_INTEGRATOR_HPP
#define FUB_AMREX_AXIAL_HYPERBOLIC_SPLIT_TIME_INTEGRATOR_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/State.hpp"

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <memory>

namespace fub::amrex {

template <typename EulerEquation, typename Tag = execution::OpenMpSimdTag>
struct AxialTimeIntegrator {
  static_assert(EulerEquation::Rank() == 1);

  template <typename F>
  AxialTimeIntegrator(EulerEquation eq, std::vector<::amrex::Geometry> geoms,
                      std::shared_ptr<const std::vector<::amrex::MultiFab>> p,
                      F&& A);

  AxialTimeIntegrator(const AxialTimeIntegrator&);
  AxialTimeIntegrator(AxialTimeIntegrator&&) noexcept = default;

  AxialTimeIntegrator& operator=(const AxialTimeIntegrator&);
  AxialTimeIntegrator& operator=(AxialTimeIntegrator&&) noexcept = default;
  // explicit AxialTimeIntegrator(Tag) {}

  void UpdateConservatively(::amrex::MultiFab& dest,
                            const ::amrex::MultiFab& src,
                            const ::amrex::MultiFab& fluxes,
                            const ::amrex::MultiFab& Ax,
                            const ::amrex::MultiFab& pressure,
                            const ::amrex::Geometry& geom, Duration dt,
                            Direction dir);

  void UpdateConservatively(IntegratorContext& context, int level, Duration dt,
                            Direction dir);

private:
  void RecomputeAxIfNeeded(int level);

  EulerEquation equation_;
  std::vector<::amrex::Geometry> geoms_;
  std::shared_ptr<const std::vector<::amrex::MultiFab>> pressure_;
  std::function<double(double)> A_;
  std::vector<::amrex::MultiFab> Ax_;
  IndexMapping<EulerEquation> indices_{equation_};
};

template <typename EulerEquation, typename F>
AxialTimeIntegrator(const EulerEquation&, std::vector<::amrex::Geometry>,
                    std::shared_ptr<const std::vector<::amrex::MultiFab>>, F &&)
    ->AxialTimeIntegrator<EulerEquation, execution::OpenMpSimdTag>;

template <typename EulerEquation, typename Tag>
template <typename F>
AxialTimeIntegrator<EulerEquation, Tag>::AxialTimeIntegrator(
    EulerEquation eq, std::vector<::amrex::Geometry> geoms,
    std::shared_ptr<const std::vector<::amrex::MultiFab>> p, F&& A)
    : equation_(std::move(eq)), geoms_(std::move(geoms)),
      pressure_(std::move(p)), A_(std::move(A)) {
  if (!pressure_) {
    throw std::runtime_error("AxialTimeIntegrator: shared_ptr 'p' is invalid.");
  }
}

template <typename EulerEquation, typename Tag>
AxialTimeIntegrator<EulerEquation, Tag>::AxialTimeIntegrator(
    const AxialTimeIntegrator& other)
    : equation_{other.equation_}, geoms_{other.geoms_},
      pressure_{other.pressure_}, A_{other.A_}, indices_{other.indices_} {}

template <typename EulerEquation, typename Tag>
AxialTimeIntegrator<EulerEquation, Tag>&
AxialTimeIntegrator<EulerEquation, Tag>::
operator=(const AxialTimeIntegrator& other) {
  return ((*this) = AxialTimeIntegrator{other});
}

template <typename EulerEquation, typename Tag>
void AxialTimeIntegrator<EulerEquation, Tag>::RecomputeAxIfNeeded(int level) {
  FUB_ASSERT(pressure_);
  const std::vector<::amrex::MultiFab>& pressure = *pressure_;
  if (pressure.size() != Ax_.size()) {
    Ax_.resize(pressure.size());
  }
  FUB_ASSERT(static_cast<std::size_t>(level) < Ax_.size());
  const std::size_t i = static_cast<std::size_t>(level);
  const ::amrex::MultiFab& p = pressure[i];
  ::amrex::MultiFab& Ai = Ax_[i];
  if (p.boxArray() != Ai.boxArray() ||
      p.DistributionMap() != Ai.DistributionMap() ||
      p.nGrowVect() != Ai.nGrowVect()) {
    Ai.define(p.boxArray(), p.DistributionMap(), 1, p.nGrowVect());
    ForEachFab(Tag(), Ai, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& fab = Ai[mfi];
      ForEachIndex(mfi.growntilebox(), [&](auto... is) {
        ::amrex::IntVect iv{int(is)...};
        double xM[AMREX_SPACEDIM];
        geoms_[i].LoFace(iv, 0, xM);
        fab(iv, 0) = A_(xM[0]);
      });
    });
  }
}

template <typename EulerEquation, typename Tag>
void AxialTimeIntegrator<EulerEquation, Tag>::UpdateConservatively(
    ::amrex::MultiFab& dest, const ::amrex::MultiFab& src,
    const ::amrex::MultiFab& fluxes, const ::amrex::MultiFab& Ax,
    const ::amrex::MultiFab& pressure, const ::amrex::Geometry& geom,
    Duration dt, [[maybe_unused]] Direction dir) {
  FUB_ASSERT(dir == Direction::X);
  const int n_cons = fluxes.nComp();
  const double dx = geom.CellSize(int(dir));
  ForEachFab(Tag(), dest, [&](::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& next = dest[mfi];
    const ::amrex::FArrayBox& prev = src[mfi];
    const ::amrex::FArrayBox& flux = fluxes[mfi];
    const ::amrex::FArrayBox& area = Ax[mfi];
    const ::amrex::FArrayBox& p = pressure[mfi];
    const ::amrex::Box tilebox = mfi.growntilebox();
    const ::amrex::Box all_faces_tilebox =
        ::amrex::surroundingNodes(tilebox, 0);
    const ::amrex::Box all_fluxes_box = flux.box();
    const ::amrex::Box flux_box = all_faces_tilebox & all_fluxes_box;
    const ::amrex::Box cell_box = enclosedCells(flux_box);
    const IndexBox<AMREX_SPACEDIM + 1> cells = Embed<AMREX_SPACEDIM + 1>(
        AsIndexBox<AMREX_SPACEDIM>(cell_box), {0, n_cons});
    auto nextv = MakePatchDataView(next).Subview(cells);
    auto prevv = MakePatchDataView(prev).Subview(cells);
    const auto faces = Embed<AMREX_SPACEDIM + 1>(
        AsIndexBox<AMREX_SPACEDIM>(flux_box), {0, n_cons});
    auto fluxv = MakePatchDataView(flux).Subview(faces);
    auto areav = MakePatchDataView(area, 0, flux_box);
    auto pv = MakePatchDataView(p, 0, flux_box);
    double dt_over_dx = dt.count() / dx;
    ForEachIndex(cell_box, [&](auto... is) {
      Index<AMREX_SPACEDIM> cell{is...};
      Index<AMREX_SPACEDIM> fL{is...};
      Index<AMREX_SPACEDIM> fR = Shift(fL, Direction::X, 1);

      const double A_left = areav(fL);
      const double A_right = areav(fR);
      const double ooA = 2.0 / (A_left + A_right);
      for (int i = 0; i < n_cons; ++i) {
        nextv(cell, i) =
            prevv(cell, i) +
            dt_over_dx * (A_left * fluxv(fL, i) - A_right * fluxv(fR, i)) * ooA;
      }
      const double p_left = pv(fL);
      const double p_right = pv(fR);
      const double p_center = 0.5 * (p_left + p_right);
      nextv(cell, indices_.momentum[0]) +=
          dt_over_dx * (A_right - A_left) * p_center * ooA;
    });
  });
}

template <typename EulerEquation, typename Tag>
void AxialTimeIntegrator<EulerEquation, Tag>::UpdateConservatively(
    IntegratorContext& context, int level, Duration dt, Direction dir) {
  RecomputeAxIfNeeded(level);
  ::amrex::MultiFab& data = context.GetScratch(level);
  const ::amrex::MultiFab& fluxes = context.GetFluxes(level, dir);
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  this->UpdateConservatively(data, data, fluxes, Ax_[level],
                             (*pressure_)[level], geom, dt, dir);
}

} // namespace fub::amrex

#endif
