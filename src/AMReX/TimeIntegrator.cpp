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

#include "fub/AMReX/TimeIntegrator.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/ForEach.hpp"
#include "fub/HyperbolicPatchIntegrator.hpp"

namespace fub::amrex {
template <typename Tag>
void ForwardIntegrator<Tag>::UpdateConservatively(
    ::amrex::MultiFab& dest, const ::amrex::MultiFab& src,
    const ::amrex::MultiFab& fluxes, const ::amrex::Geometry& geom, Duration dt,
    Direction dir) {
  const int n_cons = fluxes.nComp();
  const double dx = geom.CellSize(int(dir));
  ::amrex::IntVect fgrow = fluxes.nGrowVect();
  ForEachFab(Tag(), dest, [&](::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& next = dest[mfi];
    const ::amrex::FArrayBox& prev = src[mfi];
    const ::amrex::FArrayBox& flux = fluxes[mfi];
    const ::amrex::Box& box = mfi.tilebox();
    const IndexBox<AMREX_SPACEDIM + 1> cells = Embed<AMREX_SPACEDIM + 1>(AsIndexBox<AMREX_SPACEDIM>(::amrex::grow(box, fgrow)), {0, n_cons});
    auto nv = MakePatchDataView(next).Subview(cells);
    auto pv = MakePatchDataView(prev).Subview(cells);
    const auto faces = Grow(cells, dir, {0, 1});
    auto fv = MakePatchDataView(flux).Subview(faces);
    HyperbolicPatchIntegrator<Tag>::UpdateConservatively(nv, pv, fv, dt, dx,
                                                         dir);
  });
}

template <typename Tag>
void ForwardIntegrator<Tag>::UpdateConservatively(IntegratorContext& context,
                                                  int level, Duration dt,
                                                  Direction dir) {
  ::amrex::MultiFab& data = context.GetScratch(level);
  const ::amrex::MultiFab& fluxes = context.GetFluxes(level, dir);
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  this->UpdateConservatively(data, data, fluxes, geom, dt, dir);
}

template struct ForwardIntegrator<execution::SequentialTag>;
template struct ForwardIntegrator<execution::OpenMpTag>;
template struct ForwardIntegrator<execution::SimdTag>;
template struct ForwardIntegrator<execution::OpenMpSimdTag>;
} // namespace fub::amrex
