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

#include "fub/AMReX/cutcell/HGridFluxMethod.hpp"
#include "fub/AMReX/cutcell/TimeIntegrator.hpp"
#include "fub/ext/omp.hpp"

namespace fub::amrex::cutcell {
namespace {
template <typename FluxMethod>
void ComputeRegularFluxes(FluxMethod& fm, ::amrex::MultiFab& fluxes,
                          const ::amrex::MultiFab& cells,
                          const IntegratorContext& context, Duration dt,
                          double dx, Direction dir) {
  ForEachFab(execution::openmp, fluxes, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::Box face_tilebox = mfi.growntilebox();
    const ::amrex::Box cell_validbox = scratch[mfi].box();
    auto [cell_box, face_box] =
        GetCellsAndFacesInStencilRange(cell_validbox, face_tilebox, gcw, dir);
    auto flux =
        MakeView<Conservative<Equation>>(fluxes[mfi], *equation, face_box);
    auto states =
        MakeView<const Complete<Equation>>(src[mfi], *equation, cell_box);
    const ::amrex::FabType type = context.GetFabType(level, mfi);
    if (type == ::amrex::FabType::singlevalued) {
      CutCellData<AMREX_SPACEDIM> geom = hierarchy.GetCutCellData(level, mfi);
      fm.ComputeRegularFluxes(flux, states, geom, dt, dx, dir);
    } else {
      fm.ComputeNumericFluxes(flux, states, dt, dx, dir);
    }
  });
}

template <typename Equation> struct HGridFM {
  HGridFluxMethod<Equation>* method_;

  /// \brief This function performs an operator splitstep in the specified
  /// direction.
  void DoSplitStep(::amrex::MultiFab& dest, const ::amrex::MultiFab& src,
                   IntegratorContext& simulation_data, int level, Duration dt,
                   Direction dir) {
    OmpLocal<Equation> equation(method_->GetEquation());

    ::amrex::MultiFab& fluxes = simulation_data.GetFluxes(level, dir);
    const double dx = simulation_data.GetDx(level, dir);

    ComputeRegularFluxes(fluxes, src, simulation_data, dt, dx, dir);

    ComputeSinglyShieldedFluxes(src, simulation_data, level, dt, dir);

    TimeIntegrator2::UpdateConservatively(dest, src, simulation_data, level, dt,
                                          dir);
  };
}
} // namespace

// This function performs a half time step with a non-conservative first order
// accurate dimensionally split h-grid method. Then it stores the resulting
// tangential mass flows on the embedded boundaries as reference values.
template <typename Equation>
void HGridFluxMethod<Equation>::PreAdvanceLevel(
    IntegratorContext& simulation_data, int level, Duration dt) {
  const ::amrex::MultiFab& scratch = simulation_data.GetScratch(level);
  const std::size_t slevel = static_cast<std::size_t>(level);

  // Allocate space for values on the half time step level
  ::amrex::MultiFab scratch_half_time = AllocateLike(scratch);

  // Godunov Splitting to the half time step level
  const Duration dt_half = 0.5 * dt;
  auto directions = MakeSplitDirections<Rank>(0, {});
  HGridFM hgrid{this};
  hgrid.DoSplitStep(*this, scratch_half_time, src, simulation_data, level,
                    dt_half, dir);
  auto directions_tail = span(directions).subspan<1>();
  for (Direction dir : directions_tail) {
    hgrid.DoSplitStep(*this, scratch_half_time, scratch_half_time,
                      simulation_data, level, dt_half, dir);
  }

  // Store reference values for mass flow
  if (reference_mass_flow_per_level_.size() < slevel) {
    reference_mass_flow_per_level_.resize(slevel);
  }
  ReallocLike(reference_mass_flow_per_level_[slevel], scratch);
  StoreMassflowOnEmbeddedBoundary_(
      method_, reference_mass_flow_per_level_[slevel], scratch_half_time);
}

template <typename Equation>
void HGridFluxMethod<Equation>::ComputeNumericFluxes(InegratorContext& context,
                                                     int level, Duration dt,
                                                     Direction dir) {
  ComputeRegularFluxes_(context, level, dt, dir);
  ComputeCutCellFluxes_()
}

} // namespace fub::amrex::cutcell